module GE1p5Downwind

using AcousticAnalogies
using AcousticMetrics
using ColorSchemes: colorschemes
using GLMakie
using DelimitedFiles
using KinematicCoordinateTransformations
using LinearAlgebra: ×
using StaticArrays
using AcousticMetrics
using Statistics


function doit()

    num_blades = 3  # number of blades
    Rhub = 1.25  # meters
    Rtip = 38.5  # meters
    radii = [0, 2.3576, 5.8663, 8.1488, 10.3908, 
        12.5901, 14.7441, 16.85, 18.9047, 20.9043, 
        22.8447, 24.7209, 26.5271, 28.2561, 29.8994, 
        31.4461, 32.8823, 34.1892, 35.3404, 36.2952, 
        36.9837, 37.25	] .+ Rhub

    # Each of the num_blades blades have the same radial coordinates in 
    # the blade-fixed frame, but different angular coordinates: blade number 2
    #  will be offset 180° from blade number 1.
    θs = 2*pi/num_blades.*(0:(num_blades-1))

    # blade element's length from the element centers and the hub and tip location
    # DJI: this `get_dradii` function assumes:
    # DJI;   1. that the blade element interfaces are midway between each radial station defined by the `radii` vector.
    # DJI:   2. that the inner interface of the first element is at `Rhub`
    # DJI:   3. thta the outer interface of the last element is at `Rtip`.
    # DJI: (the `get_dradii` function is very simple, only a few lines: https://github.com/OpenMDAO/AcousticAnalogies.jl/blob/57c199bc35cd0db6d13f17831ec34e5f1d8a9eea/src/utils.jl#L10
    # DJI: But here, it looks like radii[1] == Rhub, which makes the elements near the hub look a little goofy.
    # DJI: using the `get_dradii` routine is optional---it just returns a vector of the radial spacing of each blade element, but if you know that already, you can just use that.
    dradii = get_dradii(radii, Rhub, Rtip)

    # cross-sectional area of each element
    cs_area_over_chord_squared = 0.064

    chord = [2., 2., 3., 3., 3., 3., 
        3., 3., 2., 2., 2., 2., 
        2., 1., 1., 1., 1., 1., 1., 
        1., 0.2, 0.1] # fake_chord


    fig = Figure()
    ax1 = fig[1, 1] = Axis(fig, xlabel="Span Position (m)", ylabel="Chord (m)")
    l1 = lines!(ax1, radii, chord)
    hidexdecorations!(ax1, grid=false)
    save(joinpath(@__DIR__, "chord.png"), fig)



    cs_area = cs_area_over_chord_squared.*chord.^2

    # load openfast output
    # Import the DelimitedFiles module
    # Define the file path
    file_path = joinpath(@__DIR__, "turbulent_720s.out")

    # Function to parse a line of data, converting strings to floats
    function parse_line(line)
        # Split the line by whitespace and filter out any empty strings
        elements = filter(x -> !isempty(x), split(line))
        # Convert elements to Float64
        return map(x -> parse(Float64, x), elements)
    end

    # Initialize an empty array to store the data
    data = []

    # Open the file and read the data, skipping the first 8 lines
    open(file_path) do file
        # Skip the first 8 lines (header and description)
        for i in 1:8
            readline(file)
        end

        # Read the rest of the lines and parse them
        for line in eachline(file)
            push!(data, parse_line(line))
        end
    end

    # Convert the data to an array of arrays (matrix)
    data = reduce(hcat, data)
    time = data[1, :]
    avg_wind_speed = mean(data[2, :])
    sim_length_s = time[end] - time[1] # s
    @show length(time)

    # Reopen the file and read the lines
    lines = open(file_path) do f
        readlines(f)
    end

    # Find the index of the line that contains the column headers
    header_index = findfirst(x -> startswith(x, "Time"), lines)

    # Extract the headers
    headers = split(lines[header_index], '\t')

    id_b1_Fn = findfirst(x -> x == "AB1N001Fn", headers)
    id_b2_Fn = findfirst(x -> x == "AB2N001Fn", headers)
    id_b3_Fn = findfirst(x -> x == "AB3N001Fn", headers)
    id_b1_Ft = findfirst(x -> x == "AB1N001Ft", headers)
    id_b2_Ft = findfirst(x -> x == "AB2N001Ft", headers)
    id_b3_Ft = findfirst(x -> x == "AB3N001Ft", headers)
    id_rot_speed = findfirst(x -> x == "RotSpeed", headers)
    n_elems = length(radii)
    Fn_b1 = data[id_b1_Fn:id_b1_Fn+n_elems-1,:]
    Ft_b1 = data[id_b1_Ft:id_b1_Ft+n_elems-1,:]
    Fn_b2 = data[id_b2_Fn:id_b2_Fn+n_elems-1,:]
    Ft_b2 = data[id_b2_Ft:id_b2_Ft+n_elems-1,:]
    Fn_b3 = data[id_b3_Fn:id_b3_Fn+n_elems-1,:]
    Ft_b3 = data[id_b3_Ft:id_b3_Ft+n_elems-1,:]
    omega_rpm = mean(data[id_rot_speed,:])

    # DJI: are these loadings normal/tangential to the blade's rotation?
    # DJI: or are they aligned with/normal to the chord line?
    # DJI: if they're normal/tangential to the blade's rotation, then I think everything looks good.
    # DJI: if they're aligned with/normal to the chord line, the I think they would need to be rotated by the amount of twist at each radial station.
    # DJI: could do that manually, or use one of the `SteadyRot{X,Y,Z}Transformation` things from KinematicCoordinateTransformations.
    fig = Figure()
    ax1 = fig[1, 1] = Axis(fig, xlabel="Span Position (m)", ylabel="Fn (N/m)")
    ax2 = fig[2, 1] = Axis(fig, xlabel="Span Position (m)", ylabel="Ft (N/m)")
    l1 = lines!(ax1, radii, Fn_b1[:,1], label ="b1")
    l1 = lines!(ax1, radii, Fn_b2[:,1], label ="b2")
    l1 = lines!(ax1, radii, Fn_b3[:,1], label ="b3")
    l2 = lines!(ax2, radii, Ft_b1[:,1], label ="b1")
    l2 = lines!(ax2, radii, Ft_b2[:,1], label ="b2")
    l2 = lines!(ax2, radii, Ft_b3[:,1], label ="b3")
    hidexdecorations!(ax1, grid=false)
    save(joinpath(@__DIR__, "Fn_t.png"), fig)

    # DJI: I was curious how unsteady the loading was, so I plotted each blade's normal and tangential loading as a function of span position for a bunch of different times.
    # DJI: times are indicated by the colorbar on the right of the plot.
    # DJI: only plotting 1 of every 500 timesteps (note the `for tidx in 1:500:ntimes_loading`)
    @assert size(Fn_b1) == size(Fn_b2) == size(Fn_b3) == size(Ft_b1) == size(Ft_b2) == size(Ft_b3) 
    ntimes_loading = size(Fn_b1, 2)
    fig = Figure()
    ax11 = fig[1, 1] = Axis(fig, xlabel="Span Position (m)", ylabel="Fn (N/m)", title="blade 1")
    ax21 = fig[2, 1] = Axis(fig, xlabel="Span Position (m)", ylabel="Ft (N/m)")
    ax12 = fig[1, 2] = Axis(fig, xlabel="Span Position (m)", ylabel="Fn (N/m)", title="blade 2")
    ax22 = fig[2, 2] = Axis(fig, xlabel="Span Position (m)", ylabel="Ft (N/m)")
    ax13 = fig[1, 3] = Axis(fig, xlabel="Span Position (m)", ylabel="Fn (N/m)", title="blade 3")
    ax23 = fig[2, 3] = Axis(fig, xlabel="Span Position (m)", ylabel="Ft (N/m)")
    bpp = 60/omega_rpm/num_blades
    colormap = colorschemes[:viridis]
    for tidx in 1:500:ntimes_loading
        cidx = (time[tidx] - time[1])/sim_length_s
        l1 = lines!(ax11, radii, Fn_b1[:,tidx], label ="b1", color=colormap[cidx])
        l1 = lines!(ax12, radii, Fn_b2[:,tidx], label ="b2", color=colormap[cidx])
        l1 = lines!(ax13, radii, Fn_b3[:,tidx], label ="b3", color=colormap[cidx])
        l2 = lines!(ax21, radii, Ft_b1[:,tidx], label ="b1", color=colormap[cidx])
        l2 = lines!(ax22, radii, Ft_b2[:,tidx], label ="b2", color=colormap[cidx])
        l2 = lines!(ax23, radii, Ft_b3[:,tidx], label ="b3", color=colormap[cidx])
    end

    cbar = Colorbar(fig[:, 4]; colormap=:viridis, limits=(0, (time[end] - time[1])/bpp), label="source time, blade passes")#, ticks=0.0:0.2:1.0)


    linkxaxes!(ax21, ax11)
    linkxaxes!(ax21, ax11)
    linkxaxes!(ax12, ax11)
    linkxaxes!(ax22, ax11)
    linkxaxes!(ax13, ax11)
    linkxaxes!(ax23, ax11)

    linkyaxes!(ax12, ax11)
    linkyaxes!(ax13, ax11)
    linkyaxes!(ax22, ax21)
    linkyaxes!(ax23, ax21)

    hidexdecorations!(ax11, grid=false)
    hidexdecorations!(ax12, grid=false)
    hidexdecorations!(ax13, grid=false)

    hideydecorations!(ax12, grid=false)
    hideydecorations!(ax13, grid=false)
    hideydecorations!(ax22, grid=false)
    hideydecorations!(ax23, grid=false)

    save(joinpath(@__DIR__, "Fn_t-all_time.png"), fig)

    # normal and circumferential loading as a function of radial position along 
    # the blade. The loading is in units of force per unit span (here, Newtons/meter).
    fn = cat(transpose(Fn_b1), transpose(Fn_b2), transpose(Fn_b3), dims=3)
    fc = cat(transpose(Ft_b1), transpose(Ft_b2), transpose(Ft_b3), dims=3)

    # ambient air density and speed of sound.
    rho = 1.0  # kg/m^3
    c0 = 340.0  # m/s

    # rotor motion in lateral speed (0 for wind turbines) and rotational speed in rad/s
    # DJI: freestream velocity is 10 m/s in the positive x direction.
    # DJI: but to do F1A correctly, we need to put all source elements in a coordinate system that moves with the fluid, i.e. one in which the fluid appears stationary.
    # So, to do that, we have the blades translating in the negative x direction.
    v = -avg_wind_speed  # m/s
    omega = omega_rpm * 2*pi/60  # rad/s

    # some reshaping, ses[i, j, k] holds the CompactSourceElement at src_time[i], radii[j], and blade number k
    θs = reshape(θs, 1, 1, :)
    radii = reshape(radii, 1, :, 1)
    dradii = reshape(dradii, 1, :, 1)
    cs_area = reshape(cs_area, 1, :, 1)
    src_times = reshape(time, :, 1, 1)  # This isn't really necessary.

    # source elements, with negative fn
    ses = CompactSourceElement.(rho, c0, radii, θs, dradii, cs_area, -fn, 0.0, fc, src_times)

    t0 = 0.0  # Time at which the angle between the source and target coordinate systems is equal to offest.
    offset = 0.0  # Angular offset between the source and target coordinate systems at t0.
    # steady rotation around the x axis
    rot_trans = SteadyRotXTransformation(t0, omega, offset)

    # orient the rotation axis of the blades as it is the global frame
    rot_axis = @SVector [1.0, 0.0, 0.0] # rotation axis aligned along global x-axis 
    blade_axis = @SVector [0.0, 0.0, 1.0]  # blade 1 pointing up, along z-axis 
    global_trans = ConstantLinearMap(hcat(rot_axis, blade_axis, rot_axis×blade_axis))

    # blade to move with the appropriate forward velocity, and 
    # start from the desired location in the global reference frame
    y0_hub = @SVector [0.0, 0.0, 0.0]  # Position of the hub at time t0
    v0_hub = SVector{3}(v.*rot_axis)   # Constant velocity of the hub in the global reference frame
    const_vel_trans = ConstantVelocityTransformation(t0, y0_hub, v0_hub)

    # combine these three transformations into one, and then use that on the SourceElements
    trans = compose.(src_times, Ref(const_vel_trans), compose.(src_times, Ref(global_trans), Ref(rot_trans)))

    # trans will perform the three transformations from right to left (rot_trans, global_trans, const_vel_trans)
    ses = ses .|> trans

    # ses = AcousticAnalogies.CompactSourceElement.(rho, c0, radii, θs, dradii, cs_area, -fn, 0.0, fc, src_times) .|> trans


    # The ses object now describes how each blade element source is moving through the global reference
    #  frame over the time src_time. As it does this, it will emit acoustics that can be sensed by an acoustic observer 
    # (a human, or a microphone). The exact "amount" of acoustics the observer will experience depends 
    # on the relative location and motion between each source and the observer. 
    # So we'll need to define our acoustic observer before we can calculate the noise heard by it. 
    # For this example, we'll assume that our acoustic observer is stationary in the global frame.

    # DJI: Above we decided we want the blades to translate 10 m/s in the negative x direction to account for the 10 m/s freestream velocity in the positive x direction in the ground-fixed frame.
    # DJI: We need to move the acoustic observer in the same way (with the same velocity.
    x0 = @SVector [118.5, 0.0, -79.0] # IEC location, 118.5 m downwind on the ground (hub height is 79 m)
    # obs = StationaryAcousticObserver(x0)
    # This creates an acoustic observer moving with constant velocity v0_hub that is at location `x0` at time `t0`.
    obs = ConstVelocityAcousticObserver(t0, x0, v0_hub)

    # Now, in order to perform the F1A calculation, 
    # we need to know when each acoustic disturbance emitted 
    # by the source arrives at the observer. This is referred 
    # to an advanced time calculation, and is done this way:
    # That returns an array the same size of ses of the time each acoustic disturbance reaches the observer obs
    # DJI: In the guided example the observer time is calculated explicitly like this, but you can also just go like this:
    # DJI:  apth = f1a.(ses, Ref(obs))
    # DJI: and it will be calculated internally.
    # obs_time = adv_time.(ses, Ref(obs))

    # # compact F1A calculation!
    # apth = f1a.(ses, Ref(obs), obs_time)
    apth = f1a.(ses, Ref(obs))

    # We now have a noise prediction for each of the individual source elements in ses at the acoustic observer obs. 
    # What we ultimately want is the total noise prediction at obs—we want to add all the acoustic pressures in apth together. 
    # But we can't add them directly, yet, since the observer times are not all the same. What we need to do 
    # is first interpolate the apth of each source onto a common observer time grid, and then add them up. 
    # We'll do this using the AcousticAnalogies.combine function.
    period = 2*pi/omega
    bpp = period/num_blades  # blade passing period
    obs_time_range = sim_length_s/60*omega_rpm*bpp
    @show obs_time_range/sim_length_s # should be equalt to the 1/number of blades?
    # DJI: Need to be careful to avoid extrapolation in the `combine` calculation.
    # DJI: That won't happen in this case, since obs_time_range/sim_length_s is 1/3, so the observer time range is much less than the source time range.
    # DJI: The observer time range is 1/3 of the source time range, and we're using the same number of simulation times, so that means the observer time step is 1/3 that of the source time step.
    num_obs_times = length(time)
    apth_total = combine(apth, obs_time_range, num_obs_times, 1)
    # DJI: It appears the loading data is unsteady, so may need to be careful to window the time history to avoid problems with discontinuities going from the begining/end of the pressure time history..

    # We can now have a look at the total acoustic pressure time history at the observer:
    fig = Figure()
    ax1 = fig[1, 1] = Axis(fig, xlabel="time, s", ylabel="monopole, Pa")
    ax2 = fig[2, 1] = Axis(fig, xlabel="time, s", ylabel="dipole, Pa")
    ax3 = fig[3, 1] = Axis(fig, xlabel="time, s", ylabel="total, Pa")
    # DJI: The x axis label says the time is in units of blade passes, but I think the `time` vector here is in seconds.
    # DJI: So need to divide by the blade passing period `bpp` to get blade passes.
    l1 = lines!(ax1, time./bpp, apth_total.p_m)
    l2 = lines!(ax2, time./bpp, apth_total.p_d)
    l3 = lines!(ax3, time./bpp, apth_total.p_m.+apth_total.p_d)
    hidexdecorations!(ax1, grid=false)
    hidexdecorations!(ax2, grid=false)
    save(joinpath(@__DIR__, "openfast-apth_total.png"), fig)
    # DJI: the plot shows that the monopole/thickness noise is much lower than the dipole/loading noise.
    # DJI: I think that makes sense, at least from my experience with propellers/helicopters.
    # DJI: Wind turbine blades are relatively slender, which would tend to reduce thickness noise.
    # DJI: Also the observer is downstream of the rotation plane, which is where loading noise is traditionally thought to dominate (monopole/thickness noise is more significant in the rotor rotation plane, usually).

    # Calculate the overall sound pressure level from the acoustic pressure time history.
    oaspl_from_apth = AcousticMetrics.OASPL(apth_total)
    # Calculate the narrowband spectrum.
    nbs = AcousticMetrics.MSPSpectrumAmplitude(apth_total)
    # Calculate the OASPL from the NBS.
    oaspl_from_nbs = AcousticMetrics.OASPL(nbs)
    (oaspl_from_apth, oaspl_from_nbs)

    name = joinpath(@__DIR__, "vtk", "ge1p5_vtk")
    # DJI: The WriteVTK.jl library assumes that the directory it's writing files to already exists.
    # DJI: So we need to create the directory that will contain the VTK files ahead of time.
    mkpath(dirname(name))
    outfiles = AcousticAnalogies.to_paraview_collection(name, ses)
end

# DJI: When working on a script that I'm running over and over again, I like to put it in a module, include it in a Julia REPL using Revise.jl's `includet` function, then run it from the REPL via GE1p5Downwind.doit().
# DJI: But if you want to manually run the script from a shell, this should do that automatically (sort of like Python's `if __name__ == "__main__":` that you see a the bottom of a lot of Python scripts.
if ! isinteractive()
    doit()
end

end # module
