# OpenFAST_to_AcousticAnalogies
This repository links the popular aero-servo-hydro-elastic solver OpenFAST (https://github.com/OpenFAST/openfast) and the acoustic solver AcousticAnalogies.jl (https://github.com/OpenMDAO/AcousticAnalogies.jl)

OpenFAST is released by the U.S. National Renewable Energy Laboratory (NREL) and is a popular open-source wind turbine simulation tool. AcousticAnalogies is released by the OpenMDAO team at NASA.

The example is meant to simulate the rotor of the DOE1.5 wind turbine, which is a GE1.5 SLE wind turbine owned by the U.S. Department of Energy and located at NREL Flatirons Campus in Arvada, Colorado. The proprietary nature of the turbine forced us to mask some of the data such as hub radius and chord distribution.

The code is built in the Julia programming language and heavily borrows from the guided example of AcousticAnalogies.jl, see https://github.com/OpenMDAO/AcousticAnalogies.jl/blob/main/docs/src/guided_example.md