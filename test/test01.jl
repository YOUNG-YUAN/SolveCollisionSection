using SolveCollisionSection, SolveCollisionSection.Units
using Plots

"""
This test.jl is a good example for solve collision section of H-Kr system, which potaintial is 
Lennard-Jones potaintial.
"""

SetUnits(:Atom)  # Set the Units

ρ=3.57Å
α=2*1836*Me / ħ^2

SetPotential(:LennardJonesPotential, 5.9meV, 3.57Å)
M                   = 1836Me    # Mass of Hydrogen atom, where Me is mass of electron
Rmin                = 0.5*ρ
Rmax                = 10*ρ
Step                = 0.001*ρ
L                   = [4,5,6]
Energy              = Vector(0.01meV : 0.01meV : 3.5meV)
iterativeDirection  = :right
iniValue1 = exp(- sqrt(5.9meV * α / 25) * Rmin^(-5))
iniValue2 = exp(- sqrt(5.9meV * α / 25) * (Rmin+Step)^(-5))
SetWaveFunction(M, Rmin, Rmax, Step, iterativeDirection, iniValue1, iniValue2)
(GCS, PS) = GetGrossCrossSection(5.0*ρ, L, Energy, M)
plot(Energy, GCS)

#=
The following code can draw the wave function.
=#

# WFSetting = WaveFunctionSetting(M, Rmin, Rmax, Step, iterativeDirection, iniValue1, iniValue2)
# RWF = SolveWaveFunction(4, 1meV, WFSetting)
# R = Vector{Float64}()
# Wave = Vector{Float64}()
# for i in eachindex(RWF)
#     append!(R, RWF[i].R)
#     append!(Wave, RWF[i].Value)
# end
# plot(R, Wave)

