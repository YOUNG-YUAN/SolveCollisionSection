module SolveCollisionSection

cd("C:\\Users\\Young.Y\\Documents\\Julia\\SolveCollisionSection\\src")

using Plots

export Units, SetUnits, eV, meV, keV, nm, Å, Light, Me, ħ, ℎ, AddUnits
export SetPotential, SetWaveFunction, SolveWaveFunction, GetPhaseShift, GetGrossCrossSection
export funcpoint, WaveFunctionSetting

# Basic data and data structures
include("module-Units.jl");                             using .Units
include("module-PotentialFunction.jl");                 using .PotentialFunction
include("module-WaveFunction.jl");                      using .WaveFunction
include("module-Solves.jl");                            using .Solves

greet() = print("Hello, this is a programme to solve collision section!")

end # module SolveCollisionSection
