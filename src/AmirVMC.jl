module AmirVMC

using HDF5

include("Lattice.jl")
include("Model.jl")
include("DetMatrix.jl")
include("random.jl")
include("free_states.jl")
include("WaveFunction.jl")
include("measurement.jl")
include("MonteCarlo.jl")

export chain
export kagomestrip_LC

export runVMC
end
