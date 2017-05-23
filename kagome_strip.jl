include("lattice.jl")
include("model.jl")
include("hopping.jl")
include("wavefunction.jl")
include("monte-carlo.jl")

function make_kagome_strip_model(L::Int, mu::Float64, tc::Complex128)
    sites = [
        Site(1,:s1),
        Site(2,:s2),
        Site(3,:s1)
    ]

    edges = [
        Edge1D(1, :e2, sites[1], sites[2], 0),
        Edge1D(2, :e2, sites[3], sites[2], 0),
        Edge1D(3, :e2, sites[1], sites[2], 1),
        Edge1D(4, :e2, sites[3], sites[2], 1),
        Edge1D(5, :e1, sites[1], sites[1], 1),
        Edge1D(6, :e1, sites[3], sites[3], 1),
    ]

    unitcell = UnitCell1D(sites, edges)
    lattice = Lattice1D(unitcell, L, :periodic)

    amplitudes_1site = Dict(:s1 => 0,
                            :s2 => mu)
    amplitudes_2site = Dict(:e1 => complex(1),
                            :e2 => tc)

    modelname = "kagomestrip,mu=$(mu),tc=$(tc)"
    return Model(modelname,
                 lattice,
                 amplitudes_1site,
                 amplitudes_2site)
end

model = make_kagome_strip_model(6,1.,complex(2.))
runVMC(model)
