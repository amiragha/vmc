include("include.jl")

function make_periodic_chain_model(L::Int)
    sites = [
        Site(1, :s1)
    ]

    edges = [
        Edge1D(1, :e1, sites[1], sites[1], 1)
    ]

    unitcell  = UnitCell1D(sites, edges)
    lattice = Lattice1D(unitcell, L, :periodic)

    amplitudes_1site = Dict(:s1 => 0)
    amplitudes_2site = Dict(:e1 => complex(-1))

    modelname = "chain"
    return Model(modelname,
                 lattice,
                 amplitudes_1site,
                 amplitudes_2site)
end

model = make_periodic_chain_model(24)
@time runVMC(model, 200000)
