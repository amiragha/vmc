include("lattice.jl")
include("hopping.jl")

function make_periodic_chain_model(L::Int)
    sites = [
        Site(1, 1)
    ]

    edges = [
        Edge1D(1, 1, sites[1], sites[1], 1)
    ]

    unitcell  = UnitCell1D(sites, edges)
    lattice = Lattice1D(unitcell, L, :periodic)

    onesite_interaction_brands = [0]
    twosite_interaction_brands = [-1] .* (1+0im)

    return Model(lattice,
                 onesite_interaction_brands,
                 twosite_interaction_brands)
end

model = make_periodic_chain_model(6)

energies, eigenstates = solve_hopping(model)
println(energies)
println(eigenstates)
