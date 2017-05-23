include("lattice.jl")
include("hopping.jl")

function make_kagome_strip_model(L::Int)
    sites = [
        Site(1,1),
        Site(2,2),
        Site(3,1)
    ]

    edges = [
        Edge1D(1, 2, sites[1], sites[2], 0),
        Edge1D(2, 2, sites[3], sites[2], 0),
        Edge1D(3, 2, sites[1], sites[2], 1),
        Edge1D(4, 2, sites[3], sites[2], 1),
        Edge1D(5, 1, sites[1], sites[1], 1),
        Edge1D(6, 1, sites[3], sites[3], 1),
    ]

    unitcell = UnitCell1D(sites, edges)
    lattice = Lattice1D(unitcell, L, :periodic)

    onesite_interaction_brands = [0, 1]
    twosite_interaction_brands = [1, 2] .* (1+0im)

    return Model(lattice,
                 onesite_interaction_brands,
                 twosite_interaction_brands)
end

model = make_kagome_strip_model(6)

energies, eigenstates = solve_hopping(model)
println(energies)
println(eigenstates)
