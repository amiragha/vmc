# type definition for Edge and Lattice
type Site
    brand :: Symbol
end

type Edge
    brand :: Symbol
    site1 :: Int
    site2 :: Int
end

type UCEdge1D
    brand :: Symbol
    site1 :: Int
    site2 :: Int
    offset :: Int
end

abstract Lattice

type Lattice1D <: Lattice
    Lx :: Int
    unitcell_sites :: Vector{Site}
    unitcell_edges :: Vector{UCEdge1D}
    sites :: Vector{Site}
    edges :: Vector{Edge}
    boundary :: Symbol
end

# TODO: possibly combine the next two functions
"""
make_periodic_lattice(Lx::Int, uc_sites::Vector{Site}, uc_edges::Vector{UCEdge1D})

returns a periodic lattice by making the full sites and edges vector
using the unitcell `uc_sites` and `uc_edges`.
"""
function make_periodic_lattice(Lx::Int,
                               uc_sites::Vector{Site},
                               uc_edges::Vector{UCEdge1D})

    num_uc_sites = length(uc_sites)
    num_uc_edges = length(uc_edges)

    sites = Vector{Site}()
    edges = Vector{Edge}()

    # add all the normal sites
    for l=1:Lx
        for s=1:num_uc_sites
            push!(sites, Site(uc_sites[s].brand))
        end
    end

    # add all the edges
    for l=1:Lx
        for e=1:num_uc_edges
            edge = uc_edges[e]
            brand = edge.brand
            lnn = l - 1 + edge.offset
            site1 = num_uc_sites * (l-1) + edge.site1
            site2 = num_uc_sites * (lnn%Lx)  + edge.site2
            push!(edges, Edge(brand, site1, site2))
        end
    end

    return Lattice1D(Lx, uc_sites, uc_edges, sites, edges, :periodic)
end

"""
make_open_lattice(Lx::Int, uc_sites::Vector{Site}, uc_edges::Vector{UCEdge1D}, obc_prescription::Vector{Int})

returns an open lattice by making the full sites and edges vector
using the unitcell `uc_sites` and `uc_edges`. The boundary condition
is handled using the prescription for obc in `obc_prescription`
which say how many extra sites should be added.
"""
function make_open_lattice(Lx::Int,
                           uc_sites::Vector{Site},
                           uc_edges::Vector{UCEdge1D},
                           obc_prescription::Vector{Int})

    num_uc_sites = length(uc_sites)
    num_uc_edges = length(uc_edges)
    @assert(length(obc_prescription) == num_uc_sites)
    extras = find(x->x==1, obc_prescription)

    sites = Vector{Site}()
    edges = Vector{Edge}()

    # add all the normal sites
    for l=1:Lx
        for s=1:num_uc_sites
            push!(sites, Site(uc_sites[s].brand))
        end
    end
    # add extra sites according to the prescription
    for s=extras
        push!(sites, Site(uc_sites[s].brand))
    end

    # add all the edges
    for l=1:Lx
        for e=1:num_uc_edges
            edge = uc_edges[e]
            brand = edge.brand
            lnn = l - 1 + edge.offset
            site1 = num_uc_sites * (l-1) + edge.site1
            if lnn < Lx
                site2 = num_uc_sites * (lnn)  + edge.site2
                push!(edges, Edge(brand, site1, site2))
            elseif lnn == Lx
                site2_idx = findfirst(extras, edge.site2)
                if site2_idx > 0
                    site2 = num_uc_sites * Lx + site2_idx
                    push!(edges, Edge(brand, site1, site2))
                end
            end
        end
    end
    return Lattice1D(Lx, uc_sites, uc_edges, sites, edges, :open)
end

"""
chain_lattice(L::Int, boundary::Symbol)

Definition of the 1D chain lattice. Generates lattice
for chain of size `L`. The current supported boundary
conditions are `:periodic` and `:open`.

returns Lattice1D
"""
function chain_lattice(Lx::Int, boundary::Symbol)

    # 3 sites per unit cell, the middle site brand is different
    uc_sites = [
        Site(:s1),
    ]

    # 6 edges per unit cell, :e1 is leg brand, :e2 is cross brand
    uc_edges = [
        UCEdge1D(:e1, 1, 1, 1),
    ]

    if boundary == :periodic
        return make_periodic_lattice(Lx, uc_sites, uc_edges)

    elseif boundary == :open
        return make_open_lattice(Lx, uc_sites, uc_edges,  [0])
    end
end

"""
kagomestrip_lattice(L::Int, boundary::Symbol)

Definition of the kagome strip Leg-Cross lattice. Generates lattice
for the kagome strip of size `L`. The current supported boundary
conditions are `:periodic`, `:opensym` and `:open`. `open` is just
normal OBC, which is L unitcells. `opensym` adds one extra middle
site to make the open lattice symmetric.

returns Lattice1D
"""
function kagomestrip_lattice(Lx::Int, boundary::Symbol)

    # 3 sites per unit cell, the middle site brand is different
    uc_sites = [
        Site(:s1),
        Site(:s2),
        Site(:s1)
    ]

    # 6 edges per unit cell, :e1 is leg brand, :e2 is cross brand
    uc_edges = [
        UCEdge1D(:e2, 1, 2, 0),
        UCEdge1D(:e2, 3, 2, 0),
        UCEdge1D(:e2, 1, 2, 1),
        UCEdge1D(:e2, 3, 2, 1),
        UCEdge1D(:e1, 1, 1, 1),
        UCEdge1D(:e1, 3, 3, 1)
    ]

    if boundary == :periodic
        return make_periodic_lattice(Lx, uc_sites, uc_edges)

    elseif boundary == :opensym
        return make_open_lattice(Lx, uc_sites, uc_edges,  [0, 1, 0])

    elseif boundary == :open
        return make_open_lattice(Lx, uc_sites, uc_edges,  [0, 0, 0])
    end
end
