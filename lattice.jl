type Site
    id :: Int
    brand :: Symbol
end

type Edge1D
    id :: Int
    brand :: Symbol
    site1 :: Site
    site2 :: Site
    offset :: Int
end

type UnitCell1D
    sites :: Vector{Site}
    edges :: Vector{Edge1D}
end

type Lattice1D
    unitcell :: UnitCell1D
    size :: Int
    boundary_condition :: Symbol
end
