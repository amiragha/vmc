type Configuration
    instance :: Vector{Int}
end

spin_isup(x) = x == 1
spin_isdn(x) = x == 0

function random_configuration(lattice :: Lattice1D, N_occupied)
    Lx = lattice.size
    Ly = length(lattice.unitcell.sites)

    instance = zeros(Int, Ly*Lx)

    selects = select_random(Ly*Lx, N_occupied)

    for i=selects
        instance[i] = 1
    end
    return Configuration(instance)
end

function swap2_config(config::Configuration, up_idx, dn_idx)
    @assert config.instance[up_idx] == 1
    @assert config.instance[dn_idx] == 0

    config.instance[up_idx] = 0
    config.instance[dn_idx] = 1
end
