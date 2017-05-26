type WaveFunction
    states :: Matrix{Complex128}
    configuration :: Configuration

    uplist :: Vector{Int}
    dnlist :: Vector{Int}

    slater :: SlaterDet2
end

function random_wavefunction(model::Model, states::Matrix{Complex128})
    N_occupied = size(states)[end]
    configuration = random_configuration(model.lattice, N_occupied)

    uplist = find(spin_isup, configuration.instance)
    dnlist = find(spin_isdn, configuration.instance)

    println(configuration, uplist, dnlist)

    slater =
        SlaterDet2([states[x,k] for k=1:N_occupied, x=uplist],
                   [states[x,k] for k=1:N_occupied, x=dnlist])

    return WaveFunction(states, configuration,
                        uplist, dnlist, slater)
end

function acceptance_probability(old_slater, new_slater)
    old_p = norm(old_slater.det_up * old_slater.det_dn)^2
    new_p = norm(new_slater.det_up * new_slater.det_dn)^2
    return min(1, new_p/old_p)
end
