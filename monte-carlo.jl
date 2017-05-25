function step!(wf::WaveFunction)
    # choose randomly from uplist and dnlist
    up_index = rand(1:length(wf.uplist))
    dn_index = rand(1:length(wf.dnlist))

    new_slater = slater_swap_column(wf.slater, up_index, dn_index)

    if rand(Float64) < acceptance_probability(wf.slater,
                                              new_slater)
        # update wavefunction fields
        swap2_config(wf.configuration,
                     wf.uplist[up_index], wf.dnlist[dn_index])
        wf.uplist[up_index], wf.dnlist[dn_index] =
            wf.dnlist[dn_index], wf.uplist[up_index]
        wf.slater = new_slater
        return true
    else
        return false
    end
end

function runVMC(model::Model, N_steps=20000)

    num_sites = model.lattice.size *
        length(model.lattice.unitcell.sites)

    # NOTE: assuming total number of sites are even
    N_occupied = div(num_sites, 2)

    states = meanfield_states(model, N_occupied)
    wavefunction = random_wavefunction(model, states)

    half_correlations = Measurement{Float64}(:correlation, N_steps, num_sites)
    step = 1
    total_num_of_tries = 0
    while step < N_steps + 1
        println(step)
        if step!(wavefunction)
            measure(half_correlations, wavefunction)
            step += 1
        end
        total_num_of_tries += 1
    end
    println(total_num_of_tries, N_steps)
    plot_structure(report(half_correlations), num_sites, :chain)
end
