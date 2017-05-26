function step!(wf::WaveFunction)
    #@show wf
    # choose randomly from uplist and dnlist
    #@show
    up_index = rand(1:length(wf.uplist))
    #@show
    dn_index = rand(1:length(wf.dnlist))
    #@show
    up_pos = wf.uplist[up_index]
    #@show
    dn_pos = wf.dnlist[dn_index]


    # make the new slater determinants
    #@show
    new_col_up = wf.states[dn_pos,:]
    #@show
    new_col_dn = wf.states[up_pos,:]

    #@show
    new_slater = slater_change_column(wf.slater, up_index, dn_index,
                                      new_col_up, new_col_dn)

    if rand(Float64) < acceptance_probability(wf.slater,
                                              new_slater)
        # update wavefunction fields
        swap2_config!(wf.configuration,
                      up_pos, dn_pos)
        wf.uplist[up_index], wf.dnlist[dn_index] = dn_pos, up_pos
        wf.slater = new_slater

        #@show wf
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
        if step!(wavefunction)
            println(step)
            measure!(half_correlations, wavefunction)
            step += 1
        end
        total_num_of_tries += 1
    end
    println(total_num_of_tries, N_steps)
    report(half_correlations)
end
