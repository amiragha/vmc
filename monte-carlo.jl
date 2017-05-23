function move!(wf::WaveFunction)
    # choose a up and dn site
    chosen_up = rand(find(spin_isup, wf.configuration))
    chosen_dn = rand(find(spin_isdn, wf.configuration))

    new_wf = WaveFunction(wf, chosen_up, chosen_dn)

    if rand(Float64) < acceptance_probability(wf, new_wf)
        wf = new_wf
        return true
    else
        return false
    end
end

function runVMC(model::Model)
    # NOTE: assuming total number of sites are even
    N_occupied =
        div(model.size * length(model.lattice.unitcell.sites), 2)

    states = meanfield_states(model, N_occupied)

end
