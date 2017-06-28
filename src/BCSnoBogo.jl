# Gutzwiller projected slater determinant wavefunction
# TODO: Does Gutzwiller need to carry the states?
type BCSnoBogo
    states :: Matrix{Complex128}
    configuration :: Vector{Int}

    olist :: Vector{Int}
    xlist :: Vector{Int}

    slater :: DetMatrix
end

"""
random_nobogo_half(states::Matrix{Complex128})

Generates a nobogo wavefunction with randomly chosen half of sites to
be occupied by both d1 and d2 fermions. The input, states, is a matrix
with double space orbitals (v*,-u*), corresponding to half of the
negative eigenvalues of Bogoliubov matrix, as columns. The number of
orbitals(columns) has to be exactly half of the number of
sites(rows/2).

returns a BCSnoBogo wavefunction
"""
function random_nobogo_half(states::Matrix{Complex128})
    num_states = size(states)[2]

    # even lattice size and half-filling
    @assert 4*num_states == size(states)[1]

    num_sites = div(size(states)[1], 2)

    ## TODO: can these be simplified in one function?
    configuration = random_ones(num_sites, num_states)
    olist = find(x-> x==1, configuration)
    xlist = find(x-> x==0, configuration)

    ## TODO: should be possible to slice them from states, directly!
    # NOTE: while the double states(orbitals) are originally given on
    # cols. Here we put states on rows and positions on
    # columns. Notice the order of the for expression!
    slater = DetMatrix(
        hcat(
            [states[x,           k] for k=1:num_states, x=olist],
            [states[x+num_sites, k] for k=1:num_states, x=olist]
        )
    )

    return BCSnoBogo(states, configuration,
                     uplist, slater)
end

"""
propose_step!(wf::BCSnoBogo)

propose and, if accepted, update a step for the BCSnoBogo
wavefunction. The process is to choose randomly one occupied site
and one unoccupied one.

The acceptance probability is the minimum of 1 and square of the
determinant ratio.

returns `true` if updates the `wf`, and `false` otherwise.
"""
function propose_step!(wf::BCSnoBogo)
    # TODO: split these operations into multiple functions.
    # randomly choose one occupied and one unoccupied site
    # NOTE: `o/x_index` are random indeces in `o/xlist` pointing to
    # positions, namely, `o/x_site` in the wf.configuration!
    o_index = rand(1:length(wf.olist))
    x_index = rand(1:length(wf.xlist))
    o_site = wf.uplist[o_index]
    x_site = wf.dnlist[x_index]

    # new columns, i.e. just swap up and dn
    new_up_col = wf.dnslater.matrix[:, dn_index]
    new_dn_col = wf.upslater.matrix[:, up_index]

    # find the probability of the change
    up_detrat = det_ratio_1col(wf.upslater, new_up_col, up_index)
    dn_detrat = det_ratio_1col(wf.dnslater, new_dn_col, dn_index)

    # up_detrat = exact_det_ratio_1col(wf.upslater, new_up_col, up_index)
    # dn_detrat = exact_det_ratio_1col(wf.dnslater, new_dn_col, dn_index)

    if rand(Float64) < min(1, norm(up_detrat * dn_detrat)^2)
        # update the Gutzwiller wavefunction
        wf.configuration[up_site] = 0
        wf.configuration[dn_site] = 1

        wf.uplist[up_index], wf.dnlist[dn_index] = dn_site, up_site

        update_detmatrix_1col!(wf.upslater,
                               new_up_col, up_index, up_detrat)
        update_detmatrix_1col!(wf.dnslater,
                               new_dn_col, dn_index, dn_detrat)

        # exact_update_detmatrix_1col!(wf.upslater,
        #                        new_up_col, up_index, up_detrat)
        # exact_update_detmatrix_1col!(wf.dnslater,
        #                        new_dn_col, dn_index, dn_detrat)

        return true
    else
        return false
    end
end
