"""
solve_bcs_nobogo(model::Model, delta::Float64)

solves for the eigenstates of Bogoliubov matrix. This can be
interpreted as pure hopping as in *avoiding Bogoliubov* algorithm or
simply be used for u, v function to generate Bogoliubov particles. The
pairing amplitude, delta, is considered constant for all k values for
now.

returns 2L*L matrix for states for negative eigenvalues only. each
column is a V* column followed by -U* column. If needed The
eigenstates for positive energies.
"""
function solve_bcs(model::Model, delta::Float64)
    ## TODO: find a way to present possible values for pairing
    # NOTE: because of the above for now I have to write model
    # specific code! currently only for kagome strip

    # generate the diagonal hopping part and duplicate

    # generate the off-diagonal pairing part
end
