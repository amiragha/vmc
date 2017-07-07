"""
solve_free_full(model::Model, num_filled::Int)

solves the eigen problem for the free hopping hamiltonian for
`model`. If the Hamiltonian is real (hoppings are real) this gives
real states.

returns `num_filled` eigenstates in a matrix. Each state is a column,
so the matrix is of dimensions: volume x num_filled.
"""
function solve_free_full(model::Model, num_filled::Int)

    # parameter
    sites = model.lattice.sites
    num_sites = length(sites)
    edges = model.lattice.edges

    ### Generate the one particle Hamiltonian
    hamiltonian = zeros(Complex128, num_sites, num_sites)

    # onesite chemical potentials
    for s=1:num_sites
        mui = model.amplitudes_1site[sites[s].brand]
        hamiltonian[s, s] = mui
    end

    # two site hoppings
    for e=1:length(edges)
        t = model.amplitudes_2site[edges[e].brand]

        i = edges[e].site2
        j = edges[e].site1

        hamiltonian[i, j] += t
        hamiltonian[j, i] += conj(t)
    end

    return eigfact(Hermitian(hamiltonian), 1:num_filled)[:vectors]
end

"""
planewave(Lx::Int, k::Int, x::Int)

returns the value 1D plane wave wavefunction of size `L` and
wavevector `k` at position `x` .
"""
planewave(Lx::Int, k::Int, x::Int) = exp(im * (2*pi*k/Lx) * x)/sqrt(Lx)
planewave(Lx::Int, k::Float64, x::Int) = exp(im * (2*pi*k/Lx) * x)/sqrt(Lx)

"""
solve_free_periodic(model::Model, num_filled::Int, periodicity::Symbol=:PBC)

solves the eigen problem for the free hopping hamiltonian for
`model`. This version solves in the momentum space, so it only works
for periodic lattices.

returns `num_filled` eigenstates in a matrix. Each state is a column,
so the matrix is of dimensions: volume x num_filled.
"""
function solve_free_periodic(model::Model,
                             num_filled::Int,
                             periodicity::Symbol=:PBC)

    @assert model.lattice.boundary == :periodic
    @assert (periodicity == :PBC || periodicity == :APBC)

    # NOTE: here sites and number of site refer to the unit cell!
    sites = model.lattice.unitcell_sites
    num_sites = length(sites)
    edges = model.lattice.unitcell_edges
    Lx = model.lattice.Lx

    energies = zeros(Float64, num_sites, Lx)
    eigenstates = zeros(Complex128, num_sites, num_sites, Lx)

    for k_index=0:Lx-1
        periodicity == :PBC ? (k = k_index) : (k = k_index + 0.5)

        k_hamiltonian = zeros(Complex128, num_sites, num_sites)

        # 1site chemical potentials
        for s=1:num_sites
            mui = model.amplitudes_1site[sites[s].brand]
            k_hamiltonian[s, s] = mui
        end

        # 2site hoppings
        for e=1:length(edges)
            t = model.amplitudes_2site[edges[e].brand]
            offset = edges[e].offset
            if offset != 0
                t *= exp(-im * (2*pi*k/Lx) * edges[e].offset)
            end

            i = edges[e].site2
            j = edges[e].site1

            k_hamiltonian[i, j] += t
            k_hamiltonian[j, i] += conj(t)
        end

        eigenresult = eigfact(Hermitian(k_hamiltonian))
        energies[:,k_index+1] = eigenresult[:values]
        eigenstates[:,:,k_index+1] = eigenresult[:vectors]
    end

    states = zeros(Complex128, num_sites * Lx, num_filled)

    # TODO: explain the sorting process in words
    # sort energies and fill states with lowest energy eigenvectors
    for index = 1:num_filled
        sorted_index = sortperm(reshape(energies,
                                        length(energies)))[index]

        i, j = ind2sub(energies, sorted_index)

        #println("lowest energy indeces = ", i,", ", j)
        # TODO: more elegant code! should be able to kron full states
        periodicity == :PBC ? (k = j-1) : (k = j-1 + 0.5)
        k_state = [ planewave(Lx, k, x) for x=1:Lx ]
        states[:, index] = kron(k_state, eigenstates[:, i, j])
    end
    return states
end
