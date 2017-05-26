planewave(L, k, x) = exp(im*k*x)/sqrt(L)

# NOTE: for now only works with periodic
function solve_hopping(model::Model)

    sites = model.lattice.unitcell.sites
    edges = model.lattice.unitcell.edges
    n_sites = length(sites)
    L = model.lattice.size

    energies = zeros(Float64, n_sites, L)
    eigenstates = zeros(Complex128, n_sites, n_sites, L)

    ks = 2*pi*collect(0:L-1)/L

    # NOTE: assuming the site ids are in order starting from 1
    for kidx=0:L-1
        ft_Hamiltonian = zeros(Complex128, n_sites, n_sites)
        for i=1:n_sites
            mui = model.amplitudes_1site[sites[i].brand]
            ft_Hamiltonian[i,i] = mui
        end

        for idx=1:length(edges)
            t = model.amplitudes_2site[edges[idx].brand]
            offset = edges[idx].offset
            if offset != 0
                t *= exp(-im * (2*pi*kidx/L) * edges[idx].offset)
            end
            i = edges[idx].site2.id
            j = edges[idx].site1.id

            ft_Hamiltonian[i,j] += t
            ft_Hamiltonian[j,i] += conj(t)
        end

        println("ft_Hamiltonian", ft_Hamiltonian)
        eigenresult = eigfact(Hermitian(ft_Hamiltonian))
        energies[:,kidx+1] = eigenresult[:values]
        eigenstates[:,:,kidx+1] = eigenresult[:vectors]
    end

    return energies, eigenstates
end

# return states : space volume x N_occupied
# the space part is column first
function meanfield_states(model::Model, N_occupied::Int)
    energies, eigenstates = solve_hopping(model)
    println("energies = ", energies)

    L = model.lattice.size
    states = zeros(Complex128,
                   length(model.lattice.unitcell.sites) * L, N_occupied)

    # sort energies, fill states with eigenvectors of lowest eigenvalues
    for index = 1:N_occupied
        sorted_index = sortperm(reshape(energies, length(energies)))[index]
        i, j = ind2sub(energies, sorted_index)
        println("lowest energy indeces = ", i,", ", j)
        k_state = [ planewave(L, 2*pi*(j-1)/L, x) for x=1:L ]
        states[:, index] = kron(k_state, eigenstates[:, i, j])
    end
    println(states)
    return states
end
