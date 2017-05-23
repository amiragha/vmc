# NOTE: for now only works with periodic
type Model
    lattice :: Lattice1D
    onesite :: Vector{Float64}
    twosite :: Vector{Complex128}
end

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
            ft_Hamiltonian[i,i] = model.onesite[sites[i].brand]
        end

        for idx=1:length(edges)
            amplitude = model.twosite[edges[idx].brand]
            offset = edges[idx].offset
            if offset != 0
                amplitude *= exp(-im * (2*pi*kidx/L) * edges[idx].offset)
            end
            i = edges[idx].site2.id
            j = edges[idx].site1.id
            ft_Hamiltonian[i,j] += amplitude
            ft_Hamiltonian[j,i] += conj(amplitude)
        end

        println(ft_Hamiltonian)
        eigenresult = eigfact(Hermitian(ft_Hamiltonian))
        energies[:,kidx+1] = eigenresult[:values]
        eigenstates[:,:,kidx+1] = eigenresult[:vectors]
    end

    return energies, eigenstates
end
