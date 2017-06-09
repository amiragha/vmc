# The Model type definition
type Model
    name :: String
    lattice :: Lattice1D
    amplitudes_1site :: Dict{Symbol, Float64}
    amplitudes_2site :: Dict{Symbol, Complex128}
end

"""
chain(L:Int, boundary::Symbol, mu::Float64, t:Complex128)

return 1D chain model with possible boundary conditions
`:periodic` :open`. `open` is just normal OBC, which
is L unitcells.

the parameters `mu` is the chemical potential and
`tc` is the hopping amplitude.

returns Model.
"""
function chain(Lx::Int,
               boundary::Symbol,
               mu::Float64,
               t::Complex128)

    name = "chain,Lx$Lx$boundary,mu=$mu,t=$(norm(t))"

    lattice = chain_lattice(Lx, boundary)

    amplitudes_1site = Dict(:s1 => mu)

    amplitudes_2site = Dict(:e1 => t)

    return Model(name, lattice,
                 amplitudes_1site, amplitudes_2site)
end

"""
kagomestrip_LC(L:Int, boundary::Symbol, mu::Float64, tc:Complex128)

return the kagome strip Leg-Cross model with boundary conditions
`:periodic`, `:opensym` and `:open`. `open` is just normal OBC, which
is L unitcells. `opensym` adds one extra middle site to make the open
lattice symmetric.

the parameters `mu` is the chemical potential on the middle sites and
`tc` is the hopping on the cross bonds. `tl` is set to `1.0` by
default.

returns Model.
"""
function kagomestrip_LC(Lx::Int,
                        boundary::Symbol,
                        mu::Float64,
                        tc::Complex128)

    name = "kagomestripLC,Lx=$Lx$boundary,mu=$mu,tcr=$(norm(tc)),tcphi=$(angle(tc))"

    lattice = kagomestrip_lattice(Lx, boundary)

    amplitudes_1site = Dict(:s1 => 0,
                            :s2 => mu)
    amplitudes_2site = Dict(:e1 => complex(1),
                            :e2 => tc)

    if boundary == :opensym
        return Model(name, lattice,
                     amplitudes_1site, amplitudes_2site)
    else
        return Model(name, lattice,
                     amplitudes_1site, amplitudes_2site)
    end
end
