# The Model type definition
struct Model
    name :: String
    lattice :: Lattice1D
    amplitudes_1site :: Dict{Symbol, Float64}
    amplitudes_2site :: Dict{Symbol, ComplexF64}
end

"""
    chain(Lx, boundary, mu, t)

return 1D chain model with possible boundary conditions `:periodic`
:open`. `open` is just normal OBC, which is `Lx` unitcells.  The
parameters `mu` is the chemical potential and `t` is the hopping
amplitude respectively.

"""
function chain(Lx::Int,
               boundary::Symbol,
               mu::Float64,
               t::ComplexF64)

    name = "chain,Lx$Lx$boundary,mu=$mu,t=$(norm(t))"

    lattice = chain_lattice(Lx, boundary)

    amplitudes_1site = Dict(:s1 => mu)

    amplitudes_2site = Dict(:e1 => t)

    return Model(name, lattice,
                 amplitudes_1site, amplitudes_2site)
end

"""
    kagomestrip_LC(Lx, boundary, mu, tc)

return the kagome strip Leg-Cross model with boundary conditions
`:periodic`, `:opensym` and `:open`. `open` is just normal OBC, which
is L unitcells. `opensym` adds one extra middle site to make the open
lattice symmetric.  The parameters `mu` is the chemical potential on
the middle sites and `tc` is the hopping on the cross bonds.

"""
function kagomestrip_LC(Lx::Int,
                        boundary::Symbol,
                        mu::Float64,
                        tc::ComplexF64)

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
