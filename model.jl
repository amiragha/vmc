type Model
    name :: String
    lattice :: Lattice1D
    amplitudes_1site :: Dict{Symbol, Float64}
    amplitudes_2site :: Dict{Symbol, Complex128}
end
