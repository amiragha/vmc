planewave(L, k, x) = exp(im*k*x)/sqrt(L)

spin_isup(x) = x == 1
spin_isdn(x) = x == 0

type WaveFunction

    L :: Int
    states :: Matrix{Complex128}

    configuration :: Vector{Int}

    up_slater :: Matrix{Complex{Float64}}
    dn_slater :: Matrix{Complex{Float64}}

    up_det :: Complex{Float64}
    dn_det :: Complex{Float64}

    function WaveFunction(L)
        planewaveL(k,x) = planewave(L,k,x)

        ks = [2*pi*k/L for k=-div(L,4):div(L,4)]

        # 1 represent up, 0 represent dn
        initial_configuration = [x%2 for x=1:L]

        up_slater = [planewaveL(k,x) for k=ks, x=find(spin_isup, initial_configuration)]
        dn_slater = [planewaveL(k,x) for k=ks, x=find(spin_isdn, initial_configuration)]

        up_det = det(up_slater)
        dn_det = det(dn_slater)
        return new(L, ks,
                   initial_configuration,
                   up_slater, dn_slater,
                   up_det, dn_det)
    end
    function WaveFunction(old_wf::WaveFunction, up_index, dn_index)
        planewaveL(k,x) = planewave(L,k,x)
        L = old_wf.L
        ks = old_wf.ks

        configuration = old_wf.configuration
        configuration[up_index] = 0
        configuration[dn_index] = 1

        up_slater = [planewaveL(k,x) for k=ks, x=find(spin_isup, configuration)]
        dn_slater = [planewaveL(k,x) for k=ks, x=find(spin_isdn, configuration)]

        up_det = det(up_slater)
        dn_det = det(dn_slater)
        return new(L, ks,
                   configuration,
                   up_slater, dn_slater,
                   up_det, dn_det)
    end
end

function acceptance_probability(old_wf, new_wf)
    old_p = norm(old_wf.up_det * old_wf.dn_det)
    new_p = norm(new_wf.up_det * new_wf.dn_det)
    return min(1, new_p/old_p)
end
