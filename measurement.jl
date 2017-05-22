function measure(wavefunction::WaveFunction, what::Symbol)
    if what == :correlation
        L = wavefunction.L
        correlations = zeros(div(L*(L-1),2))
        index = 1
        for i=1:L, j=i+1:L
            ni = wavefunction.configuration[i] - 1/2
            nj = wavefunction.configuration[j] - 1/2
            correlations[index] = ni * nj
            index += 1
        end
        return correlations
    end
end
