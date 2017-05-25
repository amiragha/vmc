type Measurement{T}
    num_operators :: Int
    data :: Matrix{T}
    measure_capacity :: Int
    measures_recorded :: Int

    function Measurement(what::Symbol, N::Int, num_sites::Int)
        if what == :correlation
            num_operators = div(num_sites*(num_sites-1),2)
            data = zeros(num_operators, N)
            return new(num_operators, data, N, 0)
        end
    end
end

function measure(ms::Measurement,
                 wf::WaveFunction)
    @assert ms.measures_recorded < ms.measure_capacity
    num_sites = size(wf.states)[1]
    N_tot = ms.num_operators
    # assuming N_tot is divisible by two
    correlations = zeros(Float64, ms.num_operators)
    index = 1
    for i=1:num_sites, j=i+1:num_sites
        ni = wf.configuration.instance[i] - 1/2
        nj = wf.configuration.instance[j] - 1/2
        correlations[index] = ni * nj
        index += 1
    end
    ms.measures_recorded += 1
    return ms.data[:,ms.measures_recorded] = correlations
end

function report(measurement::Measurement)
    return mean(measurement.data,2)
end
