type Measurement{T}
    num_operators :: Int
    data :: Matrix{T}
    batch_averages :: Matrix{T}
    batch_size ::Int
    measure_capacity :: Int
    current_batch :: Int
    measures_recorded :: Int

    function Measurement(what::Symbol, capacity::Int, num_sites::Int, batch_size::Int=1000)
        @assert div(capacity, batch_size) * batch_size == capacity
        if what == :correlation
            num_operators = div(num_sites*(num_sites-1),2)
            data = zeros(num_operators, batch_size)
            avgs = zeros(num_operators, div(capacity, batch_size))
            return new(num_operators, data, avgs, batch_size, capacity, 1, 0)
        end
    end
end

function measure!(ms::Measurement,
                 wf::WaveFunction)
    @assert ms.measures_recorded < ms.measure_capacity
    num_sites = size(wf.states)[1]
    N_tot = ms.num_operators

    correlations = zeros(Float64, ms.num_operators)
    index = 1
    for i=1:num_sites, j=i+1:num_sites
        ni = wf.configuration.instance[i] - 1/2
        nj = wf.configuration.instance[j] - 1/2
        correlations[index] = ni * nj
        index += 1
    end
    ms.measures_recorded += 1
    index_in_batch = (ms.measures_recorded % ms.batch_size) + 1
    ms.data[:,index_in_batch] = correlations
    if index_in_batch == ms.batch_size
        ms.batch_averages[:,ms.current_batch] =
            mean(ms.data[:,ms.batch_size],2)
        ms.current_batch += 1
        ms.data = zeros(N_tot, ms.batch_size)
    end
end

function report(measurement::Measurement)
    h5write("output.h5", "correlations", mean(measurement.batch_averages,2))
end
