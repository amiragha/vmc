mutable struct Measurement{T} where {T<:Number}
    num_operators :: Int
    accdata :: Vector{T}
    measure_capacity :: Int
    measures_recorded :: Int

    function Measurement(what::Symbol, capacity::Int, num_sites::Int)
        if what == :correlation
            num_operators = div(num_sites*(num_sites-1),2)
            accdata = zeros(num_operators)
            return new(num_operators, accdata, capacity, 0)
        end
    end
end

function measure!(ms::Measurement,
                  wf::GutzwillerSlater)
    @assert ms.measures_recorded < ms.measure_capacity
    num_sites = size(wf.states)[1]

    correlations = zeros(Float64, ms.num_operators)
    index = 1
    for i=1:num_sites, j=i+1:num_sites
        ni = wf.configuration[i] - 1/2
        nj = wf.configuration[j] - 1/2
        correlations[index] = ni * nj
        index += 1
    end
    ms.measures_recorded += 1
    ms.accdata += correlations
end

function report(measurement::Measurement, modelname::String)
    random_sim_id = rand(1000:9999)
    outputfile = string("vmc", random_sim_id, ",", modelname, ".h5")
    h5write(outputfile,
            "correlations",
            measurement.accdata/measurement.measures_recorded)
end
