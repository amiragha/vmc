include("wavefunction.jl")
include("monte-carlo.jl")
include("measurement.jl")

wf = WaveFunction(6)

println(measure(wf, :correlation))

N = 2000
window = zeros(Float64, 15, N)
index = 1
while index < (N+1)
    if move!(wf)
        window[:, index] = measure(wf, :correlation)
        index += 1
    end
end
println(mean(window,2))
println(std(window,2))
