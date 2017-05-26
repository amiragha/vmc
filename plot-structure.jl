    using Plots
using HDF5

data = h5read("output.h5", "correlations")
println(data)

function index_of(m, n, N)
    if m > n
        return index_of(n, m, N)
    else
        return div((2*N-m)*(m-1),2) + (n-m)
    end
end

function plot_structure(data, N, pattern):
    @assert length(data) == div(N*(N-1),2)
    correlations = eye(N,N) * 0.25
    if pattern == :chain
        for i = 1:N
            for j = 1:N
                #println(i, ", ",j , " = ", index_of(i,j,N))
                if i != j
                    correlations[i, j] = data[index_of(i,j,N)]
                end
            end
        end
        ys = zeros(Float64, div(N,2)+1)
        for i=1:N
            ys = rfft(circshift(correlations[:, i],1-i))
            plot!(real(ys))
        end
        savefig("plot.pdf")
    end
    if pattern == :kagome
        for i = 1:N
            for j = 1:N
                if i != j
                    correlations[i,j] = data[index_of(3*(i-1)+1,3*(j-1)+1,3*N)]
                end
            end
        end
        ys = zeros(Float64, div(N,2)+1)
        for i=1:N
            ys = rfft(circshift(correlations[:, i],1-i))
            plot!(real(ys))
        end
        savefig("plot.pdf")
    end
end

plot_structure(data, 24, :chain)
