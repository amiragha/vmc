    using Plots
using HDF5

data = h5read("vmc4274,kagomestripLC,Lx=32periodic,mu=-2.5,tcr=1.0,tcphi=0.0.h5", "correlations")
#println(data)

function index_of(m, n, N)
    if m > n
        return index_of(n, m, N)
    else
        return div((2*N-m)*(m-1),2) + (n-m)
    end
end

function plot_structure(data::Vector{Float64}, N::Int, pattern::Symbol):
    @assert length(data) == div(N*(N-1),2)
    avg_sfactor = zeros(Float64, div(N,2)+1)
    correlations = eye(Float64, N,N) * 0.25
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
                    i_middle = 3*(i-1) + 2
                    j_middle = 3*(j-1) + 2
                    correlations[i,j] = data[index_of(i_middle,j_middle,3*N)]
                end
            end
        end
        ys = zeros(Float64, div(N,2)+1)
        for i=1:N
            ys = rfft(circshift(correlations[:, i],1-i))
            avg_sfactor += real(ys)
            #plot!(real(ys))
        end
        plot(3*avg_sfactor/N)
        savefig("plot.pdf")
    end
end

plot_structure(data, 32, :kagome)
