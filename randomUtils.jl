function select_random(N::Int, m::Int)
    @assert N > 0 && m < N && m > 0

    picks_set = IntSet()
    for pivot=(N-m+1):N
        pick = rand(1:N-m+1)
        if in(pick, picks_set)
            push!(picks_set, pivot)
        else
            push!(picks_set, pick)
        end
    end
    return picks_set
end
