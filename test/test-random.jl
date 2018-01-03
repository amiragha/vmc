@testset "tools" begin
    # test random_ones
    n = 12
    m = rand(1:n)
    vector = Tmp.random_ones(n, m)
    @test length(find(x->x==1, vector)) == m
    @test length(find(x->x==0, vector)) == n - m
end
