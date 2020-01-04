@testset "tools" begin
    # test random_ones
    n = 12
    m = rand(1:n)
    vector = VMC.random_ones(n, m)
    @test length(filter(x->x==1, vector)) == m
    @test length(filter(x->x==0, vector)) == n - m
end
