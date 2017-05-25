include("sherman_morrison.jl")

mat = rand(Complex128, 20, 20)
col = rand(Complex128, 20)

det1 = new_determinant_MDL(mat, inv(mat), det(mat), col, 7)
inv1 = sherman_morrison_1col(mat, inv(mat), det(mat), col, 7)
mat[:,7] = col

det2 = det(mat)
inv2 = inv(mat)

println(det1, det2)
println(vecnorm(inv2-inv1))
