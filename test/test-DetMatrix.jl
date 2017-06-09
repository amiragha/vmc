# random matrix, vector and index1
randmat = rand(Complex128, 10, 10)
mat1 = copy(randmat)
col = rand(Complex128, 10)
col_index = rand(1:10)
mat2 = copy(mat1)
mat2[:,col_index] = col

detmat1 = AmirVMC.DetMatrix(mat1)
detmat2 = AmirVMC.DetMatrix(mat2)

# test determinant
detratio_exact = detmat2.determinant/detmat1.determinant
detratio_exact_fn  = AmirVMC.exact_det_ratio_1col(detmat1, col, col_index)
detratio = AmirVMC.det_ratio_1col(detmat1, col, col_index)
@test isapprox(detratio_exact, detratio)
@test isapprox(detratio_exact, detratio_exact_fn)

# test exact update
AmirVMC.exact_update_detmatrix_1col!(detmat1, col, col_index, detratio)
@test isapprox(detmat1.inverse, detmat2.inverse)
@test isapprox(detmat1.matrix, detmat2.matrix)
@test isapprox(detmat1.determinant, detmat2.determinant)

# restore mat1
mat1 = copy(randmat)
detmat1 = AmirVMC.DetMatrix(mat1)

# test the update
AmirVMC.update_detmatrix_1col!(detmat1, col, col_index, detratio)
@test isapprox(detmat1.inverse, detmat2.inverse)
@test isapprox(detmat1.matrix, detmat2.matrix)
@test isapprox(detmat1.determinant, detmat2.determinant)
