# random matrix, vector and index1
randmat = rand(Complex128, 10, 10)
cols = rand(Complex128, 10, 2)
c1 = rand(1:10)
c2 = rand(1:10)

### 1 column change
mat1 = copy(randmat)
mat2 = copy(mat1)
mat2[:,c1] = cols[:,1]

detmat1 = AmirVMC.DetMatrix(mat1)
detmat2 = AmirVMC.DetMatrix(mat2)

# test determinant
detratio_exact = detmat2.determinant/detmat1.determinant
detratio = AmirVMC.det_ratio_1col(detmat1, cols[:,1], c1)
@test isapprox(detratio_exact, detratio)

# test the update
AmirVMC.update_detmatrix_1col!(detmat1, cols[:,1], c1, detratio)
@test isapprox(detmat1.inverse, detmat2.inverse)
@test isapprox(detmat1.matrix, detmat2.matrix)
@test isapprox(detmat1.determinant, detmat2.determinant)

### 2 column change
mat1 = copy(randmat)
mat2 = copy(mat1)
mat2[:,c1] = cols[:,1]
mat2[:,c2] = cols[:,2]

detmat1 = AmirVMC.DetMatrix(mat1)
detmat2 = AmirVMC.DetMatrix(mat2)

# test determinant
detratio_exact = detmat2.determinant/detmat1.determinant
detratio_mat = AmirVMC.det_ratio_2cols(detmat1, cols, c1, c2)
@test isapprox(detratio_exact, det(detratio_mat))

# test the update
AmirVMC.update_detmatrix_2cols!(detmat1, cols, c1, c2, detratio_mat)
@test isapprox(detmat1.inverse, detmat2.inverse)
@test isapprox(detmat1.matrix, detmat2.matrix)
@test isapprox(detmat1.determinant, detmat2.determinant)
