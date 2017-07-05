# DetMatrix type to represent determinantal wavefunctions
type DetMatrix
    matrix      :: Matrix{Complex128}
    determinant :: Complex128
    inverse     :: Matrix{Complex128}
end

# constructor with just a matrix
function DetMatrix(mat::Matrix{Complex128})
    lu = lufact(mat)
    return DetMatrix(mat, det(mat), inv(mat))
end

"""
det_ratio_1col(dmat::DetMatrix, col:: Vector{Complex128}, c :: Int)

Making use of the Matrix Determinant Lemma (MDL), finds the ratio of
the new determinant to old one given by the change of `c`
column to `col` in the `dmat`.

returns determinant ratio.
"""
function det_ratio_1col(dmat      :: DetMatrix,
                        col   :: Vector{Complex128},
                        c :: Int)
    # row syntax return a column(Vector)!
    # need a dot product which doesn't conjugate the first vector!
    return (dmat.inverse[[c],:] * col)[1]
end

"""
function det_ratio_2col(dmat:: DetMatrix, cols:: Matrix{Complex128}, c1:: Int, c2:: Int)

Making use of the Matrix Determinant Lemma (MDL), finds the matrix 2x2
that its determinant is ratio of the new determinant to old one given
by the change of `c1, c2` column to the two columns in `cols`
respectively.

returns determinant ratio matrix (2x2)
"""
function det_ratio_2cols(dmat      :: DetMatrix,
                         cols      :: Matrix{Complex128},
                         c1        :: Int,
                         c2        :: Int)
    # row syntax return a column(Vector)!
    # need a dot product which doesn't conjugate the first vector!\
    ## TODO: find a vectorize way to do this?
    detratio_mat = zeros(Complex128, 2, 2)
    detratio_mat[1, 1] = (dmat.inverse[[c1],:] * cols[:,1])[1]
    detratio_mat[1, 2] = (dmat.inverse[[c1],:] * cols[:,2])[1]
    detratio_mat[2, 1] = (dmat.inverse[[c2],:] * cols[:,1])[1]
    detratio_mat[2, 2] = (dmat.inverse[[c2],:] * cols[:,2])[1]
    return detratio_mat
end

"""
update_detmatrix_1col!(dmat:: DetMatrix, new_col   :: Vector{Complex128}, col_index :: Int, detratio ::  Complex128)

Making use of the Sherman-Morrison formula, updates the matrix,
determinant and inverse of the new DetMatrix given by the change of
`col_index` column to `new_col`.

updates dmat.
"""
function update_detmatrix_1col!(dmat      :: DetMatrix,
                                col   :: Vector{Complex128},
                                c :: Int,
                                detratio ::  Complex128)
    # update matrix
    dmat.matrix[:,c] = col

    # update determinant
    dmat.determinant *= detratio

    # update inverse, note: column is an auxiliary vector
    # TODO: elegance!
    column = dmat.inverse * col
    column[c] -= 1

    dmat.inverse -= (column * dmat.inverse[[c],:])/detratio
    return nothing
end

"""
update_detmatrix_2cols!(dmat::DetMatrix, cols:: Matrix{Complex128}, c1:: Int, c1:: Int, detratio ::  Matrix{Complex128})

Making use of the generalized Sherman-Morrison formula (Woodbury?!),
updates the matrix, determinant and inverse of the new DetMatrix given
by the change of `c1,c2` columns to the two columns in `cols` respectively.

updates dmat.
"""
function update_detmatrix_2cols!(dmat      :: DetMatrix,
                                 cols      :: Matrix{Complex128},
                                 c1        :: Int,
                                 c2        :: Int,
                                 detratio_mat  ::  Matrix{Complex128})
    # update matrix
    dmat.matrix[:,c1] = cols[:, 1]
    dmat.matrix[:,c2] = cols[:, 2]

    # update determinant
    dmat.determinant *= det(detratio_mat)

    # update inverse, note: columns is an auxiliary 2 column matrix
    # TODO: elegance!
    columns = dmat.inverse * cols
    columns[c1, 1] -= 1
    columns[c2, 2] -= 1

    dmat.inverse -= columns * inv(detratio_mat) * dmat.inverse[[c1,c2], :]
    return nothing
end

"""
check_and_update_detmatrix!(dmat::DetMatrix)

check the inverse and determinant in the detmat and if they are off,
replace them with exact values. TODO: What should the tolerance be here?!
"""
function check_and_update_detmatrix!(dmat::DetMatrix)
    if isapprox(dmat.matrix * dmat.inverse,
                eye(Complex128, size(dmat.matrix)[1]), rtol=1.e-13)
        return nothing
    else
        print("o")
        dmat.inverse = inv(dmat.matrix)
        return nothing
    end
end
