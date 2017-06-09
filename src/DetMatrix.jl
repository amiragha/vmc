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
det_ratio_1col(dmat::DetMatrix, new_col:: Vector{Complex128}, col_index :: Int)

Making use of the Matrix Determinant Lemma (MDL), finds the ratio of
the new determinant to old one given by the change of `col_index`
column to `new_col` in the `dmat`.

returns determinant ratio.
"""
function det_ratio_1col(dmat      :: DetMatrix,
                        new_col   :: Vector{Complex128},
                        col_index :: Int)
    # row syntax return a column(Vector)!
    # need a dot product which doesn't conjugate the first vector!
    return sum(dmat.inverse[col_index,:] .* new_col)
end

"""
exact_det_ratio_1col(dmat::DetMatrix, new_col:: Vector{Complex128}, col_index :: Int)

exact version of the determinant ratio for 1 column change

returns determinant ratio.
"""
function exact_det_ratio_1col(dmat      :: DetMatrix,
                              new_col   :: Vector{Complex128},
                              col_index :: Int)
    newmat = copy(dmat.matrix)
    newmat[:,col_index] = new_col
    return det(newmat)/dmat.determinant
end

"""
update_detmatrix_1col!(dmat:: DetMatrix, new_col   :: Vector{Complex128}, col_index :: Int, detratio ::  Complex128)

Making use of the Sherman-Morrison formula, updates the matrix,
determinant and inverse of the new DetMatrix given by the change of
`col_index` column to `new_col`.

updates dmat.
"""
function update_detmatrix_1col!(dmat      :: DetMatrix,
                                new_col   :: Vector{Complex128},
                                col_index :: Int,
                                detratio ::  Complex128)
    # update matrix
    dmat.matrix[:,col_index] = new_col

    # update determinant
    dmat.determinant *= detratio

    # update inverse, note: column is an auxiliary vector
    # TODO: elegance!
    column = dmat.inverse * new_col
    column[col_index] -= 1
    dmat.inverse -= (column .* transpose(dmat.inverse[col_index,:]))/detratio
    return nothing
end

"""
exact_update_detmatrix_1col!(dmat:: DetMatrix, new_col   :: Vector{Complex128}, col_index :: Int, detratio ::  Complex128)

exact version that updates the matrix, determinant and inverse of the
new DetMatrix given by the change of `col_index` column to `new_col`.

updates dmat.
"""
function exact_update_detmatrix_1col!(dmat      :: DetMatrix,
                                      new_col   :: Vector{Complex128},
                                      col_index :: Int,
                                      detratio ::  Complex128)
    newmat = copy(dmat.matrix)
    newmat[:, col_index] = new_col

    dmat.matrix = newmat
    dmat.inverse = inv(newmat)
    dmat.determinant = det(newmat)
    return nothing
end
