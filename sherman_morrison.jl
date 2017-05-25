# using the Matrix Determinant Lemma (MDL)
function new_determinant_MDL(mat::Matrix{Complex128},
                         inverse::Matrix{Complex128},
                         determinant::Complex128,
                         new_col :: Vector{Complex128},
                         col_index :: Int)
    # I don't know why row slice syntax return a column?!
    # and there doesn't exist a dot product which doesn't conjugate
    return determinant * (transpose(inverse[col_index,:]) * new_col)[1]
end

function sherman_morrison_1col(mat::Matrix{Complex128},
                               inverse::Matrix{Complex128},
                               determinant::Complex128,
                               new_col :: Vector{Complex128},
                               col_index :: Int)
    # calculate the new determinant using Matrix Determinant Lemma
    # I don't know why row slice syntax return a column?!
    # and there doesn't exist a dot product which doesn't conjugate
    det_ratio = (transpose(inverse[col_index,:]) * new_col)[1]

    # making a vector of zeros with kth element 1
    one_k_col = zeros(Complex128, length(new_col))
    delta_k_col[col_index] = 1

    # the Sherman_morrison formula for changing "k"th colum to "col"
    return inverse - kron((inverse * new_col - one_k_col),transpose(inverse[col_index,:]))/det_ratio

end
