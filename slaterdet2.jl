type SlaterDet2
    matrix_up :: Matrix{Complex128}
    matrix_dn :: Matrix{Complex128}

    det_up :: Complex128
    det_dn :: Complex128

    function SlaterDet2(mat_up::Matrix{Complex128},
                        mat_dn::Matrix{Complex128})
        lu_up = lufact(mat_up)
        lu_dn = lufact(mat_dn)
        det_up = det(lu_up)
        det_dn = det(lu_dn)
        return new(mat_up, mat_dn, det_up, det_dn)
    end
end

function slater_swap_column(slater::SlaterDet2,
                            up_idx::Int,
                            dn_idx::Int)

    mat_up = copy(slater.matrix_up)
    mat_dn = copy(slater.matrix_dn)
    col_up = mat_up[:, up_idx]
    col_dn = mat_dn[:, dn_idx]

    mat_up[:, up_idx] = col_dn
    mat_dn[:, dn_idx] = col_up

    return SlaterDet2(mat_up, mat_dn)
end
