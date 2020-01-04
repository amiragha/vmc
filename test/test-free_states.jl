# test chain
Lx = 4
model = VMC.chain(Lx, :periodic, 0.0, Complex(-1.))

states1 = VMC.solve_free_full(model, 4)
states2 = VMC.solve_free_periodic(model, 4)

for i=1:4
    println(states1[:,i])
    println(states2[:,i])
end

# ### test periodic kagome
# Lx = 6
# num_states = div(3*Lx,2)
# model = VMC.kagomestrip_LC(Lx, :periodic, -2.5, complex(1.))

# # compare momentum space and full real space solvers
# states1 = VMC.solve_free_full(model, num_states)
# states2 = VMC.solve_free_periodic(model, num_states)
# orthogonality = states1' * states2
# @test isapprox(orthogonality, ones(ComplexF64, num_states, num_states)/orthogonality[1,1])
