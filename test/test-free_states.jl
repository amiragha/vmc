# test chain
Lx = 4
model = AmirVMC.chain(Lx, :periodic, 0.0, Complex(-1.))

states1 = AmirVMC.solve_free_full(model, 4)
states2 = AmirVMC.solve_free_periodic(model, 4)

for i=1:4
    println(states1[:,i])
    println(states2[:,i])
end

# ### test periodic kagome
# Lx = 6
# num_states = div(3*Lx,2)
# model = AmirVMC.kagomestrip_LC(Lx, :periodic, -2.5, complex(1.))

# # compare momentum space and full real space solvers
# states1 = AmirVMC.solve_free_full(model, num_states)
# states2 = AmirVMC.solve_free_periodic(model, num_states)
# orthogonality = states1' * states2
# @test isapprox(orthogonality, ones(Complex128, num_states, num_states)/orthogonality[1,1])
