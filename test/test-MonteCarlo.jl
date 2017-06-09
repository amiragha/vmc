# test small chain

# model = AmirVMC.chain(4, :periodic, 0., complex(1.0))
# AmirVMC.runVMC(model, 2, true)

# test VMC for chain

# model = AmirVMC.chain(32, :periodic, 0., complex(1.0))
# AmirVMC.runVMC(model, 10000)

# test VMC for kagome strip

model = AmirVMC.kagomestrip_LC(40, :open, -2.5, complex(1.0))
AmirVMC.runVMC(model, 800000)
