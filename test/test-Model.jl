# test chain
function test_model_chain(Lx::Int, boundary::Symbol)
    model = VMC.chain(Lx, boundary, 0., complex(1.))

    println(boundary)
    for i=1:length(model.lattice.sites)
        println(model.lattice.sites[i])
    end
    for i=1:length(model.lattice.edges)
        println(model.lattice.edges[i])
    end
end

# test kagome
function test_model_kagome(Lx::Int, boundary::Symbol)
    model = VMC.kagomestrip_LC(Lx, boundary, -2.5, complex(1.))

    println(boundary)
    for i=1:length(model.lattice.sites)
        println(model.lattice.sites[i])
    end
    for i=1:length(model.lattice.edges)
        println(model.lattice.edges[i])
    end
end

test_model_kagome(3, :periodic)
test_model_kagome(3, :opensym)
test_model_kagome(3, :open)

test_model_chain(4,:periodic)
test_model_chain(4,:open)
