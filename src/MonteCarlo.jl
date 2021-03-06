mutable struct VMCgutzwiller
    model :: Model

    wavefunction :: GutzwillerSlater

    total_steps :: Int64

    num_proposed :: Int64
    num_accepted :: Int64

    measurement :: Measurement
end

"""
    report(simulation)

reports the results and saves the measurements for VMCgutzwiller
`simulation`.

"""
function report(simulation::VMCgutzwiller)
    println("Finished with ", simulation.total_steps,
            " of VMC for Model: ", simulation.model.name)
    println("Number of proposed moves = ", simulation.num_proposed)
    println("Number of accepted moves = ", simulation.num_accepted)
    report(simulation.measurement, simulation.model.name)
    return nothing
end

"""
    runVMC(model, total_steps)

runs the Variational Monte-Carlo simulation until it moves by
`total_steps` and measures the measurements at each accepted step. It
finally reports the results.

"""
function runVMC(model::Model,
                total_steps::Int64,
                debug::Bool=false)

    num_sites = length(model.lattice.sites)

    if model.lattice.boundary == :periodic
        #states = solve_free_full(model, div(num_sites,2))
        #states = solve_free_periodic(model, div(num_sites,2))
        states = solve_free_periodic(model, div(num_sites,2), :APBC)
    else
        states = solve_free_full(model, div(num_sites,2))
    end

    wavefunction = random_gutzwiller_half(states)

    half_correlations = Measurement{Float64}(:correlation,
                                             total_steps, num_sites)

    sim = VMCgutzwiller(model, wavefunction,
                        total_steps, 0, 0,
                        half_correlations)

    while sim.num_accepted < sim.total_steps
        sim.num_proposed += 1
        if propose_step!(sim.wavefunction)
            sim.num_accepted += 1
            measure!(sim.measurement, sim.wavefunction)

            # display progress every 1000 accepted steps
            if sim.num_accepted % 1000 == 0
                check_and_update_gutzwiller!(sim.wavefunction)
                print("*")
            end

            if debug
                println(wavefunction.configuration)
                println(wavefunction.upslater.matrix)
                println(wavefunction.dnslater.matrix)
                println(sim.measurement.accdata)
            end
        end
    end
    println()
    report(sim)
end
