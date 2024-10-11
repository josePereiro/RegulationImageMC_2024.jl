using Statistics

# .- - -. - . .. .- .- - . - - -. - - -. -.- .--. 
function _net0_globals!(rab::raBlob, net0)
    
    # lep0
    lep0 = lepmodel(net0)
    
    # box
    blep0 = box(lep0, LP_SOLVER; nths = NTHREADS, verbose = true)

    # EchelonLEPModel
    eblep0 = EchelonLEPModel(blep0; verbose = true)

    # Test FBA
    bioms = Float64[]
    for m in [net0, lep0, blep0, eblep0]
        opm = FBAOpModel(m, LP_SOLVER)
        optimize!(opm)
        biom = objective_value(opm)
        push!(bioms, biom)
        @show biom
    end
    biom0 = mean(bioms)
    @assert all(isapprox.(bioms[1], bioms; atol = 1e-4))

    # store globals
    rab["net0"] = net0
    rab["net0.lep0"] = lep0
    rab["net0.blep0"] = blep0
    rab["net0.eblep0"] = eblep0
    # reference iders
    # - all stored index vector will refers to this order
    rab["net0.rxns"] = reactions(net0)
    rab["net0.eblep0.idxi"] = eblep0.idxi
    rab["net0.biom0"] = biom0

    return nothing
end

# .- - -. - . .. .- .- - . - - -. - - -. -.- .--. 
nothing