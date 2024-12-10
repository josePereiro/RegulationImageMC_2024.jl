using Statistics

# .- - -. - . .. .- .- - . - - -. - - -. -.- .--. 
function _net0_models!(net0; 
        script_id = "gen.net0",
        box_eps = 0.0, 
        box_reduce = true, 
        box_nths = NTHREADS
    )
    
    # net0
    frame = hashed_id("net0.cache.", net0)
    net0_ref = blobyio!(C, frame, "model", :getser!) do
        return net0
    end

    # lep0
    frame = hashed_id("lep0.cache.", net0)
    lep0_ref = blobyio!(C, frame, "model", :getser!) do
        return lepmodel(net0)
    end
    lep0 = C[lep0_ref]
    
    # fva_strip
    biom_id = extras(net0, "BIOM")
    frame = hashed_id("blep0.cache.", lep0, box_eps, box_reduce)
    blep0_ref = blobyio!(C, frame, "model", :getser!) do
        _blep0 = fva_strip(lep0, LP_SOLVER; 
            nths = box_nths, 
            verbose = true, 
            eps = box_eps, 
            reduce = box_reduce
        )
        linear_weights!(_blep0, biom_id, 1.0)
        return _blep0
    end
    blep0 = C[blep0_ref]

    # EchelonLEPModel
    frame = hashed_id("eblep0.cache.", blep0)
    eblep0_ref = blobyio!(C, frame, "model", :getser!) do
        _eblep0 = EchelonLEPModel(blep0; verbose = true)
        linear_weights!(net0, biom_id, 1.0)
        return _eblep0
    end
    eblep0 = C[eblep0_ref]

    # Test FBA
    bioms = Float64[]
    for m in [net0, lep0, blep0, eblep0]
        opm = FBAOpModel(m, LP_SOLVER)
        optimize!(opm)
        biom = objective_value(opm)
        push!(bioms, biom)
        @show biom
    end
    # biom0 = mean(bioms)
    biom0 = maximum(bioms)
    # @assert all(isapprox.(bioms[1], bioms; atol = 1e-4))

    # store globals
    G[script_id, "net0.ref"] = net0_ref
    G[script_id, "net0.lep0.ref"] = lep0_ref
    G[script_id, "net0.blep0.ref"] = blep0_ref
    G[script_id, "net0.eblep0.ref"] = eblep0_ref
    # reference iders
    # - all stored index vector will refers to this order
    G[script_id, "net0.rxns"] = reactions(net0)
    G[script_id, "net0.eblep0.idxi"] = eblep0.idxi
    G[script_id, "net0.biom0"] = biom0

    return nothing
end

# .- - -. - . .. .- .- - . - - -. - - -. -.- .--. 
nothing