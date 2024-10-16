using Statistics

# .- - -. - . .. .- .- - . - - -. - - -. -.- .--. 
function _net0_globals!(rab::raBlob, net0; 
        box_eps = 0.0, 
        box_reduce = true, 
        box_nths = NTHREADS
    )
    
    # lep0
    lep0 = lepmodel(net0)
    
    # fva_strip
    frame = hashed_id("blep0.cache.", lep0, box_eps, box_reduce, box_nths)
    @show frame
    blep0_ref = withblob!(rab, :get!, frame, "model") do
        _blep0 = fva_strip(lep0, LP_SOLVER; 
            nths = box_nths, 
            verbose = true, 
            eps = box_eps, 
            reduce = box_reduce
        )
        return _blep0
    end
    blep0 = rab[frame, "model"]

    # EchelonLEPModel
    frame = hashed_id("eblep0.cache.", blep0)
    eblep0_ref = withblob!(rab, :get!, frame, "model") do
        _eblep0 = EchelonLEPModel(blep0; verbose = true)
        return _eblep0
    end
    eblep0 = rab[frame, "model"]

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
    rab["net0"] = net0
    rab["net0.lep0"] = lep0
    rab["net0.blep0.ref"] = blep0_ref
    rab["net0.eblep0.ref"] = eblep0_ref
    # reference iders
    # - all stored index vector will refers to this order
    rab["net0.rxns"] = reactions(net0)
    rab["net0.eblep0.idxi"] = eblep0.idxi
    rab["net0.biom0"] = biom0

    return nothing
end

# .- - -. - . .. .- .- - . - - -. - - -. -.- .--. 
nothing