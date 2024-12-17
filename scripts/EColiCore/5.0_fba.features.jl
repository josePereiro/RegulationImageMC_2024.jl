@time begin
    using RegulationImageMC_2024
    using ProjFlows
    using Random
    using DataStructures
    using Combinatorics
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.99_sim.base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # clear
    empty!(G)

    # meta
    script_id = "fba.feasures"
    script_ver = v"0.1.0"
    ctx_hash = combhash(script_id, script_ver)

    # local cache
    S = blobbatch!(B, 
        hashed_id(
            string("cache.", script_id), 
            script_ver
        )
    )
    S["script_id"] = script_id
    S["script_ver"] = script_ver

    # Reset (uncomment to reset)
    # foreach_batch(rm, B, script_id); rm(S); return

    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # PARAMETERS

    netid = G["gen.net0", "netid"]
    net0 = G["gen.net0", "net0.ref"][]::MetNet
    net0_hash = hash(net0)
    blep0 = G["gen.net0", "net0.blep0.ref"][]::LEPModel
    lb0, ub0 = lb(blep0), ub(blep0)
    M, N = size(blep0)

    ## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    # up context
    ctx_hash = combhash(ctx_hash, net0_hash)

    # state
    opm = FBAOpModel(blep0, LP_SOLVER)
    exit_status = "ATINIT"

    # all subsets
    # powerset keep it sorted
    @time for ps_bb in eachbatch(B, "gen.powersets")
        @show ps_bb.id

        # context control
        # TODO: do context control
        
        # ctx control

        islocked(ps_bb) && continue

        ff_bb = nothing
        try; 
            lock(ps_bb)

            ps_bb_ref = blobyref(ps_bb)
            
            # Check already done
            version_hist = get!(ps_bb, script_id, "version.hist", Set()) 
            script_ver ∈ version_hist && continue

            # new BlobBatch
            ff_bb = blobbatch!(B, rbbid(script_id))
            bloblim!(ff_bb, bloblim(ps_bb))

            lock(ff_bb)

            # @show ps_bb
            for (bi, ps_b) in enumerate(ps_bb)
                if iszero(mod(bi, 100))
                    @show (bi, ps_b.uuid)
                end

                # TODO: add get(vb, key..., dflt) interface also for frames
                # get(ps_b, "duplicate.flag", false) && continue

                feaset = ps_b["cargo.feaset", "feaset"]
                D = length(feaset)
                # @show D
    
                # apply feaset
                ## reset model
                lb!(opm, lb0)
                lb!(opm, feaset, 0)
                @assert all(iszero, lb(opm, feaset))
                ub!(opm, ub0)
                ub!(opm, feaset, 0)
                @assert all(iszero, ub(opm, feaset))
    
                # Test FBA
                sol = Float64[]
                try
                    optimize!(opm)
                    sol = solution(opm)
                catch err
                    # @show err
                end 
    
                ff_vb = blob!(ff_bb, ps_b.uuid) # sync batches
                ff_vb["cargo.fba", "sol"] = sol
    
            end # for ps_b
    
            # write ps_bb 
            # (keep track of ss versions)
            push!(version_hist, script_ver)
            serialize!(ps_bb, script_id)
    
            # write ff_bb
            ## ctx
            merge!(ff_bb, script_id, @litecontext)

            ff_bb[script_id, "ps_bb_ref"] = ps_bb_ref
            ff_bb[script_id, "src"] = read(@__FILE__, String)
            
            serialize!(ff_bb)

        finally;
            unlock(ps_bb)
            isnothing(ff_bb) || unlock(ff_bb)
        end

    end # for hd_bb

    # write globals
    merge!(G, script_id, @litecontext)
    G[script_id, "src"] = read(@__FILE__, String)
    serialize!(G; lk = true)

    nothing
end


# --.-...- --. -. - -.-..- -- .-..- -. -. 