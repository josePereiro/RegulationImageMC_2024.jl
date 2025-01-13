@time begin
    using RegulationImageMC_2024
    using ProjFlows
    using Random
    using Base.Threads
    using DataStructures
    using Combinatorics
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.99_sim.base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
#MARK: # fba.feasures
let
    # clear
    empty!(G)

    # meta
    script_id = "fba.feasures"
    script_ver = v"0.1.0"
    ctx_hash = combhash(script_id, script_ver)
    done_target = 1

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
    #MARK: ## opm
    opm = FBAOpModel(blep0, QUAD_LP_SOLVER)
    fix_delta = 1e-3

    # all subsets
    # powerset keep it sorted
    #MARK: ## for ps_bb
    @time for ps_bb in eachbatch(B, "sample.feasets")
        @show ps_bb.id
        ps_bc = blobcount(ps_bb)

        # context control
        # TODO: do context control
        
        # control
        islocked(ps_bb) && continue

        ff_bb = nothing
        try; 
            lock(ps_bb)

            # TODO: Test there are one fba batch for each sampled.feaset batch
            ps_bb_ref = blobyref(ps_bb)
            
            # Check already done
            if _check_done_count!(S, ps_bb.id, done_target; lk = true)
                @info "DONE"
                continue
            end

            # new BlobBatch
            ff_bb = blobbatch!(B, hashed_id(script_id, ps_bb.id))
            # rm(ff_bb); empty!(ff_bb)
            bloblim!(ff_bb, bloblim(ps_bb))

            lock(ff_bb)

            # @show ps_bb
            #MARK: ### for (bi, ps_b)
            for (bi, ps_b) in enumerate(ps_bb)
                # info
                iszero(mod(bi, 100)) &&
                    println(
                        "bi: ", bi, "/", ps_bc,
                        " [",
                        "pid: ", getpid(), 
                        " thid: ", threadid(),
                        "]"
                    )

                # control flags
                #MARK: #### control flags
                get(ps_b, "flags", "duplicate.flag", false)  && continue

                feaset = ps_b["cargo.feaset", "feaset"]
                D = length(feaset)
    
                # apply feaset
                #MARK: #### apply feaset
                ## reset model
                op_status = "NEW"
                lb!(opm, lb0)
                lb!(opm, feaset, 0)
                @assert all(iszero, lb(opm, feaset))
                ub!(opm, ub0)
                ub!(opm, feaset, 0)
                @assert all(iszero, ub(opm, feaset))

                # new blob
                ff_vb = blob!(ff_bb, ps_b.uuid) # sync batches
    
                #MARK: #### FBA max.biom.fix
                sol = Float64[]
                try
                    # biom
                    id = extras(net0, "BIOM")
                    set_linear_obj!(opm, id, MAX_SENSE)
                    optimize!(opm)
                    
                    obj = solution(opm, id)
                    l, u = bounds(opm, id)
                    lb!(opm, id, max(obj - (obj * fix_delta), l))
                    ub!(opm, id, min(obj + (obj * fix_delta), u))
                    
                    sol = solution(opm)
                    op_status = "KO"

                catch err
                    op_status = "ERR.max.fba"
                end 
    
                ff_vb["cargo.fba.max.biom", "sol"] = sol
                ff_vb["cargo.fba.max.biom", "status"] = op_status

                #MARK: #### FBA min.glc.fix
                sol = Float64[]
                try
                    op_status != "OK" && error(op_status)

                    # glc
                    id = extras(net0, "EX_GLC")
                    set_linear_obj!(opm, id, MIN_SENSE)
                    optimize!(opm)
                    
                    obj = solution(opm, id)
                    l, u = bounds(opm, id)
                    lb!(opm, id, max(obj - (obj * fix_delta), l))
                    ub!(opm, id, min(obj + (obj * fix_delta), u))
                    
                    sol = solution(opm)
                    op_status = "KO"

                catch err
                    op_status = "ERR.min.glc"
                end 
    
                ff_vb["cargo.fba.min.glc", "sol"] = sol
                ff_vb["cargo.fba.min.glc", "status"] = op_status

                #MARK: #### FBA min.v2
                sol = Float64[]
                try
                    op_status != "OK" && error(op_status)

                    # Min v2
                    set_v2_obj!(opm, MIN_SENSE)
                    optimize!(opm)

                    sol = solution(opm)
                    op_status = "KO"

                catch err
                    op_status = "ERR.min.v2"
                end 
    
                ff_vb["cargo.fba.min.v2", "sol"] = sol
                ff_vb["cargo.fba.min.v2", "status"] = op_status

                #MARK: #### FBA max.v2
                sol = Float64[]
                try
                    op_status != "OK" && error(op_status)

                    # Max v2
                    set_v2_obj!(opm, MAX_SENSE)
                    optimize!(opm)

                    sol = solution(opm)
                    op_status = "KO"

                catch err
                    op_status = "ERR.max.v2"
                end 
    
                ff_vb["cargo.fba.max.v2", "sol"] = sol
                ff_vb["cargo.fba.max.v2", "status"] = op_status

    
            end # for ps_b
    
            # write ff_bb
            ## ctx
            merge!(ff_bb, script_id, @litecontext)

            ff_bb[script_id, "ps_bb_ref"] = ps_bb_ref
            ff_bb[script_id, "src"] = read(@__FILE__, String)
            
            serialize!(ff_bb)

            # done_target
            _up_done_count!(S, ps_bb.id, 1; lk = true)

            serialize!(S)

        finally;
            unlock(ps_bb)
            isnothing(ff_bb) || unlock(ff_bb)
        end

    end # for ps_bb

    # write globals
    merge!(G, script_id, @litecontext)
    G[script_id, "src"] = read(@__FILE__, String)
    serialize!(G; lk = true)

    nothing
end


# --.-...- --. -. - -.-..- -- .-..- -. -. 