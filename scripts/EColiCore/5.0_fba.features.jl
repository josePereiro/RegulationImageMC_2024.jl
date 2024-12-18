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
let
    # clear
    empty!(G)

    # meta
    script_id = "fba.feasures"
    script_ver = v"0.1.0"
    done_count = 1

    # local cache
    global S = blobbatch!(B, 
        hashed_id(
            string("cache.", script_id), 
            script_ver
        )
    )
    S["script_id"] = script_id
    S["script_ver"] = script_ver

    # all subsets
    # powerset keep it sorted
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
            # TODO: Fix this
            # - The problem is that you need to comunicate with the other processes
            # - similar to dups_buff
            # - TODO: Make a secure interface for data exchange 
            #   - Basically is locking and merging your version with the one on disk...
            global done_reg = get!(S, script_id, "hist", Dict())
            if get(done_reg, ps_bb.id, -1) === done_count 
                @info("DONE")
                continue
            end

            # done_count
            # get!(done_reg, ps_bb.id, 0)
            # done_reg[ps_bb.id] += 1
            # serialize!(S)

        finally;
            unlock(ps_bb)
            isnothing(ff_bb) || unlock(ff_bb)
        end

    end # for ps_bb


    nothing
end


# --.-...- --. -. - -.-..- -- .-..- -. -. 


## --.-...- --. -. - -.-..- -- .-..- -. -. 
## --.-...- --. -. - -.-..- -- .-..- -. -. 
## --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # clear
    empty!(G)

    # meta
    script_id = "fba.feasures"
    script_ver = v"0.1.0"
    ctx_hash = combhash(script_id, script_ver)
    done_count = 1

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
            done_reg = get!(S, script_id, "hist", Dict())
            get(done_reg, ps_bb.id, -1) === done_count && continue
            return done_reg

            # new BlobBatch
            ff_bb = blobbatch!(B, hashed_id(script_id, ps_bb.id))
            bloblim!(ff_bb, bloblim(ps_bb))

            lock(ff_bb)

            # @show ps_bb
            for (bi, ps_b) in enumerate(ps_bb)
                # info
                iszero(mod(bi, 100)) &&
                    println(
                        "bi: ", bi, "/", ps_bc,
                        " [",
                        " pid: ", getpid(), 
                        " thid: ", threadid(),
                        "]"
                    )

                # control duplicates
                get(ps_b, "flags", "duplicate.flag", false)  && continue

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
    
            # write ff_bb
            ## ctx
            merge!(ff_bb, script_id, @litecontext)

            ff_bb[script_id, "ps_bb_ref"] = ps_bb_ref
            ff_bb[script_id, "src"] = read(@__FILE__, String)
            
            serialize!(ff_bb)

            # done_count
            get!(done_reg, ps_bb.id, 0)
            done_reg[ps_bb.id] += 1
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