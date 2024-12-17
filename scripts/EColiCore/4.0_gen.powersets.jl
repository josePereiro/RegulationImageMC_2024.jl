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
# TODO: I need a way to know if a power set is already done
# - A two way graph of relationships
#   - hb_bb <-> ps_bb
# - It can be done by naming the ps_bb accordantly
#   - But again, DO NOT USE FILE NAMES FOR CONTEXT RESOLUTION

# --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # clear
    empty!(G)

    # meta
    script_id = "gen.powersets"
    script_ver = v"0.1.0"

    BLOBS_PER_BATCH = 5000

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

    # state
    dups_count = 0
    nondups_count = 0

    # compat
    hd_script_ver = v"0.1.0"
    exch_ids = ["EX_glc__D_e"]

    # all subsets
    # powerset keep it sorted
    @time for hd_bb in eachbatch(B, "hit.and.down")
        @show hd_bb.id

        # hd context control
        _hd_script_ver = get(hd_bb, "hit.and.down", "script.ver", "") 
        _hd_script_ver == hd_script_ver || continue
        _exch_ids = get(hd_bb, "gen.net0", "exch_ids", "") 
        _exch_ids == exch_ids || continue

        islocked(hd_bb) && continue

        hd_bb_ref = blobyref(hd_bb)
        
        # Check already done
        version_hist = get!(hd_bb, script_id, "version.hist", Set()) 
        script_ver ∈ version_hist && continue

        # new BlobBatch
        ps_bb = blobbatch!(B, rbbid(script_id))
        bloblim!(ps_bb, BLOBS_PER_BATCH)

        # dup tracker
        dups_buff = _dups_tracker(S; 
            dup_buff_size = 2_000_000
        )
        @info("HASH_SET", 
            dups_buff = length(dups_buff),
        )

        try; 
            lock(hd_bb)

            for hd_vb in hd_bb
                get(hd_vb, "duplicate.flag", false) && continue
                
                downset0 = hd_vb["cargo.downset", "downset"]
                D = length(downset0)
                # @show D
                    
                    # suppersets
                    sub_downsets = powerset(downset0, 1, D)
                    for feaset in sub_downsets
                        
                        if isfullbatch(ps_bb)
                            ## capture context
                            merge!(ps_bb, script_id, @litecontext)
                            
                            ps_bb[script_id, "hd_bb_ref"] = hd_bb_ref
                            ps_bb[script_id, "script.ver"] = script_ver
                            ps_bb[script_id, "src"] = read(@__FILE__, String)

                            # write
                            serialize!(ps_bb)
        
                            # new ps_bb
                            ps_bb = blobbatch!(B, rbbid(script_id))
                            bloblim!(ps_bb, BLOBS_PER_BATCH)

                            # dup tracker
                            dups_buff = _dups_tracker(S; 
                                dup_buff_size = 2_000_000
                            )
                            @info("HASH_SET", 
                                dups_buff = length(dups_buff),
                            )
                        end
                        
                        # store feaset
                        ## check dupplicate
                        feaset_hash = combhash(feaset)
                        if check_duplicate!(dups_buff, feaset_hash)
                            dups_count += 1 
                            continue # if dupplicate ignore
                        end
                        nondups_count += 1
        
                        # blob control
                        bc = blobcount(ps_bb)
                        # info
                        if iszero(mod(bc, BLOBS_PER_BATCH ÷ 10))
                            @info(script_id,
                                bc,
                                dups_buff_len = length(dups_buff),
                                dups_count,
                                nondups_count,
                                dups_frac = dups_count / nondups_count
                            )
                        end
        
                        ps_vb = rblob!(ps_bb)
                        ps_vb["cargo.feaset", "feaset"] = feaset
                    end
            end # for hd_b

        finally
            unlock(hd_bb)
        end

        # write hd_bb (keep track of ss versions)
        push!(version_hist, script_ver)
        serialize!(hd_bb, script_id; lk = true)

        ## capture context
        merge!(ps_bb, script_id, @litecontext)
        
        ps_bb[script_id, "hd_bb_ref"] = hd_bb_ref
        ps_bb[script_id, "script.ver"] = script_ver
        ps_bb[script_id, "src"] = read(@__FILE__, String)

        # write
        serialize!(ps_bb; lk = true)
        serialize!(S; lk = true)

    end # for hd_bb

    # write globals
    merge!(G, script_id, @litecontext)
    G[script_id, "src"] = read(@__FILE__, String)
    serialize!(G; lk = true)

    nothing
end


# --.-...- --. -. - -.-..- -- .-..- -. -. 