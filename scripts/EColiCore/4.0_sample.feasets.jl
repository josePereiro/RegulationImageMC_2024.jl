@time begin
    using RegulationImageMC_2024
    using ProjFlows
    using Random
    using DataStructures
    using Combinatorics
    using StatsBase
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.99_sim.base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# TODO: to avoid reserialization creates a serialize.flag interface
# - It should containg a symbol
# - Maybea per frame flag
#   - serializeflag!(:error) error if attempt seriallize
#   - serializeflag!(:ignore) ignore serialization calls
#   - serializeflag!(:do) do serialization if called

# --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # clear
    empty!(G)

    # meta
    script_id = "sample.feasets"
    script_ver = v"0.1.0"
    
    # done count
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
    
    # params
    refresh_dup = true
    BLOBS_PER_BATCH = 100000
    DUP_BUFF_SIZE = 10_000_000

    # all subsets
    # powerset keep it sorted
    @time for hd_bb in eachbatch(B, "hit.and.down"; sortfun = shuffle!)
        @show hd_bb.id

        # hd context control
        # TODO: control context

        islocked(hd_bb) && continue

        hd_bb_ref = blobyref(hd_bb)
        
        # Check already done
        if _check_done_count!(S, hd_bb.id, done_target; lk = true)
            @info "DONE"
            continue
        end
        
        # new BlobBatch
        ps_bb = blobbatch!(B, rbbid(script_id))
        bloblim!(ps_bb, BLOBS_PER_BATCH)

        # dup tracker
        dups_buff = _dups_tracker!(S; 
            dup_buff_size = DUP_BUFF_SIZE
        )
        @info("HASH_SET", 
            dups_buff = length(dups_buff),
        )

        try; 
            lock(hd_bb)

            for hd_vb in hd_bb
                
                # control flags
                # (if flag is missing ignore)
                get(hd_vb, "flags", "duplicate.flag", false) && continue
                get(hd_vb, "flags", "done.flag", true) || continue
                
                koset = hd_vb["cargo.koset", "koset"]
                length(koset) == 1 && continue # killing ko

                feaset0 = koset[1:end-1]
                D = length(feaset0)
                for samplei in 1:D # sample D times
                    feaset = StatsBase.sample(
                        feaset0, rand(1:D);
                        replace = false
                    )
                    
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

                        # refresh dup buf
                        # dup tracker
                        if refresh_dup 
                            dups_buff = _dups_tracker!(S; 
                                dup_buff_size = DUP_BUFF_SIZE
                            )
                        end
                        @info("HASH_SET", 
                            dups_buff = length(dups_buff),
                        )
                    end
                    
                    # store feaset
                    ## check dupplicate
                    feaset_hash = combhash(feaset)
                    if check_duplicate!(dups_buff, feaset_hash)
                        continue # if dupplicate ignore
                    end
    
                    # blob control
                    bc = blobcount(ps_bb)
                    # info
                    if iszero(mod(bc, BLOBS_PER_BATCH ÷ 10))
                        @info(script_id,
                            vblobcount = bc,
                            dups_buff_len = length(dups_buff),
                            dups_count = dups_count(dups_buff),
                            nondups_count = nondups_count(dups_buff),
                            dups_frac = dup_ratio(dups_buff)
                        )
                    end
    
                    ps_vb = rblob!(ps_bb)
                    ps_vb["cargo.feaset", "feaset"] = feaset
                
                end # for samplei
            end # for hd_b

            # done_target
            _up_done_count!(S, hd_bb.id, 1; lk = true)

            serialize!(S; lk = true)

        finally
            unlock(hd_bb)
        end

        ## capture context
        merge!(ps_bb, script_id, @litecontext)
        
        ps_bb[script_id, "hd_bb_ref"] = hd_bb_ref
        ps_bb[script_id, "script.ver"] = script_ver
        ps_bb[script_id, "src"] = read(@__FILE__, String)
        serialize!(ps_bb; lk = true)

    end # for hd_bb

    # write globals
    merge!(G, script_id, @litecontext)
    G[script_id, "src"] = read(@__FILE__, String)
    serialize!(G; lk = true)

    nothing
end


# --.-...- --. -. - -.-..- -- .-..- -. -. 