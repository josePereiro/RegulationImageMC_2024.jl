@time begin
    using RegulationImageMC_2024
    using ProjFlows
    using Random
    using DataStructures
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.99_sim.base.jl")
include("3.99_base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# NOTES:
# - downset: a downregulation set that lead to an unfeasible network
# - fea_downset: a downregulation set that lead to an  still feasible network

# - Each ensemble has a fixed enviroment

# --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # clear
    empty!(G)
    empty!(C)

    # meta
    script_id = "hit.and.down"
    script_ver = v"0.1.0"

    # local cache
    S = blobbatch!(B, 
        hashed_id(
            string("cache.", script_id), 
            script_ver
        )
    )
    S["script_id"] = script_id
    S["script_ver"] = script_ver
     
    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # PARAMETERS

    ## retrieve from net0
    netid = G["gen.net0", "netid"]::String
    down_factor = 0.0
    net0 = G["gen.net0", "net0.ref"][]
    net0_hash = hash(net0)
    objid = extras(net0, "BIOM")
    objidx = colindex(net0, objid)
    feasible_th = 1e-2
    traj_lim = 35 # max length of trajectories
    blep0 = G["gen.net0", "net0.blep0.ref"][]::LEPModel
    lb0, ub0 = lb(blep0), ub(blep0)
    M, N = size(blep0)
    BLOBS_PER_BATCH = 250
    
    ## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    
    # sim state
    opm = FBAOpModel(blep0, LP_SOLVER)
    @assert colids(opm) == colids(blep0) 
    exit_status = "ATINIT"

    # Reset (uncomment to reset)
    # foreach_batch(rm, B, script_id); rm(S); return
    
    dups_count = 0
    nondups_count = 0

    # select muatble reactions
    eblep0 = G["gen.net0", "net0.eblep0.ref"][]
    @assert colids(eblep0) == colids(blep0)
    iids_pool0 = filter(colids(eblep0, eblep0.idxi)) do id
        id == objid && return false
        startswith(id, "EX_") && return false
        return true
    end
    iidxs_pool0 = colindex(eblep0, iids_pool0)

    # for each batch
    for _batchi in 1:1
        @label BATCH_INIT
        
        @show _batchi
        
        # new BlobBatch
        hd_bb = getbatch!(B, script_id, rbbid(script_id)) do _bb
            islocked(_bb) && return false
            isfullbatch(_bb) && return false
            return true
        end
        bloblim!(hd_bb, BLOBS_PER_BATCH)

        # duplicate buffer
        dups_buff = _dups_tracker(S; 
            dup_buff_size = 1_000_000
        )
        @info("HASH_SET", 
            dups_buff = length(dups_buff),
        )

        try; lock(hd_bb)
            
            # run state
            downset = Int[]
            biomset = Float64[]
            __downfactors = ones(Float64, N)
            _downfactors = ones(Float64, N)

            # for each run
            for _runi in 1:Int(1e10)

                # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                @label RUN_INIT

                # random restarting point
                # random select a head subset of the current path
                L = length(downset) - 1
                rand_idx0 = L < 1 ? 0 : rand(0:L)
                downset = _head(downset, rand_idx0)
                biomset = _head(biomset, rand_idx0)
                _reset_idxs!(_downfactors, __downfactors, downset, 1.0)

                # runi state
                exit_status = "RUNNING"
                last_biom = isempty(biomset) ? 0.0 : last(biomset)

                # reset model
                lb!(opm, lb0 .* _downfactors)
                ub!(opm, ub0 .* _downfactors)

                # for each step
                for _stepi in 1:Int(1e10)
                    
                    @label STEP_INIT

                    # peek next rxn id
                    ridx = nothing
                    for _ridx_att in 1:500
                        ridx = rand(iidxs_pool0)
                        iszero(_downfactors[ridx]) && continue
                        _downfactors[ridx] *= down_factor
                        push!(downset, ridx)
                        break
                    end
                    
                    if length(downset) > traj_lim
                        # println("isnothing(ridx)")
                        exit_status = "ERROR.TRAJ.LIM"
                        @goto RUN_END
                    end

                    if isnothing(ridx) 
                        # println("isnothing(ridx)")
                        exit_status = "ERROR.FAILED_STEP"
                        @goto RUN_END
                    end

                    # change bounds
                    lbr = lb0[ridx] * _downfactors[ridx]
                    ubr = ub0[ridx] * _downfactors[ridx]
                    lb!(opm, ridx, lbr)
                    ub!(opm, ridx, ubr)

                    # optimize
                    try; optimize!(opm)
                        last_biom = solution(opm, objidx)
                        push!(biomset, last_biom)
                        catch err; 
                            # println("catch err;")
                            exit_status = "ERROR.ONOPTIMIZATION"
                            push!(biomset, NaN)
                            # @show(err)
                            @goto RUN_END
                    end
                    # @show ridx
                    # @show last_biom

                    # check nan
                    if !isfinite(last_biom)
                        # println("isnan(last_biom)")
                        exit_status = "OBJ.NAN.OR.INF"
                        @goto RUN_END
                    end

                    # feasibility check
                    if last_biom < feasible_th 
                        # println("last_biom < feasible_th ")
                        exit_status = "UNFEASIBLE"
                        @goto RUN_END
                    end

                    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    
                    # end trajectory computation 
                    @label RUN_END
                    
                    # continue if running
                    exit_status == "RUNNING" && continue
                    
                    # check dupplicate
                    koset_hash = combhash(downset)
                    if check_duplicate!(dups_buff, koset_hash)
                        dups_count += 1
                        @goto RUN_INIT # if dupplicate start a new run
                    end
                    nondups_count += 1

                    # Finish run (no more steps)
                    break 

                end # for _stepi 
                
                # store
                b = rblob!(hd_bb)
                b["exit_status"] = exit_status
                b["cargo.downset", "downset"] = downset
                b["cargo.biomset", "biomset"] = biomset
                
                # blob control
                bc = blobcount(hd_bb)
                # info
                if iszero(mod(bc, 100))
                    @info(script_id,
                        bc,
                        dups_count,
                        nondups_count,
                        dups_frac = dups_count / nondups_count
                    )
                end
                
                isfullbatch(hd_bb) && @goto CLOSE_BATCH

            end # for _runi

            # finish
            @label CLOSE_BATCH
            
            ## capture context
            merge!(hd_bb, "gen.net0", G[["gen.net0"]])

            merge!(hd_bb, script_id, @litecontext)
            hd_bb[script_id, "net0"] = net0
            hd_bb[script_id, "blep0"] = blep0
            hd_bb[script_id, "opm"] = opm
            hd_bb[script_id, "iidxs_pool0"] = iidxs_pool0
            hd_bb[script_id, "iids_pool0"] = iids_pool0
            
            hd_bb[script_id, "script.ver"] = script_ver
            hd_bb[script_id, "src"] = read(@__FILE__, String)

            # write
            serialize!(hd_bb)
            println("Hi")
            
            serialize!(S; lk = true)

        finally # try lock
            unlock(hd_bb)
        end

    end # for _batchi

    # write globals
    merge!(G, script_id, @litecontext)
    G[script_id, "src"] = read(@__FILE__, String)
    serialize!(G; lk = true)

    nothing
end
