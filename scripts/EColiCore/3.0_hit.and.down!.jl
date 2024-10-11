@time begin
    using RegulationImageMC_2024
    using Gurobi
    using ProjFlows
    using Random
    using DataStructures
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.0_sim.base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# NOTES:
# - ko_downset: a downregulation set that lead to an unfeasible network
# - fea_downset: a downregulation set that lead to an  still feasible network

# - Each ensemble has a fixed enviromen

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# DOING: to reduce preassure on memory, implement a Circular buffer to check 
# for very common duplications. We just ignore and do not record them on disk...

## --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # meta
    script_version = v"0.1.0"
    
    # hyper-params
    
    # globals blobs
    sim_globs = blob(B, "sim.globals")
    
    net0_globs_id = sim_globs["net0.globals.id"]::String
    net0_globs = blob(B, net0_globs_id)
    
    netid = net0_globs["net0.netid"]::String
    @show netid
    
    hnd_id = "hit.and.down" # TO SYNC
    sim_globs[string(hnd_id, ".id")] = hnd_id
    hnd_globs_id = "$(netid).$(hnd_id).globals"
    hnd_globs = blob!(B, hnd_globs_id)


    # Reset (uncomment to reset)
    # rm(hnd_globs) 
    # foreach_batch(rm, B, hnd_id)

    # duplicates
    # TODO make write/read frequent and locked
    # Make it an struct wihch creates a global (rablob) in a given Bloberia
    dups_buff = get!(hnd_globs, "local.kosets.duplicates.buffer") do
        HashTracker(50_000)
    end
    dups_count = 0
    nondups_count = 0
    
    # params
    down_factor = 0.0
    net0 = net0_globs["net0"]::MetNet
    objid = extras(net0, "BIOM")
    objidx = colindex(net0, objid)
    feasible_th = 1e-2
    traj_lim = 35 # max length of trajectories
    blep0 = net0_globs["net0.blep0"]::LEPModel
    lb0, ub0 = lb(blep0), ub(blep0)
    M, N = size(blep0)
    opm = FBAOpModel(blep0, LP_SOLVER)

    # select muatble reactions
    idxi = net0_globs["net0.eblep0.idxi"] # independent variables
    iids_pool0 = filter(colids(blep0, idxi)) do id
        id == objid && return false
        startswith(id, "EX_") && return false
        return true
    end
    iidxs_pool0 = colindex(blep0, iids_pool0)

    # for each batch
    for _batchi in 1:25
        
        # new BlobBatch
        hd_bb = headbatch!(B, hnd_id)
        if islocked(hd_bb) # ignore locked 
            hd_bb = blobbatch!(B, hnd_id)
        end
        setmeta!(hd_bb, "blobs.lim", 1000)
        lock(hd_bb)
        
        downset = Int[]
        biomset = Float64[]
        __downfactors = ones(Float64, N)
        _downfactors = ones(Float64, N)

        # for each run
        for _seti in 1:Int(1e10)

            # random restarting point
            @label RANDOM_RESTART_POINT
            L = length(downset)
            rand_idx0 = rand(0:L)
            downset = iszero(rand_idx0) ? Int[] : downset[1:rand_idx0]
            biomset = iszero(rand_idx0) ? Float64[] : biomset[1:rand_idx0]
            __downfactors .= _downfactors
            _downfactors .= 1.0
            _downfactors[downset] .= __downfactors[downset]

            # params
            sim_status = "RUNNING"
            last_biom = 0.0

            # origin model
            lb!(opm, lb0)
            ub!(opm, ub0)

            # for each step
            for _downregi in 1:Int(1e10)
                
                # peek next id
                ridx = nothing
                for _ridx_att in 1:500
                    ridx = rand(iidxs_pool0)
                    iszero(_downfactors[ridx]) && continue
                    _downfactors[ridx] *= down_factor
                    push!(downset, ridx)
                    break
                end
                if isnothing(ridx) 
                    sim_status = "ERROR.FAILED_STEP"
                    @goto END_TRAJ
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
                        sim_status = "ERROR.ONOPTIMIZATION"
                        push!(biomset, NaN)
                        @show(err)
                        @goto END_TRAJ
                end
                # @show ridx
                # @show last_biom

                # feasibility check
                if last_biom < feasible_th 
                    sim_status = "UNFEASIBLE"
                    @goto END_TRAJ
                end
                
                # end trajectory
                @label END_TRAJ

                # check dupplicate
                koset_hash = unordered_hash(downset)
                if check_duplicate!(dups_buff, koset_hash)
                    dups_count += 1
                    @goto RANDOM_RESTART_POINT
                end
                nondups_count += 1

                # break
                break

            end # for _downregi 
            
            # store
            b = blob!(hd_bb)
            merge!(b, @litescope())
            b["cargo.downset", "downset"] = downset
            b["cargo.biomset", "biomset"] = biomset
            
            # blob control
            bc = blobcount(hd_bb)
            # info
            if iszero(mod(bc, 100))
                @info("hit.and.down",
                    bc,
                    dups_count,
                    nondups_count,
                    dumps_frac = dups_count / nondups_count
                )
            end
            
            isfullbatch(hd_bb) && @goto CLOSE_BATCH

        end # for _seti

        # finish
        @label CLOSE_BATCH
        setmeta!(hd_bb, "hit.and.down.version", script_version)
        serialize(hd_bb)
        unlock(hd_bb)

    end # for _batchi

    # write globals
    hnd_globs["iidxs_pool0"] = iidxs_pool0
    serialize(hnd_globs)
    serialize(sim_globs)

    nothing
end