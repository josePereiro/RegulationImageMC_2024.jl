@time begin
    using RegulationImageMC_2024
    using Gurobi
    using ProjFlows
    using Random
    using DataStructures
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.99_sim.base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# NOTES:
# - downset: a downregulation set that lead to an unfeasible network
# - fea_downset: a downregulation set that lead to an  still feasible network

# - Each ensemble has a fixed enviromen

## --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # meta
    script_version = v"0.3.0"
    
    # globals blobs
    netid = G["netid"]::String
    @show netid

    hnd_fullid = _dot_string("hit.and.down", netid, script_version)
    
    # Reset (uncomment to reset)
    foreach_batch(rm, B, hnd_fullid)

    # duplicates
    # TODO make write/read frequent and locked
    # Make it an struct which creates a global (rablob) in a given Bloberia
    dup_buff_size = 1_000_000
    global dups_buff = get!(G, "hnd.downset.duplicates", "local.buffer") do
        HashTracker(dup_buff_size)
    end
    dups_count = 0
    nondups_count = 0
    
    return 
    # params
    down_factor = 0.0
    net0 = net0_globs["net0"]::MetNet
    objid = extras(net0, "BIOM")
    objidx = colindex(net0, objid)
    feasible_th = 1e-2
    traj_lim = 35 # max length of trajectories
    blep0 = net0_globs["net0.blep0.ref"][]
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
    for _batchi in 1:15
        @show _batchi
        
        # new BlobBatch
        hd_bb = headbatch!(B, hnd_fullid)
        if islocked(hd_bb) # ignore locked 
            hd_bb = blobbatch!(B, hnd_fullid)
        end
        setmeta!(hd_bb, "blobs.lim", 500)
        lock(hd_bb)
        
        downset = Int[]
        biomset = Float64[]
        __downfactors = ones(Float64, N)
        _downfactors = ones(Float64, N)

        # refresh hnd_globs
        lock(G) do
            # load disk version
            _dups_buff = get!(blob(G), "hnd.downset.duplicates", "local.buffer") do
                hash_set
            end
            hash_set = dups_buff.hash_set
            _hash_set = _dups_buff.hash_set
            @info("REFRESH HASH_SET", hash_set_len = length(hash_set))
            
            isempty(_hash_set) || push!(hash_set, _hash_set...)
            serialize(hnd_globs, "kosets.duplicates")
        end

        # for each run
        for _runi in 1:Int(1e10)
            
            # random restarting point
            @label RANDOM_RESTART_POINT
            L = length(downset) - 1
            rand_idx0 = L < 1 ? 0 : rand(0:L)
            downset = iszero(rand_idx0) ? Int[] : downset[1:rand_idx0]
            biomset = iszero(rand_idx0) ? Float64[] : biomset[1:rand_idx0]
            __downfactors .= _downfactors
            _downfactors .= 1.0
            _downfactors[downset] .= __downfactors[downset]

            # params
            sim_status = "RUNNING"

            last_biom = isempty(biomset) ? 0.0 : last(biomset)

            # reset model
            lb!(opm, lb0 .* _downfactors)
            ub!(opm, ub0 .* _downfactors)

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
                    # println("isnothing(ridx)")
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
                        # println("catch err;")
                        sim_status = "ERROR.ONOPTIMIZATION"
                        push!(biomset, NaN)
                        @show(err)
                        @goto END_TRAJ
                end
                # @show ridx
                # @show last_biom

                # check nan
                if isnan(last_biom)
                    # println("isnan(last_biom)")
                    sim_status = "OBJ.NAN"
                    @goto END_TRAJ
                end

                # feasibility check
                if last_biom < feasible_th 
                    # println("last_biom < feasible_th ")
                    sim_status = "UNFEASIBLE"
                    @goto END_TRAJ
                end
                
                # end trajectory
                @label END_TRAJ
                
                # continue if running
                sim_status == "RUNNING" && continue
                
                # check dupplicate
                koset_hash = unordered_hash(downset)
                if check_duplicate!(dups_buff, koset_hash)
                    # println("check_duplicate!(dups_buff, koset_hash)")
                    dups_count += 1
                    @goto RANDOM_RESTART_POINT
                end
                nondups_count += 1

                break

            end # for _downregi 
            
            # store
            b = blob!(hd_bb)
            b["sim_status"] = sim_status
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
                    dups_frac = dups_count / nondups_count
                )
            end
            
            isfullbatch(hd_bb) && @goto CLOSE_BATCH

        end # for _runi

        # finish
        @label CLOSE_BATCH
        setmeta!(hd_bb, "hit.and.down.version", script_version)
        setmeta!(hd_bb, "lite.scope", @litescope)
        serialize(hd_bb)
        unlock(hd_bb)

    end # for _batchi

    # write globals
    hnd_globs["iidxs_pool0"] = iidxs_pool0
    serialize(hnd_globs)

    G["lite.scopes", basename(@__FILE__)] = @litescope
    serialize(G)

    nothing
end