# --.-...- --. -. - -.-..- -- .-..- -. -. 
# Creates a DataFrame with FBA features

# --.-...- --. -. - -.-..- -- .-..- -. -. 
@time begin
    using RegulationImageMC_2024
    using Gurobi
    using ProjFlows
    using Random
    using DataFrames
    using Combinatorics
    using CSV
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.0_sim.base.jl")

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

    hnd_id = sim_globs["hit.and.down.id"]
    hnd_globs_id = "$(netid).$(hnd_id).globals"
    hnd_globs = blob!(B, hnd_globs_id)
    iidxs_pool0 = hnd_globs["iidxs_pool0"]::Vector{Int}
    blep0 = net0_globs["net0.blep0"]::LEPModel

    
    objid = extras(blep0, "BIOM")
    glc_id = extras(blep0, "EX_GLC")
    glc_idx = colindex(blep0, glc_id)
    objidx = colindex(blep0, objid)
    lb0, ub0 = lb(blep0), ub(blep0)
    # TODO: use to globals
    exch_ids = ["EX_glc__D_e"]
    exch_idxs = colindex(blep0, exch_ids)
    
    nrxns = size(blep0, 2)
    nirxns = length(iidxs_pool0)
    opm = FBAOpModel(blep0, LP_SOLVER)
    regmagv = ones(nrxns)
    fbasol = ones(nrxns)
    lb!(opm, lb0)
    ub!(opm, ub0)

    # data frame
    coldef = [
        ("downset.hash", UInt); 
        ("down_factor", Float64); 
        [("downfactor.$r", Float64) for r in colids(blep0, iidxs_pool0)];
        [("fba.$r", Float64) for r in colids(blep0, iidxs_pool0)]
        ("fba.biom", Float64)
        [("fba.$r", Float64) for r in exch_ids]
    ]
    df0 = DataFrame([name => T[] for (name, T) in coldef])
    df = deepcopy(df0)

    # duplication handling
    downsets_hashes = Set{UInt64}()
    
    # df batch creation
    nrows = 0
    nrows_per_df = 50_000
    max_nrows = Inf

    # reset
    df_dir = _procdir(["sklearn"])
    rm(df_dir; force = true, recursive = true)
    mkpath(df_dir)

    bb_ch = eachbatch(B, "hit.and.down"; sortfun = shuffle!)
    @time for bb in bb_ch
        # batch filters
        islocked(bb) && continue

        for b in bb
            # TODO: rename "cargo.downset" => "cargo.downset.ko"
            down_factor = b[Float64, "down_factor"]
            downset = b[Vector{Int}, "cargo.downset", "downset"]
            sort!(downset)
            D = length(downset)

            # blob filters
            hash(downset) ∈ downsets_hashes && continue

            # all subsets
            # powerset keep it sorted
            sub_downsets = powerset(downset, 1, D)
            for sub_downset in sub_downsets
                
                subhash = hash(sub_downset)

                # skip duplicates
                subhash ∈ downsets_hashes && continue
                # or add
                push!(downsets_hashes, subhash)

                # FBA
                # prepare opm
                lb!(opm, lb0) 
                ub!(opm, ub0)
                for rxni in sub_downset
                    lb, ub = bounds(opm, rxni)
                    lb!(opm, rxni, lb * down_factor)
                    ub!(opm, rxni, ub * down_factor)
                end
                # fbasol
                try; optimize!(opm) 
                    fbasol .= solution(opm)
                    
                    # if fbasol[glc_idx] < -10.0
                    #     @show fbasol[glc_idx]
                    #     @show fbasol[objidx]
                    #     @show bounds(opm, glc_idx)
                    # end

                    catch err; 
                        continue # skip unfeasible
                end

                # downfactor vec
                _downfactor_vec!(regmagv, sub_downset, down_factor)

                # push df
                push!(df, [
                    subhash; 
                    down_factor;
                    regmagv[iidxs_pool0];
                    fbasol[iidxs_pool0];
                    fbasol[objidx];
                    fbasol[exch_idxs]
                ])

                # info
                iszero(mod(nrows, 100)) && @show nrows

                # write
                dowrite = (nrows >= nrows_per_df) && iszero(mod(nrows, nrows_per_df))
                if dowrite
                    output = joinpath(df_dir, string(
                        "fba.features.df.", nrows, ".csv"
                    ))
                    @show output
                    CSV.write(output, df)
                    df = deepcopy(df0)
                end

                # end
                nrows += 1
                nrows < max_nrows || @goto BB_LOOP_END
                
            end; @label DOWNi_LOOP_END
            
        end; @label B_LOOP_END

    end; @label BB_LOOP_END

    # write last
    output = joinpath(df_dir, string(
        "fba.features.df.", nrows, ".csv"
    ))
    @show output
    CSV.write(output, df)

    # some info
    @show sizeof(df)
    describe(df)

    nothing
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
