@time begin
    using RegulationImageMC_2024
    using Random
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.99_sim.base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # globals blobs
    # duplicates
    DUP_BUFF_SIZE = 20_000_000
    dups_buff = HashTracker(DUP_BUFF_SIZE)

    sortfun = shuffle!
    dups_count = 0
    nondups_count = 0
    for bb in eachbatch(B, "fba.feasures"; sortfun)
        islocked(bb) && continue
        lock(bb) do
            @show bb.id
            ps_bb_ref = bb["fba.feasures", "ps_bb_ref"]
            flag = get!(bb, "fba.feasures", "duplicate.flag", false)
            ref_hash = hash(ps_bb_ref)
            flag = flag || check_duplicate!(dups_buff, ref_hash) 
            bb["fba.feasures", "duplicate.flag"] = flag
            flag ? (dups_count += 1) : (nondups_count += 1)

            # delete!
            flag && rm(bb)

            # info
            @info("DUPS", 
                dups_count,
                nondups_count,
                dips_frac = dups_count / (dups_count + nondups_count)
            )
        end
    end
end