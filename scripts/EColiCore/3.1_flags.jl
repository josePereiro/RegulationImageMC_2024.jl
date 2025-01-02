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
    for bb in eachbatch(B, "hit.and.down"; sortfun)
        islocked(bb) && continue
        lock(bb) do
            @show bb.id
            for b in bb
                # dupplicate
                flag = get!(b, "flags", "duplicate.flag", false)
                downset = b["cargo.koset", "koset"]::Vector{Int}
                koset_hash = combhash(downset)
                flag = flag || check_duplicate!(dups_buff, koset_hash)
                b["flags","duplicate.flag"] = flag

                # hnd.finished
                flag = get!(b, "flags", "done.flag", false)
                status = b["exit_status"]
                flag = flag || status != "RUNNING"
                b["flags", "done.flag"] = flag
            end
            serialize!(bb, "flags")

            # info
            @info("DUPS", 
                dups_count = dups_count(dups_buff),
                nondups_count = nondups_count(dups_buff),
                dups_frac = dup_ratio(dups_buff)
            )

        end
    end
end