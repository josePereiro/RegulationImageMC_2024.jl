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
    for bb in eachbatch(B, "sample.feasets"; sortfun)
        islocked(bb) && continue
        lock(bb) do
            @show bb.id
            for b in bb
                flag = get!(b, "flags", "duplicate.flag", false)
                feaset = b["cargo.feaset", "feaset"]::Vector{Int}
                koset_hash = combhash(feaset)
                flag = flag || check_duplicate!(dups_buff, koset_hash) 
                b["flags","duplicate.flag"] = flag
                flag ? (dups_count += 1) : (nondups_count += 1)
            end
            serialize!(bb, "flags")

            # info
            @info("DUPS", 
                dups_count,
                nondups_count,
                dips_frac = dups_count / (dups_count + nondups_count)
            )
        end
    end
end