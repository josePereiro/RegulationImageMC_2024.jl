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