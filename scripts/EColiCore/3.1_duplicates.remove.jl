@time begin
    using RegulationImageMC_2024
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.99_sim.base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # globals blobs
    hnd_globs = getframe(G, "hnd")
    hnd_globs_id = hnd_globs["hnd_fullid"]

    # duplicates
    dups_buff = HashTracker(5050)

    sortfun = shuffle
    for bb in eachbatch(B, hnd_globs_id; sortfun)
        lock(bb) do
            @show bb.group
            for b in bb
                flag = get(b, "duplicate.flag", false)
                downset = b["cargo.downset", "downset"]::Vector{Int}
                koset_hash = combhash(downset)
                b["duplicate.flag"] = flag || check_duplicate!(dups_buff, koset_hash)
            end
            serialize(bb)
        end
    end
end