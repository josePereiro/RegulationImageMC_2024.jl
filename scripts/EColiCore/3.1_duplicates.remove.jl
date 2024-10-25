@time begin
    using Base.Threads
    using CairoMakie
    using NDHistograms
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.0_sim.base.jl")

# --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # globals blobs
    sim_globs = blob(B, "sim.globals")
    hnd_globs_id = sim_globs["hnd.globals.id"]
    hnd_globs = blob!(B, hnd_globs_id)

    # duplicates
    dups_buff = HashTracker(5050)

    sortfun = shuffle
    for bb in eachbatch(B, hnd_globs_id; sortfun)
        lock(bb) do
            @show bb.group
            for b in bb
                flag = get(b, "duplicate.flag", false)
                downset = b["cargo.downset", "downset"]::Vector{Int}
                koset_hash = unordered_hash(downset)
                b["duplicate.flag"] = flag || check_duplicate!(dups_buff, koset_hash)
            end
            serialize(bb)
        end
    end
end