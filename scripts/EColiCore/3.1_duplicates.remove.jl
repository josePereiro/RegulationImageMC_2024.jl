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
    
    # # TODO: rename to account for sim.global
    # hnd_id = "hit.and.down" # TO SYNC
    # hnd_globs_id = "$(hnd_id).$(netid).globals"
    # hnd_globs = blob!(B, hnd_globs_id)
    sim_globs
end