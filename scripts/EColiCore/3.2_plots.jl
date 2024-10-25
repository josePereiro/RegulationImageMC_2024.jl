@time begin
    using Base.Threads
    using CairoMakie
    using NDHistograms
    const At =  Base.Threads.Atomic
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.0_sim.base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# Trajectories
let
    # globals blobs
    sim_globs = blob(B, "sim.globals")
    hnd_globs_id = sim_globs["hnd.globals.id"]
    hnd_globs = blob!(B, hnd_globs_id)
    
    f = Figure()
    ax = Axis(f[1,1]; xlabel = "idx", ylabel = "biom")

    for bb in eachbatch(B, hnd_globs_id)
        @show bb.group
        for b in bb
            sim_status = b["sim_status"]
            # downset = b["cargo.downset", "downset"]::Vector{Int}
            biomset = b["cargo.biomset", "biomset"]::Vector{Float64} 
            idxs = eachindex(biomset)
            lines!(ax, idxs, biomset;
                color = :black, alpha = 0.1
            )
            scatter!(ax, [last(idxs)], [last(biomset)];
                color = :black, alpha = 1.0
            )
        end
    end
    f
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# lenght
let
    # globals blobs
    sim_globs = blob(B, "sim.globals")
    hnd_globs_id = sim_globs["hnd.globals.id"]
    hnd_globs = blob!(B, hnd_globs_id)
    
    h0 = NDHistogram("koset.len" => 0:1000)

    @time for bb in eachbatch(B, hnd_globs_id)
        for b in bb
            get(b, "duplicate.flag", false) && continue
            biomset = b["cargo.biomset", "biomset"]
            biomset_len = length(biomset)
            count!(h0, (biomset_len, ), 1)
        end
    end

    # Plots
    @time xs, ws = hist_series(h0, "koset.len")

    f = Figure()
    ax = Axis(f[1,1]; xlabel = "koset.len", ylabel = "count")
    scatter!(ax, xs, ws)
    f
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # globals blobs
    sim_globs = blob(B, "sim.globals")
    hnd_globs_id = sim_globs["hnd.globals.id"]

    dup_count = At{Int}(0)
    nondup_count = At{Int}(0)
    @time foreach_batch(B, hnd_globs_id) do bb
        for b in bb
            if get(b, "duplicate.flag", false)
                dup_count[] += 1
            else
                nondup_count[] += 1
            end
        end
    end

    @show dup_count[]
    @show nondup_count[]
    nothing
end
