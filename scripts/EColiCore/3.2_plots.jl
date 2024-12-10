@time begin
    using Base.Threads
    using CairoMakie
    using NDHistograms
    const At =  Base.Threads.Atomic
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.99_sim.base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# Biomass trajectory
let
    # globals blobs
    
    f = Figure()
    ax = Axis(f[1,1]; xlabel = "idx", ylabel = "biom")

    for bb in eachbatch(B, hnd_globs_id)
        @show bb.group
        for b in bb
            exit_status = b["exit_status"]
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
# koset len
let
    # globals blobs
    hnd_globs = getframe(G, "hnd")
    hnd_globs_id = hnd_globs["hnd_fullid"]
    
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
# biom dist
let
    # globals blobs
    hnd_globs = getframe(G, "hnd")
    hnd_globs_id = hnd_globs["hnd_fullid"]
    
    h0 = NDHistogram("biom" => 0:0.05:1000)

    @time for bb in eachbatch(B, hnd_globs_id)
        for b in bb
            get(b, "duplicate.flag", false) && continue
            biomset = b["cargo.biomset", "biomset"]
            for (bi, biom) in enumerate(biomset)
                isfinite(biom) || continue
                count!(h0, (biom, ), bi)
            end
        end
    end

    # Plots
    @time xs, ws = hist_series(h0, "biom")

    f = Figure()
    ax = Axis(f[1,1]; xlabel = "biom", ylabel = "count")
    scatter!(ax, xs, ws)
    barplot!(ax, xs, ws)
    f
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # globals blobs
    hnd_globs = getframe(G, "hnd")
    hnd_globs_id = hnd_globs["hnd_fullid"]

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
