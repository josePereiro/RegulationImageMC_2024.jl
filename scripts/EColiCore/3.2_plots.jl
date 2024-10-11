@time begin
    using Base.Threads
    using CairoMakie
    using NDHistograms
    using Random
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
# DOING: 
# - make hit and run version
# - reduce risk of duplicates


# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.0_sim.base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# down set lenght
@time let

    netid = "ecoli_core"
    biom0 = B["$(netid).globals"]["net0.biom0"]
    @show biom0

    bbs = eachbatch(B, "hit.and.down")
    _h0 = NDHistogram(
        "downL" => 0:1000,
        "biom" => range(0.0, 2.2; step = biom0 * 1e-2),
        # "biom" => Float64,
    )
    ts = map(1:nthreads()) do _
        @spawn begin
            t_h = deepcopy(_h0)
            for bb in bbs, b in bb
                downset = b["cargo.downset", "downset"]
                biomset = b["cargo.biomset", "biomset"]::Vector{Float64}
                L = length(downset)
                for idx in 1:L
                    biom = biomset[idx]
                    isfinite(biom) || continue
                    count!(t_h, (idx, biom))
                end
            end
            return t_h
        end # @spawn
    end

    # reduce
    global h0 = deepcopy(_h0)
    for t in ts
        hi = fetch(t)
        merge!(h0, hi)
    end
    
    nothing
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    id = "biom"
    # id = "downL"

    h1 = marginal(h0, id)
    f = Figure()
    ax = Axis(f[1,1]; xlabel = id, ylabel = "count")
    xs = collect(keys(h1, id))
    # ys = map(log10, values(h1))
    ys = collect(values(h1))
    # @show xs[1:10]
    barplot!(ax, xs, ys)
    scatter!(ax, xs, ys)
    f
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# downregulation trajectories
let
    netid = "ecoli_core"
    net0 = B["$(netid).globals"][MetNet, "net0"]
    M, N = size(net0)

    f = Figure()
    ax = Axis(f[1,1]; 
        xlabel = "downreg index",
        ylabel = "rxn index",
        limits = (nothing, nothing, 0, N)
    )
    
    bb_ch = eachbatch(B, "hit.and.down"; sortfun = shuffle!)
    ntrajs = Inf
    c = 0
    @time for bb in bb_ch, b in bb
        # islocked(bb) && continue
        downset = b[Vector{Int}, "cargo.downset", "downset"]
        sort!(downset)
        D = length(downset)
        
        # filters
        D > 3 || continue
        # 73 ∈ downset || continue
        # downset[4] == 63 || continue

        lines!(ax, 1:D, downset; 
            color = :black, 
            alpha = 0.05
        )
        # scatter!(ax, [D], [downset[end]]; 
        #     color = :black, 
        #     alpha = 0.01,
        #     marker = :xcross,
        #     markersize = 18
        # )
        c += 1
        c < ntrajs || break
    end
    f
end