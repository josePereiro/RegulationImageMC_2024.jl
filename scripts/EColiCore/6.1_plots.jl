# --.-...- --. -. - -.-..- -- .-..- -. -. 
@time begin
    using Base.Threads
    using CairoMakie
    using NDHistograms
    using Base.Threads
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.99_sim.base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# PLOT
let
    blep0 = G["gen.net0", "net0.blep0.ref"][]
    h0 = C["feasets.hist", "h0"]
    
    # Plots
    # "ko.indx", "feaset.len"
    did = "feaset.len"
    # did = "ko.indx"
    h1 = marginal(h0, did)
    @time x0s, ws = hist_series(h1, did)
    si = sortperm(ws)
    xs = eachindex(x0s)
    ws = ws[si]

    f = Figure()
    ax = Axis(f[1,1]; 
        title = G["gen.net0", "netid"],
        xlabel = did, 
        ylabel = "count"
    )
    scatter!(ax, xs, ws)
    barplot!(ax, xs, ws)
    # ax.xticks = (collect(xs), colids(blep0, x0s))
    # ax.xticklabelrotation = 45
    f
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# fba 1D hist
let
    h0 = C["fba.sol.hist", "h0"]
    fid = "BIOM/InCmol"
    m0, m1 = extrema(keys(h0, fid))
    h1 = filter(h0) do v, w
        biom = v[dimindex(h0, fid)]
        biom > m1 * 0.0 || return false
        biom < m1 * 0.9 || return false
        return true
    end

    # Plots
    # "EX_O2", "EX_CO2", "EX_GLC", "EX_NH4", "EX_GLU", "BIOM", "ATPM"
    id = "BIOM/InCmol"
    h1 = marginal(h1, id)
    h1 = rebin(h1, id => -50:0.001:50)
    @time xs, ws = hist_series(h1, id)

    f = Figure()
    ax = Axis(f[1,1]; 
        title = G["gen.net0", "netid"],
        xlabel = id, 
        ylabel = "count"
    )
    scatter!(ax, xs, ws)
    barplot!(ax, xs, ws)
    # scatter!(ax, xs, log10.(ws))
    # barplot!(ax, xs, log10.(ws))
    # ax.yticks = ([1:10;], string.(1:10))
    # ax.xticks = ([0:0.001:10;], string.(0:0.001:10))
    # ax.xticklabelrotation = 45
    f
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# corr
let
    h0 = C["fba.sol.hist", "h0"]
    m0, m1 = extrema(keys(h0, "BIOM"))
    @show m0, m1
    h1 = filter(h0) do v, w
        biom = v[dimindex(h0, "BIOM")]
        biom > m1 * 0.0 || return false
        # biom > m1 * 0.8 && return false
        return true
    end
    @show length(h1)

    # Plots
    # "EX_O2", "EX_CO2", "EX_GLC", "EX_NH4", "EX_GLU", "BIOM", "ATPM"
    # id1 = "BIOM/InCmol"
    id1 = "EX_NH4"
    # h1 = marginal(h0, id1)
    # h2 = rebin(h1, id1 => -500:0.001:500)
    h2=h1
    x1s = collect(keys(h2, id1))
    @show length(x1s)
    # @time x1s, w1s = hist_series(h1, id1)

    id2 = "BIOM/InCmol"
    # h2 = rebin(h1, id2 => -500:0.001:500)
    h2=h1
    x2s = collect(keys(h2, id2))
    # @time x2s, w2s = hist_series(h2, id2)
    @show length(x2s)

    f = Figure()
    ax = Axis(f[1,1]; 
        title = G["gen.net0", "netid"],
        xlabel = id1, 
        ylabel = id2, 
    )
    scatter!(ax, x1s, x2s)
    # barplot!(ax, xs, ws)
    f
end
