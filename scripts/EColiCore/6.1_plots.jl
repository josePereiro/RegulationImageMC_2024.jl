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
#MARK: DEV
let
    _bb = B[r"\Afba.feasures"][]
    # _sol = ["cargo.fba.max.v2", "sol"] 
    f = Figure()
    ax = Axis(f[1,1])
    
    for _b in eachblob(_bb)
        
        frame = "cargo.fba.max.biom"
        sol = _b[frame, "sol"]
        scatter!(ax, eachindex(sol), sol; label = frame)
        
        frame = "cargo.fba.min.glc"
        sol = _b[frame, "sol"]
        scatter!(ax, eachindex(sol), sol; label = frame)
        
        frame = "cargo.fba.max.v2"
        sol = _b[frame, "sol"]
        scatter!(ax, eachindex(sol), sol; label = frame)
        
        
        frame = "cargo.fba.min.v2"
        sol = _b[frame, "sol"]
        scatter!(ax, eachindex(sol), sol; label = frame)
        break
    end
    f
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
#MARK: koset.hist.len
let
    blep0 = G["gen.net0", "net0.blep0.ref"][]
    h0 = C["koset.hist", "h0"]
    
    # Plots
    did = "koset.len"
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
    scatter!(ax, xs, ws ./ xs)
    barplot!(ax, xs, ws ./ xs)
    # ax.xticks = (collect(xs), colids(blep0, x0s))
    # ax.xticklabelrotation = 45
    f
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
#MARK: koset.hist.ko.indx
let
    blep0 = G["gen.net0", "net0.blep0.ref"][]
    h0 = C["koset.hist", "h0"]
    
    # Plots
    did = "ko.indx"
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
    ax.xticks = (collect(xs), colids(blep0, x0s))
    ax.xticklabelrotation = 45
    f
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
#MARK: feasets.hist
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

    # normalize
    ws = ws ./ xs
    scatter!(ax, xs, ws)
    barplot!(ax, xs, ws)
    ax.xticks = (collect(xs), colids(blep0, x0s))
    ax.xticklabelrotation = 45
    f
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
#MARK: feasets.fba.1D.hist
let
    h0 = C["fba.sol.hist", "h0"]

    fid = "BIOM/InCmol"
    m0, m1 = extrema(keys(h0, fid))
    h1 = filter(h0) do v, w
        biom = v[dimindex(h0, fid)]
        biom > m1 * 0.0 || return false
        # biom < m1 * 0.98 || return false
        return true
    end

    # Plots
    # "EX_O2", "EX_CO2", "EX_GLC", "EX_NH4", "EX_GLU", "BIOM", "ATPM"
    # id = "EX_GLC"
    id = "BIOM/InCmol" 
    h1 = marginal(h1, id)
    h1 = rebin(h1, id => -50:0.001:50)
    @time xs, ws = hist_series(h1, id)
    xs = xs * 6^(-1) * 1e3 * 180.156^(-1)

    f = Figure()
    ax = Axis(f[2,1]; 
        title = G["gen.net0", "netid"],
        xlabel = id, 
        ylabel = "count",
        # limits = (0.0, 0.03, 0.0, 5e5)
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
#MARK: feasets.fba.corrs
let
    h0 = C["fba.sol.hist", "h0"]
    m0, m1 = extrema(abs, keys(h0, "EX_GLC"))
    @show m0, m1
    h1 = filter(h0) do v, w
        val = v[dimindex(h0, "EX_GLC")]
        abs(val) > m1 * 0.98 || return false
        abs(val) < m1 * 1.02 || return false
        val = v[dimindex(h0, "BIOM")]
        abs(val) > 1e-4 || return false
        return true
    end
    @show length(h1)

    # Plots
    # "EX_O2", "EX_CO2", "EX_GLC", "EX_NH4", "EX_GLU", "BIOM", "ATPM"
    # id1 = "EX_O2"
    id1 ="EX_CO2"
    # id2 = "BIOM/InCmol"
    id2 = "feaset.len"
    h2 = rebin(h1, 
        id1 => -500:2.5:500,
        # id1 => support(h1, id1),
        id2 => -500:0.001:500,
        # id2 => support(h1, id2),
    )
    # h2 = h1

    ws = collect(values(h2))
    sidx = sortperm(ws)
    ws = ws[sidx]

    x1s = collect(keys(h2, id1))
    x1s = x1s[sidx]
    @show length(x1s)
    
    x2s = collect(keys(h2, id2))
    x2s = x2s[sidx]
    @show length(x2s)

    f = Figure()
    ax = Axis(f[1,1]; 
        title = G["gen.net0", "netid"],
        xlabel = id1, 
        ylabel = id2, 
        # limits = (0.0, 10.3, nothing, nothing)
    )
    p = scatter!(ax, x1s, x2s;
        # color = log.(ws ./ maximum(ws)),
        color = log.(ws),
        colormap = :viridis,  # Mappa colori adatta a istogrammi
        # marker = "rect",
        markersize = 30
    )
    
    Colorbar(f[1, 2], p, label = "~ log(count)")
    f
end

