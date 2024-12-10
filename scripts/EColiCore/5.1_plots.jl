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
# feasets hist
let
    # clear
    h0_ref = blobyio!(G, :get!, 
        "plots.feasets.hist",
        "cached"
    ) do
        _h0 = NDHistogram("feaset.len" => 0:1000)
        bb_count = 0
        @time for bb in eachbatch(B, "gen.powersets")
            @show bb.id
            for b in bb
                # get(b, "duplicate.flag", false) && continue
                feaset = b["cargo.feaset", "feaset"]
                feaset_len = length(feaset)
                # @show feaset
                count!(_h0, (feaset_len, ), 1)
            end
            bb_count += 1
            bb_count < 100 || break
        end
        return _h0
    end
    lock(G) do
        serialize!(G)
    end
    nothing
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    
    h0 = G["plots.feasets.hist", "cached"]

    # Plots
    @time xs, ws = hist_series(h0, "feaset.len")

    f = Figure()
    ax = Axis(f[1,1]; 
        title = G["gen.net0", "netid"],
        xlabel = "feaset.len", 
        ylabel = "count"
    )
    scatter!(ax, xs, ws)
    barplot!(ax, xs, ws)
    f
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# fba.v hist
let
    return
    # params
    net0 = G["gen.net0", "net0.ref"][]

    h0_ref = blobyio!(G, :set!, 
        "plots.fba.sol.hist", 
        "cached"
    ) do
        _h0 = NDHistogram(
            "EX_O2"  => -1000.0:0.01:1000.0,
            "EX_CO2" => -1000.0:0.01:1000.0,
            "EX_GLC" => -1000.0:0.01:1000.0,
            "EX_NH4" => -1000.0:0.01:1000.0,
            "EX_GLU" => -1000.0:0.01:1000.0,
            "BIOM"   => -1000.0:0.01:1000.0,
            "ATPM"   => -1000.0:0.01:1000.0,
        )
        bb_count = 0
        @time for bb in eachbatch(B, "fba.feasures")
            @show bb.id
            for b in bb
                # get(b, "duplicate.flag", false) && continue
                sol = b["cargo.fba", "sol"]
                isempty(sol) && continue

                function _sol(cid) 
                    return sol[colindex(net0, extras(net0, cid))]
                end

                count!(_h0, (
                        _sol("EX_O2"),
                        _sol("EX_CO2"),
                        _sol("EX_GLC"),
                        _sol("EX_NH4"),
                        _sol("EX_GLU"),
                        _sol("BIOM"),
                        _sol("ATPM" ),
                    ), 
                    1
                )
            end
            bb_count += 1
            bb_count < 1500 || break
        end
        return _h0
    end
    lock(G) do
        serialize!(G)
    end
   
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# fba 1D hist
let
    h0 = G["plots.fba.sol.hist", "cached"]

    # Plots
    # "EX_O2", "EX_CO2", "EX_GLC", "EX_NH4", "EX_GLU", "BIOM", "ATPM"
    id = "BIOM"
    h1 = marginal(h0, id)
    h1 = rebin(h1, id => -50:0.02:50)
    @time xs, ws = hist_series(h1, id)

    f = Figure()
    ax = Axis(f[1,1]; 
        title = G["gen.net0", "netid"],
        xlabel = id, 
        ylabel = "count"
    )
    scatter!(ax, xs, ws)
    barplot!(ax, xs, ws)
    f
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# corr
let
    h0 = G["plots.fba.sol.hist", "cached"]

    # Plots
    # "EX_O2", "EX_CO2", "EX_GLC", "EX_NH4", "EX_GLU", "BIOM", "ATPM"
    id1 = "EX_GLC"
    # h1 = marginal(h0, id1)
    h1 = rebin(h0, id1 => -500:0.01:500)
    x1s = collect(keys(h1, id1))
    # @time x1s, w1s = hist_series(h1, id1)

    id2 = "BIOM"
    # h2 = marginal(h0, id2)
    h2 = rebin(h0, id2 => -500:0.01:500)
    x2s = collect(keys(h2, id2))
    # @time x2s, w2s = hist_series(h2, id2)

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
