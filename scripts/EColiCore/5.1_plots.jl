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
# feasets hist
let
    # clear
    h0_ref = blobyio!(C, 
        "plots.feasets.hist",
        "cached", 
        :setser!, 
    ) do

        n_tasks = 2 * nthreads()
        ch_size = 10

        bbch = eachbatch(B, "gen.powersets"; 
            n_tasks, ch_size
        )
        
        tasks = map(1:n_tasks) do _
            @spawn let
                _h0 = NDHistogram(
                    "feaset.len" => 0:1000,
                    "ko.indx" => 0:1000,
                )
                bb_count = 0
                @time for bb in bbch
                    @show (bb.id, threadid())

                    for b in bb
                        feaset = b["cargo.feaset", "feaset"]
                        feaset_len = length(feaset)
                        # @show feaset
                        for idx in feaset
                            count!(_h0, (feaset_len, idx), 1)
                        end
                    end # for b

                    bb_count += 1
                    bb_count < Inf || break
                end # for bb
                return _h0
            end # @spawn let
        end # map task
        
        return merge!(map(fetch, tasks)...)
    end # blobyio!
    nothing
end


## --.-...- --. -. - -.-..- -- .-..- -. -. 
# PLOT
let
    h0 = C["plots.feasets.hist", "cached"]
    blep0 = G["gen.net0", "net0.blep0.ref"][]
    
    # Plots
    # "ko.indx", "feaset.len"
    # did = "feaset.len"
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
# fba.v hist
let
    # params
    blep0 = G["gen.net0", "net0.blep0.ref"][]
    net0 = G["gen.net0", "net0.ref"][]

    _INDEX_MAP = Dict(
        eid => colindex(blep0, extras(net0, eid))
        for eid in [
            "EX_O2", "EX_CO2", "EX_GLC", "EX_NH4", "EX_GLU", "BIOM", "ATPM"
            ]
    )

    n_tasks = 2 * nthreads()
    ch_size = 10

    @time h0_ref = blobyio!(C, 
        "plots.fba.sol.hist", 
        "cached", 
        :getser!
    ) do
        bbch = eachbatch(B, "fba.feasures"; 
            ch_size, n_tasks
        )

        bb_count = 0
        tasks = map(1:n_tasks) do _
            @spawn let
                _h0 = NDHistogram(
                    "EX_O2"  => -1000.0:0.01:1000.0,
                    "EX_CO2" => -1000.0:0.01:1000.0,
                    "EX_GLC" => -1000.0:0.01:1000.0,
                    "EX_NH4" => -1000.0:0.01:1000.0,
                    "EX_GLU" => -1000.0:0.01:1000.0,
                    "BIOM"   => -1000.0:0.001:1000.0,
                    "ATPM"   => -1000.0:0.01:1000.0,
                    "EX_O2/InCmol"  => -1000.0:0.001:1000.0,
                    "EX_CO2/InCmol" => -1000.0:0.001:1000.0,
                    "EX_GLC/InCmol" => -1000.0:0.001:1000.0,
                    "EX_NH4/InCmol" => -1000.0:0.001:1000.0,
                    "EX_GLU/InCmol" => -1000.0:0.001:1000.0,
                    "BIOM/InCmol"   => -1000.0:0.0001:1000.0,
                    "ATPM/InCmol"   => -1000.0:0.001:1000.0,
                )
                
                for bb in bbch
                    @show (bb.id, threadid())
                    for b in bb
                        # get(b, "duplicate.flag", false) && continue
                        sol = b["cargo.fba", "sol"]
                        isempty(sol) && continue

                        InCmol = abs(6 * sol[_INDEX_MAP["EX_GLC"]])

                        count!(_h0, (
                                sol[_INDEX_MAP["EX_O2"]],
                                sol[_INDEX_MAP["EX_CO2"]],
                                sol[_INDEX_MAP["EX_GLC"]],
                                sol[_INDEX_MAP["EX_NH4"]],
                                sol[_INDEX_MAP["EX_GLU"]],
                                sol[_INDEX_MAP["BIOM"]],
                                sol[_INDEX_MAP["ATPM"]],
                                sol[_INDEX_MAP["EX_O2"]]/InCmol,
                                sol[_INDEX_MAP["EX_CO2"]]/InCmol,
                                sol[_INDEX_MAP["EX_GLC"]]/InCmol,
                                sol[_INDEX_MAP["EX_NH4"]]/InCmol,
                                sol[_INDEX_MAP["EX_GLU"]]/InCmol,
                                sol[_INDEX_MAP["BIOM"]]/InCmol,
                                sol[_INDEX_MAP["ATPM"]]/InCmol,
                            ), 
                            1
                        )
                    end # for b 
                    bb_count += 1
                    bb_count < Inf || break
                end # for bb
                return _h0
            end # @spawn let
        end # tasks = map

        return merge!(map(fetch, tasks)...)
        
    end # blobyio!
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# fba 1D hist
let
    h0 = C["plots.fba.sol.hist", "cached"]
    m0, m1 = extrema(keys(h0, "BIOM"))
    h1 = filter(h0) do v, w
        biom = v[dimindex(h0, "BIOM")]
        biom > m1 * 0.0 || return false
        biom < m1 * 0.8 || return false
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
        ylabel = "log10 count"
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
    h0 = C["plots.fba.sol.hist", "cached"]
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
    id1 = "EX_O2"
    # h1 = marginal(h0, id1)
    h2 = rebin(h1, id1 => -500:0.001:500)
    x1s = collect(keys(h2, id1))
    # @time x1s, w1s = hist_series(h1, id1)

    id2 = "EX_CO2"
    h2 = rebin(h1, id2 => -500:0.001:500)
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
