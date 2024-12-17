@time begin
    using RegulationImageMC_2024
    using ProjFlows
    using Base.Threads
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.99_sim.base.jl")
include("2.99_get.net0.base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # reset
    empty!(G)

    # script meta
    script_id = "gen.net0"
    script_ver = v"0.1.0"
    
    # load net
    netid = "ecoli_core"
    net0 = pull_net(netid)
    @show netid
    @show size(net0)
    
    # set exchanges
    exch_ids = ["EX_glc__D_e"]
    for id in exch_ids
        bounds!(net0, id, -10.0, 0.0)
    end
    biom_id = extras(net0, "BIOM")
    linear_weights!(net0, biom_id, 1.0) 
    
    # create netid globals
    box_eps = 0.0
    box_reduce = true
    _net0_models!(net0; script_id, box_eps, box_reduce)
    
    # finish/write
    merge!(G, script_id, @litecontext)
    G[script_id, "src"] = read(@__FILE__, String)
    G[script_id, "exch_ids"] = exch_ids
    serialize!(G; lk = true)

    nothing
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
