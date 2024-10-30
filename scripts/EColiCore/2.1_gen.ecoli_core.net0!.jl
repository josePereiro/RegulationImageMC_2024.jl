@time begin
    using RegulationImageMC_2024
    using Gurobi
    using ProjFlows
    using Base.Threads
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.99_sim.base.jl")
include("2.99_get.net0.base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # net0 globals
    netid = "ecoli_core"

    # load net
    net0 = pull_net(netid)
    @show netid
    @show size(net0)
    
    # set glc exchange
    # It should be the only carbon source open
    exch_ids = ["EX_glc__D_e"]
    for id in exch_ids
        bounds!(net0, id, -10.0, 0.0)
    end
    G["net0", "exch_ids"] = exch_ids
    biom_id = extras(net0, "BIOM")
    linear_weights!(net0, biom_id, 1.0) 

    # create netid globals
    _net0_globals!(G, net0)
    merge!(G, @litescope)
    merge!(G, "net0", @litescope)
    
    serialize(G)

    nothing
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 