# @time begin
#     using RegulationImageMC_2024
#     using Gurobi
#     using Clp
#     using ProjFlows
#     using Base.Threads
# end

# # --.-...- --. -. - -.-..- -- .-..- -. -. 
# include("0.0_proj.jl")
# include("1.99_sim.base.jl")
# include("2.99_get.net0.base.jl")

# ## --.-...- --. -. - -.-..- -- .-..- -. -. 
# # DOING: create config global
# ## - This setup the current arguments for the scripts
# ## - Ex: net0, objid, etc
# ## - That is, the simulation is stateful 

# ## --.-...- --. -. - -.-..- -- .-..- -. -. 
# let
#     # net0 globals
#     netid = "iJR904"
#     net0_globs_id = string(netid, ".globals")
#     net0_globs = blob!(B, net0_globs_id)
#     net0_globs["net0.netid"] = netid

#     # load net
#     net0 = pull_net(netid)
#     @show netid
#     @show size(net0)

#     # set glc exchange
#     # It should be the only carbon source open
#     exch_ids = [extras(net0, "EX_GLC")]
#     for id in exch_ids
#         bounds!(net0, id, -10.0, 0.0)
#     end
#     net0_globs["exch_ids"] = exch_ids
#     biom_id = extras(net0, "BIOM")
#     net0_globs["net0.biom_id"] = biom_id
#     linear_weights!(net0, biom_id, 1.0) 

#     # create netid globals
#     _net0_globals!(net0_globs, net0; 
#         box_eps = 1e-2, 
#         box_reduce = false, 
#         box_nths = NTHREADS
#     )
#     serialize(net0_globs)
    
#     # up sim globals
#     merge!(G, @litescope)
#     serialize(G)

#     nothing
# end

# ## --.-...- --. -. - -.-..- -- .-..- -. -. 