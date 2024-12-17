# @time begin
#     using RegulationImageMC_2024
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
#     netid = "ecoli_core_Beg2007"
#     net0_globs_id = string(netid, ".globals")
#     net0_globs = blob!(B, net0_globs_id)
#     net0_globs["net0.netid"] = netid

#     # load net
#     net0 = pull_net(netid)
#     @show netid
#     @show size(net0)

#     # set exchange pattern
#     # Set exchanges, all nutrinets allowed
#     exch_ids = [ "EX_glc__D_e", "EX_lac__D_e", "EX_malt_e", "EX_gal_e", "EX_glyc_e"]
#     for id in exch_ids
#         bounds!(net0, id, -10.0, 0.0)
#     end
#     # TODO: try with no acetate
#     bounds!(net0, "EX_ac_e", 0.0, 1000.0) # allow production
#     net0_globs["exch_ids"] = exch_ids
#     biom_id = extras(net0, "BIOM")
#     net0_globs["net0.biom_id"] = biom_id
#     linear_weights!(net0, biom_id, 1.0) 

#     # net types
#     _net0_models!(net0_globs, net0)
#     serialize(net0_globs)
    
    
#     # up sim globals
#     merge!(G, @litescope)
#     serialize(G)

#     nothing
# end

# ## --.-...- --. -. - -.-..- -- .-..- -. -. 