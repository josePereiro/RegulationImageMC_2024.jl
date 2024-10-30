@time begin
    using RegulationImageMC_2024
    using Gurobi
    using ProjFlows
    using Random
    using DataStructures
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.99_sim.base.jl")

# --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # globals blobs
    hnd_globs_id = G["hnd.globals.id"]
    hnd_globs = blob(B, hnd_globs_id)
    
    psets_globs_id = string()
    psets_globs = blob!(B, psets_globs_id)
    G["psets.globals.id"] = psets_globs_id

    @time for bb in eachbatch(B, hnd_globs_id)
        for b in bb
            get(b, "duplicate.flag", false) && continue
            
            count!(h0, (biomset_len, ), 1)
        end
    end
end


# --.-...- --. -. - -.-..- -- .-..- -. -. 