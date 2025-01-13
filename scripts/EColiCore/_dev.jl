@time begin
    using RegulationImageMC_2024
    using ProjFlows
    using Random
    using DataStructures
    using Combinatorics
    using StatsBase
    using CairoMakie
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.99_sim.base.jl")

# --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    net0 = pull_net("ecoli_core")
    solver = QUAD_LP_SOLVER
    opm = FBAOpModel(net0, solver)
    
    # biom
    delta = 1e-3
    id = extras(net0, "BIOM")
    @show id
    set_linear_obj!(opm, id, MAX_SENSE)
    optimize!(opm)
    obj = solution(opm, id)
    @show obj
    l, u = bounds(opm, id)
    lb!(opm, id, max(obj - (obj * delta), l))
    ub!(opm, id, min(obj + (obj * delta), u))

    # glc
    delta = 1e-3
    id = extras(net0, "EX_GLC")
    @show id
    set_linear_obj!(opm, id, MIN_SENSE)
    try
        optimize!(opm)
    catch err
        # @show err
    end 
    obj = solution(opm, id)
    @show obj
    l, u = bounds(opm, id)
    lb!(opm, id, max(obj - (obj * delta), l))
    ub!(opm, id, min(obj + (obj * delta), u))

    f = Figure()
    ax = Axis(f[1,1])

    # parsimonio
    ids = colids(opm)
    # Max
    set_v2_obj!(opm, MIN_SENSE)
    optimize!(opm)
    v2_min_sol = solution(opm, ids)

    # Max
    set_v2_obj!(opm, MAX_SENSE)
    optimize!(opm)
    v2_max_sol = solution(opm, ids)
    
    
    # diff = clamp.(v2_max_sol .- v2_min_sol, -10.0, 10.0)
    # scatter!(ax, eachindex(sol), sol)
    scatter!(ax, v2_max_sol, v2_min_sol;
        scale = :
    )
    
    f
    
    # set_v2_obj!(opm, MIN_SENSE)
    # 
    # vec = solution(opm, ids)
    
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 


# --.-...- --. -. - -.-..- -- .-..- -. -. 

## -. -.- -. -. - ..-.... - - . . .- .- .- -. -...
# DONE: add blob inspection
let

    n_tasks = 40
    ch_size = 40

    bbch = eachbatch(B, "fba.feasures"; 
        n_tasks = 1, ch_size
    )

    bb_count = 100
    @time tasks = map(1:n_tasks) do _
        @spawn let
            
            for bb in bbch
                @show (bb.id, threadid())
                for b in bb
                    
                    sol = b["cargo.fba", "sol"]::Vector{Float64}
                    isempty(sol) && continue

                end # for b 
                bb_count += 1
                bb_count < Inf || break
            end # for bb
            return rand()
        end # @spawn let
    end # tasks = map

    return map(fetch, tasks)
end

## -. -.- -. -. - ..-.... - - . . .- .- .- -. -...
# DOING Sampling downset
let
    downset0 = [1]
    feaset = downset0[1:end-1]
    isempty(feaset) && continue
    D = length(feaset)
    subfeaset = StatsBase.sample(feaset, rand(1:D), replace = false)
end

## -. -.- -. -. - ..-.... - - . . .- .- .- -. -...
## -. -.- -. -. - ..-.... - - . . .- .- .- -. -...
## -. -.- -. -. - ..-.... - - . . .- .- .- -. -...
## -. -.- -. -. - ..-.... - - . . .- .- .- -. -...
using HTTP
using JSON

function retrieve_data(url::String)
    response = HTTP.get(url)
    if response.status == 200
        try
            data = JSON.parse(String(response.body))
            return data
        catch err
            return response.body
        end
    else
        error("Failed to retrieve data: HTTP status code $(response.status)")
    end
end

## -. -.- -. -. - ..-.... - - . . .- .- .- -. -...
let
    rxn = "ICDHyr"
    rxn_url = "http://bigg.ucsd.edu/api/v2/universal/reactions"
    res = retrieve_data(
        joinpath(rxn_url, rxn)
    )
    dblinks = get(res, "database_links", nothing)
    dblinks = get(dblinks, "KEGG Reaction", nothing)
    for linkdat in dblinks
        @show linkdat["id"]
        @show linkdat["link"]
        res = retrieve_data(
            # joinpath(linkdat["link"], "json")
            linkdat["link"]
        )
        return res
    end
end

## -. -.- -. -. - ..-.... - - . . .- .- .- -. -...
# https://rest.kegg.jp/get/<dbentries>[/<option>]
# R00267
# retrieve_data("https://rest.kegg.jp/get/r:R00267/json")

let
    # Bigg
    Bigg_rxn = "GLUN"
    rxn_url = "http://bigg.ucsd.edu/api/v2/universal/reactions"
    Bigg_res = retrieve_data(joinpath(rxn_url, Bigg_rxn))
    dblinks = get(Bigg_res, "database_links", nothing)
    dblinks = get(dblinks, "KEGG Reaction", nothing)
    Kegg_rxn = dblinks[1]["id"]

    # Kegg
    kegg_reaction = KEGGAPI.link("enzyme", Kegg_rxn)
    kegg_reaction_info = KEGGAPI.kegg_get([kegg_reaction.data[2][1]])
    split(kegg_reaction_info[2][1], "\n")
end

## -. -.- -. -. - ..-.... - - . . .- .- .- -. -...
## -. -.- -. -. - ..-.... - - . . .- .- .- -. -...
## -. -.- -. -. - ..-.... - - . . .- .- .- -. -...
let
    script_id = "_dev"
    script_ver = v"0.1.0"

    # local cache
    S = blobbatch!(B, 
        hashed_id(
            string("cache.", script_id), 
            script_ver
        )
    )
    S["script_id"] = script_id
    S["script_ver"] = script_ver

    for it in 1:10

        # duplicate buffer
        dups_buff = _dups_tracker!(S; 
            dup_buff_size = 1_000_000
        )
        @info("HASH_SET", 
            dups_buff = length(dups_buff),
        )

        for s in 1:100
            check_duplicate!(dups_buff, rand(UInt64))
        end
                    
    end
end

## -. -.- -. -. - ..-.... - - . . .- .- .- -. -...
let
   1 
end


# -. -.- -. -. - ..-.... - - . . .- .- .- -. -...