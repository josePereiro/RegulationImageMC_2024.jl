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
#MARK: # koset hist
let
    n_tasks = min(nthreads(), 40)
    ch_size = min(nthreads(), 40)
    
    # clear
    h0_ref = blobio!(C, 
        "koset.hist",
        "h0", 
        :setser!, 
    ) do

        bbch = eachbatch(B, "hit.and.down"; 
            n_tasks = 1, ch_size
        )
        
        tasks = map(1:n_tasks) do _
            @spawn let
                _h0 = NDHistogram(
                    "koset.len" => 0:1000,
                    "ko.indx" => 0:1000,
                )
                bb_count = 0
                for bb in bbch
                    @show (bb.id, threadid())
                    islocked(bb) && continue

                    for b in bb
                        get(b, "flags", "duplicate.flag", true) && continue
                        koset = b["cargo.koset", "koset"]
                        koset_len = length(koset)
                        # @show feaset
                        for idx in koset
                            count!(_h0, (koset_len, idx), 1)
                        end
                    end # for b

                    bb_count += 1
                    bb_count < Inf || break
                end # for bb
                return _h0
            end # @spawn let
        end # map task
        
        @show length(tasks)
        return merge(map(fetch, tasks))

    end # blobio!
    @show h0_ref
    return nothing
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
#MARK: # feasets hist
let
    n_tasks = min(nthreads(), 40)
    ch_size = min(nthreads(), 40)
    
    # clear
    h0_ref = blobio!(C, 
        "feasets.hist",
        "h0", 
        :setser!, 
    ) do

        bbch = eachbatch(B, "sample.feasets"; 
            n_tasks = 1, ch_size
        )
        
        tasks = map(1:n_tasks) do _
            @spawn let
                _h0 = NDHistogram(
                    "feaset.len" => 0:1000,
                    "ko.indx" => 0:1000,
                )
                bb_count = 0
                for bb in bbch
                    @show (bb.id, threadid())
                    islocked(bb) && continue

                    for b in bb
                        get(b, "flags", "duplicate.flag", true) && continue
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
        
        @show length(tasks)
        return merge(map(fetch, tasks))

    end # blobio!
    @show h0_ref
    return nothing
end


## --.-...- --. -. - -.-..- -- .-..- -. -. 
#MARK: # fba.v hist
let
    # params
    blep0 = G["gen.net0", "net0.blep0.ref"][]
    net0 = G["gen.net0", "net0.ref"][]

    _INDEX_MAP = Dict(
        eid => colindex(blep0, extras(net0, eid))
        for eid in [
                "EX_O2", "EX_CO2", "EX_GLC", "EX_NH4", "EX_GLU", "BIOM", "ATPM"
            ] if hascolid(blep0, extras(net0, eid))
    )

    function _getsol(sol, id)
        try
            return sol[_INDEX_MAP[id]]
        catch err
            return 0.0
        end
    end

    n_tasks = min(nthreads(), 40)
    ch_size = min(nthreads(), 40)

    h0_ref = blobio!(C, 
        "fba.sol.hist", 
        "h0", 
        :setser!
    ) do
        bbch = eachbatch(B, "fba.feasures"; 
            n_tasks = 1, ch_size
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
                        
                        sol = b["cargo.fba", "sol"]::Vector{Float64}
                        isempty(sol) && continue

                        InCmol = abs(6 * sol[_INDEX_MAP["EX_GLC"]])

                        try
                            count!(_h0, (
                                    _getsol(sol, "EX_O2"),
                                    _getsol(sol, "EX_CO2"),
                                    _getsol(sol, "EX_GLC"),
                                    _getsol(sol, "EX_NH4"),
                                    _getsol(sol, "EX_GLU"),
                                    _getsol(sol, "BIOM"),
                                    _getsol(sol, "ATPM"),
                                    _getsol(sol, "EX_O2")/InCmol,
                                    _getsol(sol, "EX_CO2")/InCmol,
                                    _getsol(sol, "EX_GLC")/InCmol,
                                    _getsol(sol, "EX_NH4")/InCmol,
                                    _getsol(sol, "EX_GLU")/InCmol,
                                    _getsol(sol, "BIOM")/InCmol,
                                    _getsol(sol, "ATPM")/InCmol,
                                ), 
                                1
                            )
                        catch err;
                            println("\n", sprint(showerror, err, catch_backtrace()))
                        end
                    end # for b 
                    bb_count += 1
                    bb_count < Inf || break
                end # for bb
                return _h0
            end # @spawn let
        end # tasks = map

        return merge(map(fetch, tasks))
        
    end # blobio!

    @show h0_ref
    return nothing
end

