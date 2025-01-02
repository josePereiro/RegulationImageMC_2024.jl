@time begin
    using RegulationImageMC_2024
    using Random
end

# --.-...- --. -. - -.-..- -- .-..- -. -. 
include("0.0_proj.jl")
include("1.99_sim.base.jl")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
let
    # globals blobs
    sortfun = shuffle!
    empty_count = 0
    nonempty_count = 0
    for bb in eachbatch(B, "fba.feasures"; sortfun)
        islocked(bb) && continue
        lock(bb) do
            @show bb.id

            for b in bb
                flag = get!(b, "flags", "sol.empty.flag", false)
                sol = b["cargo.fba", "sol"] 
                flag = flag || isempty(sol)
                b["flags", "sol.empty.flag"] = flag
                flag ? (empty_count += 1) : (nonempty_count += 1)
            end

            # info
            @info("EMPTY", 
                empty_count,
                nonempty_count,
                dips_frac = empty_count / (empty_count + nonempty_count)
            )
        end
    end
end