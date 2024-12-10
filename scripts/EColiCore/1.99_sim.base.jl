using Base.Threads

# ..- .- -. -. -. - . -. --. .-.-......- - 
# Constants/config

# Gurobi
import Gurobi
# GRB_ENV = Gurobi.Env()
# LP_SOLVER = () -> Gurobi.Optimizer(GRB_ENV)

LP_SOLVER = Clp.Optimizer

# Gurobi
NTHREADS = nthreads()

# Bloberia
BLOBS_PER_BATCH = 1000

# ..- .- -. -. -. - . -. --. .-.-......- - 
# Utils

function _isfeasible!(opm, objidx, feath)
    try; optimize!(opm) 
        sol = solution(opm, objidx)
        return sol > feath
        catch err; @show(err)
    end
    return false
end

function _downcount(vec, th = 1.0)
    return count(vec) do el
        el < th
    end
end

function _downfactor_vec!(downvec::Vector, downset::Vector, down_factor)
    downvec .= 1.0
    for idx in downset
        downvec[idx] *= down_factor
    end
    return downvec
end
_downfactor_vec(N::Int, downset::Vector, down_factor) = 
    _downfactor_vec!(ones(N), downset, down_factor)

## --.-...- --. -. - -.-..- -- .-..- -. -. 
# Hash tracker
struct HashTracker
    hash_set::Set{UInt}
    lim::Int
    function HashTracker(n::Int)
        hash_set = Set{UInt}()
        sizehint!(hash_set, n)
        new(hash_set, n)
    end

end

function check_duplicate!(t::HashTracker, h::UInt64)
    h in t.hash_set && return true
    if length(t.hash_set) >= t.lim
        # delete random
        delete!(t.hash_set, rand(t.hash_set))
    end
    push!(t.hash_set, h)
    return false
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
function _dups_tracker(
        script_id, script_ver; 
        dup_buff_size = 1_000_000
    )
    t0 = HashTracker(dup_buff_size)
    for bb in eachbatch(B, script_id)
        get(bb, script_id, "script.ver", "") == script_ver || continue
        t = get(bb, script_id, "dups_buff", nothing)
        isnothing(t) && continue
        length(t0.hash_set) == dup_buff_size && break
        merge!(t0, t)
    end
    return t0
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
import Base.merge!
Base.merge!(t0::HashTracker, t1::HashTracker) = isempty(t1.hash_set) || push!(t0.hash_set, t1.hash_set...)

import Base.length
Base.length(t::HashTracker) = length(t.hash_set)


## --.-...- --. -. - -.-..- -- .-..- -. -. 

_dot_string(a, as...) = join([a; as...], ".")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
nothing