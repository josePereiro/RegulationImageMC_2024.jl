using Base.Threads

# ..- .- -. -. -. - . -. --. .-.-......- - 
# Constants/config

# Gurobi
# import Gurobi
# GRB_ENV = Gurobi.Env()
# LP_SOLVER = () -> Gurobi.Optimizer(GRB_ENV)

LP_SOLVER = Clp.Optimizer

# Gurobi
NTHREADS = nthreads()

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
# TODO: Add check memory:
# - dup counts
# - nondup counts
struct HashTracker
    hash_set::Set{UInt}
    lim::Int
    function HashTracker(n::Int)
        hash_set = Set{UInt}()
        sizehint!(hash_set, n)
        new(hash_set, n)
    end

end

import Base.push!
function Base.push!(t::HashTracker, h::UInt64)
    h in t.hash_set && return true
    if length(t.hash_set) >= t.lim
        # delete random
        delete!(t.hash_set, rand(t.hash_set))
    end
    push!(t.hash_set, h)
    return false
end

check_duplicate!(t::HashTracker, el::UInt64) = push!(t, el)
    
import Base.merge!
function Base.merge!(t0::HashTracker, t1::HashTracker) 
    isempty(t1) && return
    for el in t1.hash_set
        push!(t0, el)
    end
end

import Base.length
Base.length(t::HashTracker) = length(t.hash_set)

import Base.isempty
Base.isempty(t::HashTracker) = isempty(t.hash_set)

## --.-...- --. -. - -.-..- -- .-..- -. -. 
function _dups_tracker(S; 
        dup_buff_size = 1_000_000
    )

    try; lock(S)

        # globals
        script_ver = S["script_ver"]
        
        # cache
        frame = hashed_id("dups.buff.cache", script_ver)
        
        # ram version
        buff1 = get!(S, frame, "buff") do
            HashTracker(dup_buff_size)
        end

        # disk version
        resetframe!(S, frame)
        buff0 = get!(S, frame, "buff") do
            HashTracker(dup_buff_size)
        end
        merge!(buff0, buff1)

        serialize!(S, frame)
        
        return buff0
    finally
        unlock(S)
    end
    
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 

_dot_string(a, as...) = join([a; as...], ".")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
nothing