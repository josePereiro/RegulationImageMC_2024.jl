using Base.Threads

# ..- .- -. -. -. - . -. --. .-.-......- - 
# Constants/config

# Gurobi
# import Gurobi
# GRB_ENV = Gurobi.Env()
# LP_SOLVER = () -> Gurobi.Optimizer(GRB_ENV)

import Clp
LP_SOLVER = Clp.Optimizer
import Ipopt
QUAD_LP_SOLVER = Ipopt.Optimizer

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
# TODO: Move to ProjFlows
# TODO: Add check memory:
# - dup counts
# - nondup counts
struct HashTracker
    hash_set::Set{UInt}
    lim::Int
    dup_count::Vector{UInt128}
    function HashTracker(n::Int)
        hash_set = Set{UInt}()
        sizehint!(hash_set, n)
        new(hash_set, n, UInt128[0, 0])
    end

end

import Base.push!
function Base.push!(t::HashTracker, h::UInt64)
    if h in t.hash_set 
        t.dup_count[1] += UInt128(1)
        return true
    end
    
    # TODO: check better caching strategies
    # - Maybe track the most frequent
    #   - Delete the less frequent first
    if length(t.hash_set) >= t.lim
        # delete random
        delete!(t.hash_set, rand(t.hash_set))
    end
    push!(t.hash_set, h)
    
    t.dup_count[2] += UInt128(1)
    return false
end

check_duplicate!(t::HashTracker, el::UInt64) = push!(t, el)

function dup_ratio(t::HashTracker)
    rat = t.dup_count[1] / (t.dup_count[1] + t.dup_count[2])
    return convert(Float64, rat)
end

function check_count(T::DataType, t::HashTracker) 
    return floor(T, t.dup_count[1] + t.dup_count[2])
end
check_count(t::HashTracker) = check_count(Float64, t) 

function dups_count(T::DataType, t::HashTracker)
    return floor(T, t.dup_count[1])
end
dups_count(t::HashTracker) = dups_count(Int, t)

function nondups_count(T::DataType, t::HashTracker)
    return floor(T, t.dup_count[2])
end
nondups_count(t::HashTracker) = nondups_count(Int, t)
    
import Base.merge!
function Base.merge!(t0::HashTracker, t1::HashTracker) 
    isempty(t1) && return
    for el in t1.hash_set
        push!(t0, el)
    end
    t0.dup_count[1] = t1.dup_count[1]
    t0.dup_count[2] = t1.dup_count[2]
    return t0
end

import Base.length
Base.length(t::HashTracker) = length(t.hash_set)

import Base.isempty
Base.isempty(t::HashTracker) = isempty(t.hash_set)

function Base.empty!(t::HashTracker)
    empty!(t.hash_set)
    r.dup_count[1] = UInt128(0)
    r.dup_count[2] = UInt128(0)
    return nothing
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
function _done_tracker!(f!::Function, S; 
        lk = true
    )
    script_id = S["script_id"]
    return mergeblobs!(S, script_id; lk) do rblob, dblob
        # ram
        hist0 = get!(rblob, "hist") do
            Dict()
        end
        # disk
        isnothing(dblob) && return hist0
        buff1 = get!(dblob, "hist") do
            Dict()
        end
        merge!(hist0, buff1)

        # custom
        val = f!(hist0)

        return val
    end
end
_done_tracker!(S; lk = true) = _done_tracker!(identity, S; lk)

function _check_done_count!(S, id, done_target; lk = true)
    _done_tracker!(S; lk) do done_reg
        done_count = get(done_reg, id, -1) 
        @show done_count
        return done_count >= done_target
    end
end

function _up_done_count!(S, id, up = 1; lk = true)
    _done_tracker!(S; lk) do done_reg
        get!(done_reg, id, 0)
        done_reg[id] += up
    end
end


## --.-...- --. -. - -.-..- -- .-..- -. -. 
function _dups_tracker!(S; 
        dup_buff_size = 1_000_000, 
        lk = true
    )

    script_ver = S["script_ver"]
    script_id = S["script_id"]
    frame = hashed_id("dups.buff.cache", script_ver, script_id)
    mergeblobs!(S, frame; lk) do rblob, dblob

        # ram version
        buff0 = get!(rblob, "buff") do
            HashTracker(dup_buff_size)
        end
        
        # disk version
        # if missing on disk return
        isnothing(dblob) && return nothing
        buff1 = get!(dblob, "buff") do
            HashTracker(dup_buff_size)
        end
                
        # merge
        merge!(buff0, buff1)
        
        return nothing
    end
    return S[frame, "buff"]
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 

_dot_string(a, as...) = join([a; as...], ".")

## --.-...- --. -. - -.-..- -- .-..- -. -. 
nothing