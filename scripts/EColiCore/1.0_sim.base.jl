using Base.Threads

# ..- .- -. -. -. - . -. --. .-.-......- - 
# Constants/config

# Gurobi
import Gurobi
GRB_ENV = Gurobi.Env()
LP_SOLVER = () -> Gurobi.Optimizer(GRB_ENV)

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

# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.
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

# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.
function unordered_hash(vec)
    combined_hash = zero(UInt)
    for x in vec
        combined_hash ⊻= hash(x)
    end
    return combined_hash
end

# ..- .- -. -. -. - . -. --. .-.-......- - 
# TODO: Move to ProjFlows
function hashed_id(s::String, args...)
    h0 = hash(0)
    for a in args
        h0 = hash(a, h0)
    end
    return string(s, h0)
end


# ..- .- -. -. -. - . -. --. .-.-......- - 
nothing