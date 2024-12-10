## --.-...- --. -. - -.-..- -- .-..- -. -. 
# Utils
function _head(vec::Vector{T}, i0) where T
    return i0 < 1 ? T[] : vec[1:i0]
end

function _reset_idxs!(vec, aux, idxs, v0)
    aux .= v0
    @inbounds for i in idxs
        aux[i] = vec[i]
    end
    vec .= aux
    return vec
end

## --.-...- --. -. - -.-..- -- .-..- -. -. 
nothing