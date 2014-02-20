using Stats
include("nsb.jl")

function entropy{T<:Integer}(x::Array{T,1})
    N = sum(x)
    K = length(x)
    P = x/N #convert to probabilitiy
    return -dot(P,log2(P + (P.==0)))
end

#functions operating on tress
function entropy(X::Dict{Int64, Dict{Int64, Dict{Int64,Int64}}})
    E = 0
    N = 0
    for (k,v) in X
        ee,n = entropy(v)
        E += ee*n
        N += n
    end
    return E/N,N
end

function entropy(X::Dict{Int64, Dict{Int64,Int64}})
    E = 0
    N = 0
    for (k,v) in X
        ee,n = entropy(v)
        E += ee*n
        N += n
    end
    return E/N,N
end

function entropy(X::Dict{Int64,Int64})
    N = sum(values(X))
    P = float(map(v -> v/N, values(X)))
    return (-sum(P.*log2(P + float(P.==0))),N)
end
