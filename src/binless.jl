import Distance
import Base.Test

function entropy{T<:Real}(D::Dict{Int64,Array{T,2}})
    E = Dict{Int64,Float64}()
    Et = 0 
    S = 0
    Ec = 0 
    N = length(D)
    for (k,v) in D
        n = size(v,2) #get the number of trials
        ee = entropy(v)
        E[k] = ee
        Et += n*ee
        S += n
        Ec -= log2(n)
    end
    Ec += N*log2(S)
    Ec /= S
    Et /= S
    return Et,Ec,E
end

function entropy{T<:Real}(x::Array{T,1})
    n = length(x)
    d = zeros(T,n)
    dd = 0.0
    dm = Inf 
    for i=1:n
        dm = Inf 
        for j=1:n
            if i==j
                continue
            end
            dd = x[i]-x[j]
            dd = sqrt(dd*dd)
            dm = dd < dm ? dd : dm
        end
        d[i] = dm
    end
    H = (1/n)*sum(log2(d)) + log2(2*(n-1)) + eulergamma/log(2)
    return H
end

function entropy{T<:Real}(X::Array{T,2})
    #find the minimum distance between pairs of columns
    ndims = size(X,1)
    N = size(X,2)
    d = Distance.pairwise(Distance.Euclidean(),X)
    #find the minimum distance for each pair
    d[d.==0] = Inf
    md = minimum(d,2)
    S = spherical_volume(ndims)
    H = (ndims/N)*sum(log2(md)) + log2(S*(N-1)/ndims) + eulergamma/log(2)
    return H
end

function spherical_volume(r::Real)
    return r*pi^(r/2)/gamma(r/2+1)
end
