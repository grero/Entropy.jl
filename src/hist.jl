module Histograms
import StatsBase

type Hist1d
    P::Integer
    C::Integer
    N::Integer
    wordlist::Array{Int32,1}
    wordcnt::Array{Float64,1}
    entropy::Float64
end

function hist1d{T<:Integer}(X::Array{T,1})
    ntrials = length(X)
    N = float(StatsBase.counts(X, 0:maximum(X)))
    W = int32(sort(unique(X)))
    nW = length(W)
    N = N[N.>0]
    return Hist1d(ntrials,nW,1,W,N,0.)
end

function hist1d{T<:Integer}(X::Array{T,2})
    ntrials = size(X,2)
    mx = maximum(X,2)+1
    M = cumprod(mx)
    M = [1;M[1:end-1]]
    Y = X'*M
    return hist1d(Y)
end

function hist1d{T<:Integer}(X::Array{T,3})
    nbins = size(X,3)
    H = Array(Hist1d,nbins)
    for i=1:nbins
        H[i] = hist1d(squeeze(X[:,:,i],3))
    end
    return H
end

end


