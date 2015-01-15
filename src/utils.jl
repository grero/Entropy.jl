import Base.hash
import StatsBase.counts

hash{T<:Integer}(X::Array{T,1}) = X

function hash{T<:Integer}(X::Array{T,2})
	ndims = size(X,1)
	base = maximum(X)+1
	x = base.^[0:ndims-1]	
	return X'*x 
end

function counts{T<:Integer}(X::Array{T,2})
	y = hash(X)
	n = StatsBase.counts(y,0:maximum(y))
	return n
end

function counts{T<:Integer}(X::Array{T,1})
	n = StatsBase.counts(X,0:maximum(X))
	return n
end

function stratify(n::Array{Int64,1})
	k = StatsBase.counts(n,1:maximum(n))
	kx = k[k.>0]'
	nx = sort(unique(n))
	nx = nx[nx.>0]'
	return nx,kx
end

