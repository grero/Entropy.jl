#Kozachenko-Leonenko estimator
import NearestNeighbors
import Distances
const NN = NearestNeighbors

function nn_entropy{T<:Real}(X::Array{T,2})
	dims = size(X,1)
	N = size(X,2)
    S = spherical_volume(dims)
	t = NN.NaiveNeighborTree(X,NN.Euclidean())
	ee = 0.0
	for i=1:N
		didx,dd = NN.nearest(t,X[:,i],1,i)
		ee += log2(dd[end])
	end
    H = (dims/N)*ee + log2(S*(N-1)/dims) + eulergamma/log(2)
end

function ksg_information{T<:Real}(X1::Array{T,2}, X2::Array{T,2};k::Integer=1)
	
	dims = size(X1,1)
	N1 = size(X1,2)
	N2 = size(X2,2)
	N = N1 + N2
	#contruct the full matrix
	X = cat(2,X1,X2)
    S = spherical_volume(dims)

	#compute the joint entropy as well as the conditional entropy using the Kozachenko-Leonenko estimator
	#get the distances
	t = NN.NaiveNeighborTree(X,NN.Euclidean())
	eex = 0.0
	eey = 0.0
	ee = 0.0
	for i=1:N
		didx,dd = NN.nearest(t,X[:,i],k,i)
		ee += log2(dd[end])
		#get the number of points within this distance for x and y separately
		qidx,bdd = NN.inball(t,X[:,i],dd[end],1)
		if i <= N1
			n = sum(qidx .<= N1) #x
			eex += digamma(n+1)
		else
			n = sum(qidx .> N1) #y
			eey += digamma(n+1)
		end
	end
	ee = digamma(k) - (eex/N1 + eey/N2) + digamma(N) #Eq. 8
    H = (dims/N)*ee + log2(S*(N-1)/dims) + eulergamma/log(2)
	return H 
end

function information(::Type{NSBEntropy},X::Array{Int64,2}, S::Array{Int64,1})
	@assert size(X,ndims(X)) == length(S)
	K = (maximum(X)+1)^size(X,1)
	labels = sort(unique(S))
	Hs = zeros(length(labels))
	y = hash(X)
	n = counts(y)
	H = nsb_entropy(n,K)
	for (i,l) in enumerate(labels)
		n = counts(hash(y[S.==l]))
		Hs[i] = nsb_entropy(n,K)
	end
	return H,Hs
end

function GroupedTemporalEntropy(N::Array{Int64,2}, bins::Array{Float64,1}, trial_labels::Array{Int64,1},word_size::Integer)
	ulabels = sort(unique(trial_labels))
	nlabels = length(ulabels)
	H = NaN*zeros(size(N,1)-word_size)
	dH = NaN*zeros(size(N,1)-word_size)
	Hc = NaN*zeros(size(N,1)-word_size,nlabels)
	dHc = NaN*zeros(size(N,1)-word_size,nlabels)
	trials_per_label = zeros(Int64,nlabels)
	for (k,l) in enumerate(ulabels)
		trials_per_label[k] = sum(trial_labels.==l)
	end
	for i=1:size(N,1)-word_size
		y = Entropy.hash(N[i:i+word_size-1,:])
		n = Entropy.counts(y)
		try
			H[i],dH[i] = Entropy.nsb_entropy(n,maximum(y))
		catch ee
			println("Problem computing entropy for bin $i")
		end

		H[i] /= log(2)
		dH[i] /= log(2)
		for (k,l) in enumerate(ulabels)
			n = Entropy.counts(y[trial_labels.==l])
			try
				Hc[i,k],dHc[i,k] = Entropy.nsb_entropy(n,maximum(y))
				Hc[i,k] /= log(2)
				dHc[i,k] /= log(2)
			catch ee
				println("Problem computing entropy for bin $i location $l")
			end
		end
	end
	return GroupedTemporalEntropy(Hc,dHc,H,dH,bins,word_size, ulabels,trials_per_label)
end

function GroupedCountEntropy(N::Array{Int64,2}, bins::Array{Float64,1}, trial_labels::Array{Int64,1},word_size::Integer)
	ntrials = size(N,2)
	ulabels = sort(unique(trial_labels))
	nlabels = length(ulabels)
	H = NaN*zeros(size(N,1)-word_size)
	dH = NaN*zeros(size(N,1)-word_size)
	Hc = NaN*zeros(size(N,1)-word_size,nlabels)
	dHc = NaN*zeros(size(N,1)-word_size,nlabels)
	trials_per_label = zeros(Int64,nlabels)
	for (k,l) in enumerate(ulabels)
		trials_per_label[k] = sum(trial_labels.==l)
	end
	y = zeros(Int64,ntrials)
	for i=1:size(N,1)-word_size
		fill!(y,0)
		for j=1:ntrials
			for k=1:word_size
				y[j] += N[i+k-1,j]
			end
		end
		n = Entropy.counts(y)
		try
			H[i],dH[i] = Entropy.nsb_entropy(n,maximum(y))
		catch ee
			println("Problem computing entropy for bin $i")
		end
		H[i] /= log(2)
		dH[i] /= log(2)
		for (k,l) in enumerate(ulabels)
			n = Entropy.counts(y[trial_labels.==l])
			try
				Hc[i,k],dHc[i,k] = Entropy.nsb_entropy(n,maximum(y))
				Hc[i,k] /= log(2)
				dHc[i,k] /= log(2)
			catch ee
				println("Problem computing entropy for bin $i location $l")
			end
		end
	end
	return GroupedCountEntropy(Hc,dHc,H,dH,bins,word_size, ulabels,trials_per_label)
end
