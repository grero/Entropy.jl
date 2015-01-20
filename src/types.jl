
abstract ShannonEntropy
abstract GroupedShannonEntropy <: ShannonEntropy
abstract Hist

type Hist1d <: Hist
    wordcnt::Array{Float64,1}
    bins::Array{Float64,1}
end

type Hist2d <: Hist
    wordcnt::Array{Float64,2}
    binsx::Array{Float64,1}
    binsy::Array{Float64,1}
end

type BinnedEntropy <: ShannonEntropy
    counts::Hist
    E::Float64
end

type BinnedInformation <: ShannonEntropy
    counts::Hist
    I::Float64
    function BinnedInformation(NN::Hist2d)
        Nc = NN.wordcnt
        N = sum(Nc[:])
        H1 = log2(Nc + (Nc.==0))
        H2 = log2(sum(Nc,1))
        H3 = log2(sum(Nc,2))
        H4 = log2(N)
        H = Nc.*(broadcast(-,broadcast(-, H1, H2),H3) + H4)
        H[isnan(H)] = 0
        H = 1/N*sum(sum(H))
        new(A,H)
    end
end

type TemporalEntropy <: ShannonEntropy
	H::Array{Float64,1}
	bins::Array{Float64,1}
	word_size::Integer
end

type GroupedTemporalEntropy <: GroupedShannonEntropy
	Hc::Array{Float64,2}
	dHc::Array{Float64,2}
	H::Array{Float64,1}
	dH::Array{Float64,1}
	I::Array{Float64,1}
	dI::Array{Float64,1}
	bins::Array{Float64,1}
	word_size::Integer
	group_labels::Array{Int64,1}
	trials_per_group::Array{Int64,1}
end
type GroupedCountEntropy <: GroupedShannonEntropy
	Hc::Array{Float64,2}
	dHc::Array{Float64,2}
	H::Array{Float64,1}
	dH::Array{Float64,1}
	I::Array{Float64,1}
	dI::Array{Float64,1}
	bins::Array{Float64,1}
	word_size::Integer
	group_labels::Array{Int64,1}
	trials_per_group::Array{Int64,1}
end



function GroupedTemporalEntropy(Hc::Array{Float64,2}, dHc::Array{Float64,2},H::Array{Float64,1},dH::Array{Float64,1}, bins::Array{Float64,1}, word_size::Integer, group_labels::Array{Int64,1},trials_per_group::Array{Int64,1})
	ps = trials_per_group/sum(trials_per_group)
	_hc = zeros(size(Hc,1))
	_dhc = zeros(size(Hc,1))
	for i=1:size(Hc,2)
		for j=1:size(Hc,1)
			_hc[j] += ps[i]*Hc[j,i]
			_dhc[j] += ps[i]*ps[i]*dHc[j,i]*dHc[j,i] #variance
		end
	end
	I = H - _hc
	dI = sqrt(dH.*dH + _dhc)
	return GroupedTemporalEntropy(Hc,dHc,H,dH,I,dI,bins,word_size,group_labels,trials_per_group)
end

function GroupedTemporalEntropy(fname::String;redo::Bool=false)
	GTE = JLD.load(fname,"GTE")
	if redo
		GTE = GroupedTemporalEntropy(GTE.Hc, GTE.dHc, GTE.H, GTE.dH, GTE.bins, GTE.word_size, GTE.group_labels,GTE.trials_per_group)
	end
	return GTE
end

function GroupedCountEntropy(fname::String)
	GCE = JLD.load(fname,"GCE")
	return GCE
end

function GroupedCountEntropy(Hc::Array{Float64,2}, dHc::Array{Float64,2},H::Array{Float64,1},dH::Array{Float64,1}, bins::Array{Float64,1}, word_size::Integer, group_labels::Array{Int64,1},trials_per_group::Array{Int64,1})
	ps = trials_per_group/sum(trials_per_group)
	_hc = zeros(size(Hc,1))
	_dhc = zeros(size(Hc,1))
	for i=1:size(Hc,2)
		for j=1:size(Hc,1)
			_hc[j] += ps[i]*Hc[j,i]
			_dhc[j] += ps[i]*ps[i]*dHc[j,i]*dHc[j,i] #variance
		end
	end
	I = H - _hc
	dI = sqrt(dH.*dH + _dhc)
	return GroupedCountEntropy(Hc,dHc,H,dH,I,dI,bins,word_size,group_labels,trials_per_group)
end

function save(GTE::GroupedTemporalEntropy, fname::String)
	JLD.save(fname, {"GTE" => GTE})
end





