module Binned
import Base.Test
import Distributions
import StatsBase
const D = Distributions

using Entropy

function entropy()
    srand(1234)
    x = rand(1:5,1000)
    E = Entropy.entropy(x)
    Base.Test.@test_approx_eq E 2.3202534717146475
    println("Test passed. E = $E")
end

function nsb_entropy(ntrials::Integer,nvars::Integer)
	srand(1234)
	#Poisson variables
	X = rand(D.Poisson(0.2),(ntrials,nvars))
	E_true  = nvars*D.entropy(D.Poisson(0.2))
	M = (maximum(X)+1).^[0:nvars-1]
	N = X*M
	S = StatsBase.counts(N,minimum(N):maximum(N))
	return E_true,Entropy.NSBEntropy.simplefunc(N)	
end

end


