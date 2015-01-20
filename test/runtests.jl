import Entropy
import Base.Test

#Testing nsb_entropy

function init()
	srand(1234)
	x = rand(1:100,50)
	n = Entropy.counts(x)
	S_nsb = Entropy.NSBEntropy.NSBEntropy(n,100)
	return S_nsb
end

function test_B_xiK()
	S_nsb = init()
	B =  Entropy..B_xiK(1.2, S_nsb)
	return Base.Test.@test_approx_eq B 1.3853866920075761
end

function test_find_nsb_entropy()
	S_nsb = init()
	Entropy.find_nsb_entropy(S_nsb,1e-5)
	return Base.Test.@test_approx_eq S_nsb.S_nsb 4.484256684719789
end

function test_mlog_evidence()
	S_nsb = init()
	nsb_mlog = Entropy.mlog_evidence(200*S_nsb.K,S_nsb)
	return Base.Test.@test_approx_eq nsb_mlog 46.05806318661056
end

S_nsb = init()
B =  Entropy.NSBEntropy.B_xiK(1.2, S_nsb)
Base.Test.@test_approx_eq B 1.3853866920075761

nsb_mlog = Entropy.mlog_evidence(200*S_nsb.K,S_nsb)
Base.Test.@test_approx_eq nsb_mlog 59.85364833069893

Entropy.find_nsb_entropy(S_nsb,1e-5)
Base.Test.@test_approx_eq S_nsb.S_nsb 4.3952671710926685


