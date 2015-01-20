module NSBEntropy
import GSL
import Base.Test
import StatsBase
using Distributions,StatsBase

s_0(B,K) = polygamma(0,K*B+1) - polygamma(0,B+1) # <S(n=0)>
sigma_0(B,K) = (B+1)/(B*K+1)*polygamma(1,B+1) - polygamma(1,B*K+1)
d_s_0(B,K) = K*polygamma(1,K*B+1) - polygamma(1,B+1)
#s_1(N,K,n) = sum((n+1)./(N+K).*(map(x-> polygamma(1,x+2),n)-polygamma(1,N+K+1))) #<S(n)>/<S(n=0)>
#S_1(N,K,n,B) = s_0(B,K)*s_1(N,K,n) # <S(n)>
#rhoB(B,n) = gamma(BigFloat(length(n)*B))/gamma(BigFloat(sum(n)+length(n)*B)).*prod(map(x->gamma(BigFloat(x+B))/gamma(B),n))
#rho(B,n) = gamma(length(n)*B)/gamma(sum(n)+length(n)*B).*prod(map(x->gamma(x+B)/gamma(B),n))
#I(B,K,N,n) = rho(B,n)*S1(N,K,n,B)*d_s_0(B,K)

I{T<:Integer}(B::Real,n::Array{T,1}) = I(B,n,length(n))

function series{T<:Real}(sz::Integer, a::Array{T,1},x::Real,n::Integer)
	results = 0.0
	for i in sz:-1:1
		results += a[i]*x^(i+n)
	end
	return results
end

function I{T<:Integer}(B::Real,n::Array{T,1},K::Integer)
    N = sum(n)
    r = rho(B,n)
    s = s1(n,B)
    d = d_s_0(B,K) #
    return r*s*d
end

function simplefunc{T<:Real}(n::Array{T,1})
    #gives the entropy, s1 and the variance of the estimate, s2
    #these are valid up to zeroth order in 1/K and 1/N, i.e. for both N and K large
    N = sum(n)
	K1 = sum(n.>=1)
    Delta = N-K1 #number of coincidences
    s1 = (eulergamma - log(2)) + 2*log(N) - polygamma(0,Delta) #Nemenman(2002, Eq. 30)
    s2 = polygamma(1,Delta) #Nemenman(2002, Eq. 31)
    return s1/log(2), s2/log(2)
end
    

function saddlepoint2{T}(n::Array{T,1}, K::Integer)
    N = sum(n)
    K1 = sum(n.>=1) #number of coincidences
    K2 = sum(n.>=2)
    d = (N - K1)/N #relative number of coincidences
	#println(d)
	order = 10
	N2 = N*N
	N3 = N2*N
	N4 = N3*N
	N5 = N4*N
	N6 = N5*N
	N7 = N6*N
	N8 = N7*N
	N9 = N8*N
	N10 = N9*N
	N11 = N10*N
	N12 = N11*N

	Nm = N-1.0
	overN = 1/N
	Nm2 = Nm^2
	Nm3 = Nm^3
	Nm4 = Nm^4
	Nm5 = Nm^5
	Nm6 = Nm^6
	Nm7 = Nm^7
	Nm8 = Nm^8
	Nm9 = Nm^9
	Nm10 = Nm^10

    #b_1 =  (N-1.0)/(2*N)
    #b0 = (-2.0+1.0/N)/3
    #b1 = (N^2.0 - N - 2.0)/(9.0*(N^2 - N))
	b = [b_1 =  (N-1.0)/(2*N),
        (-2.0+1.0/N)/3,
        (N^2.0 - N - 2.0)/(9.0*(N^2 - N)),
	    (2.0*(2.0 - 3.0*N - 3.0*N2 + 2.0*N3))/(135.0*Nm2*N),
	    (4.0*(22.0 + 13.0*N - 12*N2 - 2.0*N3 + N4))/(405.0*Nm3*N),
	    (4.0*(-40.0 + 58.0*N + 65.0*N2 - 40.0*N3 - 5.0*N4 +	
	      2.0*N5))/(1701.0*Nm4*N),
	    (4.0*(-9496.0 - 6912.0*N + 5772.0*N2 + 2251.0*N3 - 
	      1053.0*N4 - 87.0*N5 + 29.0*N6))/(42525.0*Nm5*N),
	    (16.0*(764.0 - 1030.0*N - 1434.0*N2 + 757.0*N3 + 295.0*N4 
	       - 111.0*N5 - 7.0*N6 + 2.0*N7))/(18225.0*Nm6*N),
	    (16.0*(167000.0 + 142516.0*N - 108124.0*N2 - 66284.0*N3 
	       + 26921.0*N4 + 7384.0*N5 - 2326.0*N6 - 116.0*N7 
	       + 29.0*N8))/(382725*Nm7*N),
	    (16.0*(-17886224.0 + 22513608.0*N + 37376676.0*N2 
	       - 17041380.0*N3 - 11384883.0*N4 + 3698262.0*N5 
	       + 846930.0*N6 - 229464.0*N7 - 9387.0*N8 + 2086*N9))/
			(37889775.0*Nm8*N),
	    (16.0*(-4166651072.0 - 3997913072.0*N + 2783482560.0*N2 +
	       2290151964.0*N3 - 803439834.0*N4 - 395614251.0*N5 +
	       108055443.0*N6 + 20215218.0*N7 - 4805712.0*N8 - 165395.0*N9
	       + 33079.0*N10))/(795685275*Nm9*N),
	    (32.0*(52543486208.0 - 62328059360.0*N - 118489458160.0*N2 +
	       47185442088.0*N3 + 44875379190.0*N4 - 12359832987.0*N5 -
	       5400540075.0*N6 + 1272974916.0*N7 + 200644800.0*N8 - 
	       42495955.0*N9 - 1255067.0*N10 +
	       228194.0*N11))/(14105329875.0*Nm10*N)]

	B0 = N*series(order+2, b, (N-K1)/N, -1)
	@assert B0 > 0
	b_1 = b[1]
	b0 = b[2]
	b1 = b[3]
end

function mlog_evidence(n::Array{Int64,1},Β::Real)
    N = sum(n)
    n1 = n[n.>=1] #coincidences
	K1 = length(n1)
	f = 0.0
	if K1 > 0.0
		f += lgamma(n1)
	end
end

function saddlepoint{T}(n::Array{T,1}, K::Integer)
    N = sum(n)
    K1 = sum(n.>=1) #number of coincidences
    K2 = sum(n.>=2)
	b_1 =  (N-1.0)/(2*N)
    b0 =   (-2.0+1.0/N)/3
    b1 =   (N^2.0 - N - 2.0)/(9.0*(N^2 - N))
    d = (N - K1)/N #relative number of coincidences
    k0 = N*(b_1/d + b0 + b1*d)
    idx = find(n.>1)
    k1 = 0.0
    for i=1:length(idx)
        k1 += (polygamma(0,n[idx[i]]) - polygamma(0,1))/(K1/k0^2 - polygamma(1,k0) + polygamma(1,k0 + N))
    end
    k2 = (K1/k0^3 + (polygamma(2,k0) - polygamma(2,k0+N))/2)*k1^2
    for i=1:length(idx)
        k2 += k0*(polygamma(1,n[idx[i]]) - polygamma(1,1))
    end
    k2 /= (K1/k0^2 - polygamma(1,k0) + polygamma(1,k0 + N))
    ks = k0 + k1/K1 + k2/K2
    println("k₀ = $k0, k₁ = $k1, k₂ = $k2, d = $d, b_1 = $b_1, b0 = $b0, b1 = $b1")
    return ks/K #return Beta
end

rho{T<:Integer}(B::Real, n::Array{T,1}) = rho(B,n,length(n))
rho{T<:Integer}(B::Real, n::Array{T,1}) = exp(logrho(B,n,length(n)))

logrho{T<:Integer}(B::Real, n::Array{T,1}) = logrho(B, n,length(n))

function logrho{T<:Integer}(B::Real, n::Array{T,1},K::Integer)
    N = sum(n)
    C1 = GSL.sf_lngamma(B*K) - GSL.sf_lngamma(N+B*K)
    C2 = 0.0
    for i=1:K
        C2 += GSL.sf_lngamma(n[i] + B) - GSL.sf_lngamma(B)
    end
    return C1+C2
end

s1{T<:Integer}(n::Array{T,1},B::Real) = s1(n,B,length(n))

function s1{T<:Integer}(n::Array{T,1},B::Real, K::Integer)
    N = sum(n)
    C = 0.0
    for i=1:K
        C += (n[i]+B)*polygamma(0,n[i]+B+1)
    end
    C -= (N+B*K)*polygamma(0,N+K*B+1)
    return -C/(N+B*K)
end

function s1{T<:Integer}(n::Array{T,1})
    N = sum(n)
    K = length(n)
    C = 0.0
    for i=1:K
        C += (n[i]+1)*polygamma(0,n[i]+2)
    end
    C -= (K+N)*polygamma(0,N+K+1)
    return -C/(N+K)
end

function s2{T<:Integer}(n::Array{T,1})
    N = sum(n)
    K = length(n)
    C = 0.0
    deltapg(k,x,y) = polygamma(k, x) - polygamma(k,y)
    rr = (N+K)*(N+K+1)
    for i=1:K
        C += (n[i]+1)*(n[i] + 2)*(deltapg(0,n[i]+3, N+K+2)^2 + deltapg(1,n[i]+3, N+K+2))/rr
        for j=1:K
            if i != j
                C += (n[i] + 1)*(n[j]+1)*(deltapg(0,n[i] + 2, N+K +2)*deltapg(0,n[j]+1, N+K+2) - polygamma(1,N+K+2))/rr
            end
        end
    end
    return C
end

function s2{T<:Integer}(n::Array{T,1},B::Real)
    N = sum(n)
    K = length(n)
    C = 0.0
    deltapg(k,x,y) = polygamma(k, x) - polygamma(k,y)
    rr = (N+B*K)*(N+B*K+1)
    for i=1:K
        C += (n[i]+B)*(n[i] + B + 1)*(deltapg(0,n[i]+B+2, N+B*K+2)^2 + deltapg(1,n[i]+B+2, N+B*K+2))/rr
        for j=1:K
            if i != j
                C += (n[i]+ B)*(n[j]+B)*(deltapg(0,n[i] +B + 1, N+B*K +2)*deltapg(0,n[j]+B, N+B*K+2) - polygamma(1,N+B*K+2))/rr
            end
        end
    end
    return C
end

function S1{T<:Integer}(B::Real, n::Array{T,1})
    C1 = s1(n,B)
    K = length(n)
    C2 = 1#s_0(B,K)
    return C1*C2
end

entropy{T<:Integer}(n::Array{T,1}) = getEntropy(n,1e-7,3.0,length(n))
entropy{T<:Integer}(n::Array{T,1},bim::Real, bmax::Real) = getEntropy(n,bmin,bmax,length(n))

function getEntropy{T<:Integer}(n::Array{T,1},bmin::Real,bmax::Real,K::Integer)
    #compute the normalization constant
    Z = quadgk(x -> rho(x,n)*d_s_0(x,K),bmin,bmax)
    #compute the integral over S
    Q = quadgk(x -> I(x,n),bmin,bmax)
    E = Q[1]/Z[1]
    #convert from natural logarithm to bits
    E /= log(2)
    dE = sqrt(Q[2]^2 + Z[2]^2)
    return E,dE
end

function test1()
    n = rand(1:5,100)
    q = counts(n) 
    E = S_1(length(n),length(q),q,1.0)
    return E
end

function testEntropy()
    n = rand(1:10,100)
    q = counts(n)
    E = getEntropy(q)
    return E
end

function test_s_0()
    B = 0.1
    K = 100
    ss = s_0(B,K)
    Base.Test.@test_approx_eq ss 2.775507529477798
    println("Test passed. S₀ = $ss")
end

function test_sigma_0()
    B = 0.1
    K = 100
    ss = sigma_0(B,K)
    Base.Test.@test_approx_eq ss 0.04816357939759014
    println("Test passed. ơ₀ = $ss")
end

function test_s1()
    s = s1([0,2])
    Base.Test.@test_approx_eq s 0.45833333333333326
    println("Test passed. s₁ = $s")
end

function test_simplefunc()
    srand(1234)
    n = rand(1:1000,100)
    q = StatsBase.counts(n,1:1000)
    s,ds = simplefunc(q)
    st = log2(1000)
    Base.Test.@test_approx_eq s 12.510509693214221
    Base.Test.@test_approx_eq ds 0.930443179942287
    println("Test passed. s = $s, ds² = $ds. N/K = 0.1. True entropy: $st")
end

end


