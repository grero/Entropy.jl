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

function I{T<:Integer}(B::Real,n::Array{T,1},K::Integer)
    N = sum(n)
    r = rho(B,n)
    s = s1(n,B)
    d = d_s_0(B,K) #
    return r*s*d
end

function simplefunc{T<:Real}(n::Array{T,1})
    #gives the entropy, s1 and the variance of the estimate, s2
    #these are valid up to zeroth order in 1/K and !/N, i.e. for both N and K large
    N = sum(n)
    K1 = sum(n.>=1)
    Delta = N-K1
    s1 = (eulergamma - log(2)) + 2*log(N) - polygamma(0,Delta)
    s2 = polygamma(1,Delta)
    return s1/log(2), s2/log(2)
end
    

function saddlepoint{T}(n::Array{T,1}, K::Integer)
    N = sum(n)
    K1 = sum(n.>=1)
    K2 = sum(n.>=2)
    d = (N - K1)/N #relative number of conincidences

    b_1 =  (N-1)/(2*N)
    b0 = (-2*N+1)/(3*N)
    b1 = (N^2 - N - 2)/(9*(N^2 - N))

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
    ks = k0 + k1/K + k2/K2
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
    Base.Test.@test_approx_eq s 10.65905105740672
    Base.Test.@test_approx_eq ds 0.2615937290412652
    println("Test passed. s = $s, ds² = $ds. N/K = 0.1. True entropy: $st")
end

end


