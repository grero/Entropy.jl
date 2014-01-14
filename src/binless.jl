import Distance
import Base.Test

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

function test_entropy()
    srand(1234)
    x = randn(1000)
    E = entropy(x)
    Base.Test.@test_approx_eq E 2.0082391501948837
    println("Test passed. E = $E")
end

function compteEntropy(X::Array{Float64,2})
    #find the minimum distance between pairs of columns
    nvars = size(S,1)
    ntrains = size(S,2)
    d = Distance.pairwise(Distance.Euclidian,X)
    #find the minimum distance for each pair
    d[d.==0] = inf
    md = zeros(ntrains)
    for i=1:ntrains
        md[i] = min(d[i])
    end

end
