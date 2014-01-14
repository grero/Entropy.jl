module Binless
include("../src/binless.jl")
using Distributions

function entropy()
    srand(1234)
    x = randn(1000)
    E = entropy(x)
    Base.Test.@test_approx_eq E 2.0082391501948837
    println("Test passed. E = $E")
end

function entropy2()
    srand(1234)
    MR = MultivariateNormal(zeros(5),eye(5))
    X = rand(MR,1000)
    E = entropy(X)
    Base.Test.@test_approx_eq E 10.211657210460308
    println("Test passed. E = $E")
end

end #Module
