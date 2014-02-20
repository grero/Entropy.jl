module Binned
import Base.Test
include("../src/binned.jl")

function entropy()
    srand(1234)
    x = rand(1:5,1000)
    E = entropy(x)
    Base.Test.@test_approx_eq E 2.3202534717146475
    println("Test passed. E = $E")
end

end


