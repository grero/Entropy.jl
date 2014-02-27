module Entropy

include("hist.jl")
using .Histograms
include("binless.jl")
include("binned.jl")

export Hist2d

end # module
