module Entropy
using HDF5, JLD
include("types.jl")
include("nsb2.jl")
include("utils.jl")
include("bias.jl")
include("hist.jl")
using .Histograms
include("binless.jl")
include("binned.jl")
include("algorithms.jl")
include("process.jl")
import GUICheck
if GUICheck.hasgui()
	include("plot.jl")
end

export Hist2d

end # module
