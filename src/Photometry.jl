module Photometry

using Reexport

include("aperture/Aperture.jl")
include("background/Background.jl")
include("detection/Detection.jl")

@reexport using .Aperture
@reexport using .Background
@reexport using .Background

end
