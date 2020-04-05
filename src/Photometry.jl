module Photometry

using Reexport

include("Aperture/Aperture.jl")
include("Background/Background.jl")
include("Detection/Detection.jl")

@reexport using .Aperture
@reexport using .Background
@reexport using .Detection

end
