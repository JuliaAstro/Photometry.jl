module AperturePhotometry

export area,
       mask

abstract type Aperture end

"""
    area(::Aperture)

Returns the geometric area of the aperture
"""
area(::Aperture)

"""
    mask(::Aperture)

Return a `BitArray` of the masked aperture values in the minimum bounding box
"""
mask(::Aperture)

"""
    bbox(::Aperture)

Return the (`xlow`, `xhigh`, `ylow`, `yhigh`) bounds for a given Aperture
"""
bbox(::Aperture)

include("overlap.jl")
include("circular.jl")

end
