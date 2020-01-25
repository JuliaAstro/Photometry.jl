module AperturePhotometry

export area,
       mask,
       Exact,
       Center,
       Subpixel

abstract type Aperture end

abstract type OverlapMethod end
struct Exact <: OverlapMethod end
struct Center <: OverlapMethod end
struct Subpixel <: OverlapMethod 
    n::Integer
end

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

include("circular.jl")

end
