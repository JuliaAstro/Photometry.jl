module AperturePhotometry

export area

abstract type Aperture end

"""
    area(::Aperture)

Returns the geometric area of the aperture
"""
area(::Aperture) = 0

include("circular.jl")

end
