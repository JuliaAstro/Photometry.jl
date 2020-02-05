export EllipticalAperture,
       EllipticalAnnulus


struct EllipticalAperture{T <: Number} <: AbstractAperture
    x::T
    y::T
    a::T
    b::T
    theta::T
    cxx::T
    cxy::T
    cyy::T
end

EllipticalAperture(center::AbstractVector, a, b, theta) = EllipticalAperture(center..., a, b, theta)
EllipticalAperture(x, y, a, b, theta) = EllipticalAperture(promote(x, y, a, b, theta)...)

function Base.show(io::IO, e::EllipticalAperture)
    print(io, "EllipticalAperture($(e.x), $(e.y), a=$(e.a), b=$(e.b), theta=$(e.theta))")
end
