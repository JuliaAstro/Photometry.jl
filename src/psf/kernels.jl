
module Kernels

using CoordinateTransformations
using Distances
using SpecialFunctions
using StaticArrays

"""
    Kernels.PSFKernel{T} <: AbstractMatrix{T}

A kernel used for defining a PSF.
"""
abstract type PSFKernel{T} <: AbstractMatrix{T} end

Base.checkbounds(::Type{Bool}, ::PSFKernel, idx...) = true
Base.checkbounds(::Type{Bool}, ::PSFKernel, idx::CartesianIndex) = true


function indices_from_extent(pos, fwhm, maxsize)
    halfextent = @. maxsize * fwhm / 2
    lower = @. floor(Int, pos - halfextent)
    upper = @. ceil(Int, pos + halfextent)
    return (first(lower):first(upper), last(lower):last(upper))
end


# TODO
# function indices_from_extent(pos, fwhm::AbstractMatrix, maxsize)
#     halfextent = maxsize .* fwhm ./ 2
#     lower = @. floor(Int, pos - halfextent)
#     upper = @. ceil(Int, pos - halfextent)
# end

"""
    Kernels.Gaussian(fwhm; maxsize=ceil(Int, 3fwhm))
    Kernels.Normal(fwhm; maxsize=ceil(Int, 3fwhm))

A Gaussian kernel

```math
K(d) = \\exp\\left[-4\\ln{2}\\left(\\frac{d}{\\Gamma}\\right)^2\\right]
```

```
  ┌────────────────────────────────────────┐
1 │                  .'::                  │ Gaussian(x)
  │                 .: :'.                 │
  │                 :  : :                 │
  │                :'  : ':                │
  │                :   :  :                │
  │               .'   :  ':               │
  │               :    :   :               │
  │              .:    :   '.              │
  │              :     :    :              │
  │             .'     :    ':             │
  │             :      :     :             │
  │            .'      :     ':            │
  │           .'       :      ':           │
  │          .'        :       ':          │
0 │.........''         :         ':........│
  └────────────────────────────────────────┘
  -2fwhm                               2fwhm
```
"""
struct Gaussian{T,FT,VT<:AbstractVector,IT<:Tuple} <: PSFKernel{T}
    pos::VT
    fwhm::FT
    indices::IT

    Gaussian{T}(pos::VT, fwhm::FT, indices::IT) where {T,VT<:AbstractVector,FT,IT<:Tuple} = new{T,FT,VT,IT}(pos, fwhm, indices)
end

## constructors
# default type is Float64
Gaussian(args...; kwargs...) = Gaussian{Float64}(args...; kwargs...)
# parse indices from maxsize
Gaussian{T}(pos::AbstractVector, fwhm; maxsize=3) where {T} = Gaussian{T}(pos, fwhm, indices_from_extent(pos, fwhm, maxsize))
# default position is [0, 0]
Gaussian{T}(fwhm; kwargs...) where {T} = Gaussian{T}(SA[0, 0], fwhm; kwargs...)
# # parse position to vector
Gaussian{T}(x::Number, y::Number, fwhm; kwargs...) where {T} = Gaussian{T}(SA[x, y], fwhm; kwargs...)
Gaussian{T}(xy::Tuple, fwhm; kwargs...) where {T} = Gaussian{T}(SVector(xy), fwhm; kwargs...)
# # translate polar coordinates to cartesian, optionally recentering
Gaussian{T}(p::Polar, fwhm; origin=SA[0, 0], kwargs...) where {T} = Gaussian(CartesianFromPolar()(p) .+ origin, fwhm; kwargs...)

Base.size(g::Gaussian) = map(length, g.indices)
Base.axes(g::Gaussian) = g.indices

# fallback, also covers scalar case
function Base.getindex(g::Gaussian{T}, idx::Vararg{<:Integer,2}) where T
    Δ = sqeuclidean(SVector(idx), g.pos)
    val = exp(-4 * log(2) * Δ / g.fwhm^2)
    return convert(T, val)
end
# vector case
function Base.getindex(g::Gaussian{T,<:Union{Tuple,AbstractVector}}, idx::Vararg{<:Integer,2}) where T
    weights = SA[1/first(g.fwhm)^2, 1/last(g.fwhm)^2] # manually invert
    Δ = wsqeuclidean(SVector(idx), g.pos, weights)
    val = exp(-4 * log(2) * Δ)
    return convert(T, val)
end

# matrix case
function Base.getindex(g::Gaussian{T,<:AbstractMatrix}, idx::Vararg{<:Integer,2}) where T
    R = SVector(idx) - g.pos
    Δ = R' * ((g.fwhm .^2) \ R)
    val = exp(-4 * log(2) * Δ)
    return convert(T, val)
end

# Alias Normal -> Gaussian
const Normal = Gaussian
@doc (@doc Gaussian) Normal


############

"""
    Kernels.AiryDisk(pos, fwhm)
"""
struct AiryDisk{T,FT,VT<:AbstractVector,IT<:Tuple} <: PSFKernel{T}
    pos::VT
    fwhm::FT
    indices::IT

    AiryDisk{T}(pos::VT, fwhm::FT, indices::IT) where {T,VT<:AbstractVector,FT,IT<:Tuple} = new{T,FT,VT,IT}(pos, fwhm, indices)
end

## constructors
# default type is Float64
AiryDisk(args...; kwargs...) = AiryDisk{Float64}(args...; kwargs...)
# parse indices from maxsize
AiryDisk{T}(pos::AbstractVector, fwhm; maxsize=3) where {T} = AiryDisk{T}(pos, fwhm, indices_from_extent(pos, fwhm, maxsize))
# default position is [0, 0]
AiryDisk{T}(fwhm; kwargs...) where {T} = AiryDisk{T}(SA[0, 0], fwhm; kwargs...)
# # parse position to vector
AiryDisk{T}(x::Number, y::Number, fwhm; kwargs...) where {T} = AiryDisk{T}(SA[x, y], fwhm; kwargs...)
AiryDisk{T}(xy::Tuple, fwhm; kwargs...) where {T} = AiryDisk{T}(SVector(xy), fwhm; kwargs...)
# # translate polar coordinates to cartesian, optionally recentering
AiryDisk{T}(p::Polar, fwhm; origin=SA[0, 0], kwargs...) where {T} = AiryDisk(CartesianFromPolar()(p) .+ origin, fwhm; kwargs...)

Base.size(a::AiryDisk) = map(length, a.indices)
Base.axes(a::AiryDisk) = a.indices

const rz = 3.8317059702075125 / π

function Base.getindex(a::AiryDisk{T}, idx::Vararg{<:Integer,2}) where T
    radius = a.fwhm * 1.18677
    Δ = euclidean(SVector(idx), a.pos)
    r = Δ / (radius / rz)
    val = ifelse(iszero(r), one(T), 2 * besselj1(π * r) / (π * r))
    return convert(T, val)
end

function Base.getindex(a::AiryDisk{T,<:Union{AbstractVector,Tuple}}, idx::Vararg{<:Integer,2}) where T
    weights = SA[(rz / (first(a.fwhm) * 1.18677))^2, (rz / (last(a.fwhm) * 1.18677))^2]
    r = weuclidean(SVector(idx), a.pos, weights)
    val = ifelse(iszero(r), one(T), 2 * besselj1(π * r) / (π * r))
    return convert(T, val)
end

############

"""
    Kernels.Moffat(pos, fwhm)
"""
struct Moffat{T,FT,VT<:AbstractVector,IT<:Tuple} <: PSFKernel{T}
    pos::VT
    fwhm::FT
    indices::IT

    Moffat{T}(pos::VT, fwhm::FT, indices::IT) where {T,VT<:AbstractVector,FT,IT<:Tuple} = new{T,FT,VT,IT}(pos, fwhm, indices)
end

## constructors
# default type is Float64
Moffat(args...; kwargs...) = Moffat{Float64}(args...; kwargs...)
# parse indices from maxsize
Moffat{T}(pos::AbstractVector, fwhm; maxsize=3) where {T} = Moffat{T}(pos, fwhm, indices_from_extent(pos, fwhm, maxsize))
# default position is [0, 0]
Moffat{T}(fwhm; kwargs...) where {T} = Moffat{T}(SA[0, 0], fwhm; kwargs...)
# # parse position to vector
Moffat{T}(x::Number, y::Number, fwhm; kwargs...) where {T} = Moffat{T}(SA[x, y], fwhm; kwargs...)
Moffat{T}(xy::Tuple, fwhm; kwargs...) where {T} = Moffat{T}(SVector(xy), fwhm; kwargs...)
# # translate polar coordinates to cartesian, optionally recentering
Moffat{T}(p::Polar, fwhm; origin=SA[0, 0], kwargs...) where {T} = Moffat(CartesianFromPolar()(p) .+ origin, fwhm; kwargs...)

Base.size(m::Moffat) = map(length, m.indices)
Base.axes(m::Moffat) = m.indices

# scalar case
function Base.getindex(m::Moffat{T}, idx::Vararg{<:Integer,2}) where T
    hwhm = m.fwhm / 2
    Δ = sqeuclidean(SVector(idx), m.pos)
    val = inv(1 + Δ / hwhm^2)
    return convert(T, val)
end

# vector case
function Base.getindex(m::Moffat{T,<:Union{AbstractVector,Tuple}}, idx::Vararg{<:Integer,2}) where T
    weights = SA[(2 / first(m.fwhm))^2, (2 / last(m.fwhm))^2]
    Δ = wsqeuclidean(SVector(idx), m.pos, weights)
    val = inv(1 + Δ)
    return convert(T, val)
end

end # module Kernels
