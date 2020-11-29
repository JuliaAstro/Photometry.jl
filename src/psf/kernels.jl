
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

    # hello, it's me: type soup
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
# Gaussian{T}(x, y, fwhm; kwargs...) where {T} = Gaussian(SA[x, y], fwhm; kwargs...)
Gaussian{T}(xy::Tuple, fwhm; kwargs...) where {T} = Gaussian{T}(SVector(xy), fwhm; kwargs...)
# # translate polar coordinates to cartesian, optionally recentering
# Gaussian{T}(p::Polar, fwhm; origin=SA[0, 0], kwargs...) where {T} = Gaussian(CartesianFromPolar()(p) .+ origin, fwhm; kwargs...)

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

Base.size(g::Gaussian) = map(length, g.indices)
Base.axes(g::Gaussian) = g.indices

(g::Gaussian)(dist) = exp(-4 * log(2) * (dist / g.fwhm)^2)

# fallback, also covers scalar case
function Base.getindex(g::Gaussian{T}, idx::Vararg{<:Integer,2}) where T
    Δ = sqeuclidean(SVector(idx), g.pos)
    return convert(T, exp(-4 * log(2) * Δ / g.fwhm^2))
end
# vector case
function Base.getindex(g::Gaussian{T,<:Union{Tuple,AbstractVector}}, idx::Vararg{<:Integer,2}) where T
    weights = SA[1/first(g.fwhm)^2, 1/last(g.fwhm)^2] # manually invert
    Δ = wsqeuclidean(SVector(idx), g.pos, weights)
    return convert(T, exp(-4 * log(2) * Δ))
end

# matrix case
function Base.getindex(g::Gaussian{T,<:AbstractMatrix}, idx::Vararg{<:Integer,2}) where T
    R = SVector(idx) - g.pos
    Δ = R' * ((g.fwhm .^2) \ R)
    return convert(T, exp(-4 * log(2) * Δ))
end

# Alias Normal -> Gaussian
const Normal = Gaussian
@doc (@doc Gaussian) Normal


# """
#     Kernels.Moffat(fwhm)

# A Moffat kernel

# ```math
# K(d) = \\left[1 + \\left(\\frac{d}{\\Gamma/2}\\right)^2 \\right]^{-1}
# ```

# ```
#   ┌────────────────────────────────────────┐
# 1 │                  .::.                  │ Moffat(x)
#   │                  : ::                  │
#   │                 :' :':                 │
#   │                 :  : :.                │
#   │                :   :  :                │
#   │               .:   :  :.               │
#   │               :    :   :               │
#   │              .:    :   '.              │
#   │              :     :    '.             │
#   │             :      :     :.            │
#   │            :       :      :.           │
#   │          .'        :       ':          │
#   │       ..''         :         ':.       │
#   │ ....:''            :            ''.... │
# 0 │''                  :                  '│
#   └────────────────────────────────────────┘
#   -2fwhm                               2fwhm
# ```
# """
# struct Moffat{T} <: PSFKernel
#     fwhm::T
# end

# function (m::Moffat)(dist)
#     hwhm = m.fwhm / 2
#     return inv(1 + (dist / hwhm)^2)
# end

# """
#     Kernels.AiryDisk(fwhm)

# An Airy disk. Guaranteed to work even at `r=0`.

# ```math
# K(r) = \\left[\\frac{2J_1\\left(\\pi r \\right)}{\\pi r}\\right]^2 \\quad r \\approx \\frac{d}{0.97\\Gamma}
# ```

# ```
#   ┌────────────────────────────────────────┐
# 1 │                  .'::                  │ AiryDisk(x)
#   │                 .' :':                 │
#   │                 :  : :                 │
#   │                :'  : ':                │
#   │                :   :  :                │
#   │               .'   :  ':               │
#   │               :    :   :               │
#   │              .'    :   ':              │
#   │              :     :    :              │
#   │             .:     :    '.             │
#   │             :      :     :             │
#   │            .'      :     ':            │
#   │            :       :      :.           │
#   │           :        :       :           │
# 0 │..........:'        :        '..........│
#   └────────────────────────────────────────┘
#   -2fwhm                               2fwhm
# ```
# """
# struct AiryDisk{T} <: PSFKernel
#     fwhm::T
# end

# const rz = 3.8317059702075125 / π

# function (a::AiryDisk)(dist)
#     radius = a.fwhm * 1.18677
#     r = dist / (radius / rz)
#     return iszero(r) ? 1.0 : (2besselj1(π * r) / (π * r))^2
# end

end # module Kernels