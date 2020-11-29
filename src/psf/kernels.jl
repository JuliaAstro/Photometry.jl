
"""
    Kernels

Statistical kernels for constructing point-spread functions (PSFs). These kernels act like matrices but without allocating any memory, which makes them efficient to fit and apply.

# Kernels

The following kernels are currently implemented
* [`Kernels.Gaussian`](@ref)
* [`Kernels.AiryDisk`](@ref)
* [`Kernels.Moffat`](@ref)

# Usage

Using the kernels should feel just like an array. In fact, `Kernels.PSFKernel <: AbstractMatrix`. However, no data is stored and no allocations have to be made. In other words, representing the kernels as matrices is merely a convenience, since typically astronomical data is stored in dense arrays.

```jldoctest kernel
julia> k = Kernels.Gaussian(5); # fwhm of 5 pixels, centered at (0, 0)

julia> k[0, 0]
1.0
```
because the kernel is a matrix, it needs to have a size. In this case, the size is `maxsize * FWHM` pixels, centered around the origin, and rounded up. We can see how this alters the indices from a typical `Matrix`

```jldoctest kernel
julia> size(k)
(17, 17)

julia> axes(k)
(-8:8, -8:8)
```

these axes are merely a convenience for bounding the kernel, since they accept any real number as input. 

```jldoctest kernel
julia> k[100, 10000]
0.0
```

By bounding the kernel, we get a cutout which can be applied to arrays with much larger dimensions without having to iterate over the whole matrix

```jldoctest
julia> big_mat = ones(101, 101);

julia> small_kern = Kernels.Gaussian(51, 51, 2); # center of big_mat

julia> ax = map(intersect, axes(big_mat), axes(small_kern))
(48:54, 48:54)

julia> cutout = big_mat[ax...]
7×7 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0

julia> photsum = sum(cutout .* small_kern[ax...]))
4.5322418212890625
```
Nice- we only had to reduce ~50 pixels instead of ~10,000 to calculate the aperture sum, and with some care we could make it allocation-free.

Since the kernels are lazy, that means the type of the output can be specified, as long as it can be converted to from a real number (so no integer types).

```jldoctest
julia> kbig = Kernels.Gaussian{BigFloat}(12);

julia> sum(kbig)
163.07467408408593790971336380361822460116627553361468017101287841796875
```

!!! tip "Tip: Automatic Differentation"
    Forward-mode AD libraries tend to use dual numbers, which can cause headaches getting the types correct. We recommend using the primal vector's element type to avoid these headaches
    ```julia
    # example generative model for position and scalar fwhm
    kernel(X::AbstractVector{T}) where {T} = Kernels.Gaussian{T}(X...)
    ```

# Examples

Here is a brief example which shows how to construct a loss function for fitting a `PSFKernel` to some data.

```julia
data = ones(101, 101)
# generative model
function kernel(X::AbstractVector{T}) where T
    position = X[1:2] # x, y position
    fwhm = X[3:4]     # fwhm_x, fwhm_y
    return Kernels.Gaussian{T}(position, fwhm)
end
# objective function
function loss(X::AbstractVector, data)
    k = kernel(X)
    amp = X[5]
    # l2-distance loss (χ² loss)
    return sum(abs2, data .- amp .* k[axes(data)...])
end

test_params = Float64[51, 51, 10, 10, 1]
loss(test_params, data)

# output

10031.03649468148
```

The objective function can then be used with an optimization library like [Optim.jl](https://github.com/JuliaOpt/Optim.jl) to find best-fitting parameters
```julia
using Optim

# Fit our data using test_params as a starting point
# uses Nelder-Mead optimization
res = optimize(P -> loss(P, data), test_params)

# utilize automatic differentiation (AD) to enable
# advanced algorithms, like Newton's method
res_ad = optimize(P -> loss(P, data), test_params, Newton(); autodiff=:forward)
```
we can see which result has the better loss, and then use the generative model to create a kernel that we can use elsewhere
```julia
best_res = minimum(res) < minimum(res_ad) ? res : res_ad
best_fit_params = Optim.minimizer(best_res)
model = kernel(best_fit_params)
```
"""
module Kernels

using CoordinateTransformations
using Distances
using SpecialFunctions
using StaticArrays

"""
    Kernels.PSFKernel{T} <: AbstractMatrix{T}
"""
abstract type PSFKernel{T} <: AbstractMatrix{T} end

# always inbounds
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

@doc raw"""
    Kernels.Gaussian(fwhm; maxsize=3)
    Kernels.Gaussian(position, fwhm; maxsize=3)
    Kernels.Gaussian(x, y, fwhm; maxsize=3)
    Kernels.Gaussian(::Polar, fwhm; maxsize=3, origin=[0, 0])
    Kernels.Gaussian{T}(args...; kwargs...)

An unnormalized bivariate Gaussian distribution. The position can be specified in `(x, y)` coordinates as a `Tuple`, `AbstractVector`, or as separate arguments. By default the kernel is placed at the origin. The position can also be given as a `CoordinateTransformations.Polar`, optionally centered around `origin`.

The `fwhm` can be a scalar (isotropic), vector/tuple (diagonal), or a matrix (correlated). For efficient calculations, we recommend using [StaticArrys](https://github.com/JuliaArrays/StaticArrays.jl). Here, `maxsize` is a multiple of the fwhm, and can be given as a scalar or as a tuple for each axis.

The output type can be specified, and will default to `Float64`. The amplitude is unnormalized, meaning the maximum value will always be 1. This is distinct from the probability distribution (pdf) of a bivariate Gaussian which assures the kernel *sums* to 1. This means the kernels act like a transmission weighting instead of a probability weighting.

# Functional form
```
f(x̂ | x, FWHM) = exp(-4ln2 * ||x̂ - x|| / FWHM^2)
```
where `x̂` and `x` are position vectors (indices) `||⋅||` represents the square-distance, and `FWHM` is the full width at half-maximum.

The FWHM can be unique for each axis, or include a cross-term. In this case, the functional form becomes
```
f(x̂ | x, Q) = exp(-4ln2 * (x̂ - x)ᵀ Q (x̂ - x))
```
where `Q` is the inverse covariance matrix (or precision matrix). This is equivalent to the inverse of the FWHM matrix after squaring each element.
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
Gaussian{T}(p::Polar, fwhm; origin=SA[0, 0], kwargs...) where {T} = Gaussian{T}(CartesianFromPolar()(p) .+ origin, fwhm; kwargs...)

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

"""
    Kernels.Normal

An alias for [`Kernels.Gaussian`](@ref)
"""
const Normal = Gaussian


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
