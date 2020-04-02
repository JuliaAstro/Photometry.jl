#=
Part of this work is derived from astropy/photutils and astropy/astropy. The relevant derivations
are considered under a BSD 3-clause license. =#

using Statistics
using StatsBase

########################################################################
# Location estimators

"""
    SourceExtractorBackground()

This estimator returns the background of the input using the SourceExtractorBackground algorithm.

The background is calculated using a mode estimator of the form `(2.5 * median) - (1.5 * mean)`.

If `(mean - median) / std > 0.3` then the median is used and if `std = 0` then the mean is used.

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> SourceExtractorBackground()(data)
1.0

julia> SourceExtractorBackground()(data, dims=1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```
"""
struct SourceExtractorBackground <: LocationEstimator end

""" Utility function for SourceExtractorBackground algorithm"""
function validate_SE(background, _mean, _median, _std)
    _std ≈ 0 && return _mean
    abs(_mean - _median) / _std > 0.3 && return _median
    return background
end

function (::SourceExtractorBackground)(data; dims = :)
    _mean = mean(data, dims = dims)
    _median = median(data, dims = dims)
    _std = std(data, dims = dims)

    background = @. 2.5 * _median - 1.5 * _mean
    return validate_SE.(background, _mean, _median, _std)
end

"""
    MMMBackground(median_factor=3, mean_factor=2)

Estimate the background using a mode estimator of the form `median_factor * median - mean_factor * mean`.
This algorithm is based on the `MMMBackground` routine originally implemented in DAOPHOT. `MMMBackground` uses factors of `median_factor=3` and `mean_factor=2` by default.
This estimator assumes that contaminated sky pixel values overwhelmingly display positive departures from the true value.

# Examples
```jldoctest
julia> x = ones(3, 5);

julia> MMMBackground()(x)
1.0

julia> MMMBackground(median_factor=4, mean_factor=3)(x, dims = 1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```

# See Also
* [`SourceExtractorBackground`](@ref)
"""
Base.@kwdef struct MMMBackground{T <: Number} <: LocationEstimator
    median_factor::T = 3
    mean_factor::T = 2
end

(alg::MMMBackground)(data; dims = :) = alg.median_factor .* median(data, dims = dims) .- alg.mean_factor .* mean(data, dims = dims)

"""
    BiweightLocationBackground(c = 6.0, M = nothing)

Estimate the background using the robust biweight location statistic.

``\\xi_{biloc}=M + \\frac{\\sum_{|u_i|<1}{(x_i - M)(1 - u_i^2)^2}}{\\sum_{|u_i|<1}{(1-u_i^2)^2}}``

``u_i = \\frac{(x_i - M)}{c\\cdot\\text{MAD}(x)}``

Where ``\\text{MAD}(x)`` is median absolute deviation of `x`.

# Examples
```jldoctest
julia> x = ones(3,5);

julia> BiweightLocationBackground()(x)
1.0

julia> BiweightLocationBackground(c=5.5)(x; dims = 1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```
"""
Base.@kwdef struct BiweightLocationBackground{T <: Number} <: LocationEstimator
    c::T = 6.0
    M::Union{Nothing,T} = nothing
end

# Consider PR in StatsBase.jl for biweight statistics
function biweight_location(data::AbstractArray, c , M)
    M = M === nothing ? median(data) : M
    MAD = mad(data, normalize = false)
    MAD ≈ 0 && return M
    u = @. (data - M) / (c * MAD)

    num = den = zero(eltype(u))
    for ui in u
        abs(ui) < 1 || continue
        num += ui * (1 - ui^2)^2
        den += (1 - ui^2)^2
    end

    den ≈ 0 && return M
    return M + (c * MAD * num) / den
end

_biweight_location(data::AbstractArray, c, M, ::Colon) = biweight_location(data, c, M)
_biweight_location(data::AbstractArray, c, M, dims) = mapslices(X->biweight_location(X, c, M), data, dims = dims)

(alg::BiweightLocationBackground)(data; dims = :) = _biweight_location(data, alg.c, alg.M, dims)


########################################################################
# RMS Estimators

"""
    StdRMS()

Uses the standard deviation statistic for background RMS estimation.

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> StdRMS()(data)
0.0

julia> StdRMS()(data, dims=1)
1×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
```
"""
struct StdRMS <: RMSEstimator end

(::StdRMS)(data; dims = :) = std(data; corrected = false, dims = dims)

"""
    MADStdRMS()

Uses the standard median absolute deviation (MAD) statistic for background RMS estimation.

This is typically given as

``\\sigma \\approx 1.4826 \\cdot \\text{MAD}``

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> MADStdRMS()(data)
0.0

julia> MADStdRMS()(data, dims=1)
1×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
```
"""
struct MADStdRMS <: RMSEstimator end

_mad(data, ::Colon) = mad(data, normalize=true)
_mad(data, dims) = mapslices(x->mad(x, normalize = true), data; dims = dims)
(::MADStdRMS)(data; dims = :) = _mad(data, dims)


"""
    BiweightScaleRMS(c=9.0, M=nothing)

Uses the robust biweight scale statistic for background RMS estimation.

The biweight scale is the square root of the biweight midvariance. The biweight midvariance uses a tuning constant, `c`, and an optional initial guess of the central value `M`.

``\\zeta^2_{biscl}= \\frac{n\\sum_{|u_i|<1}{(x_i - M)^2(1 - u_i^2)^4}}{\\left[\\sum_{|u_i|<1}{(1-u_i^2)(1-5u_i^2)}\\right]^2}``

``u_i = \\frac{(x_i - M)}{c\\cdot\\text{MAD}(x)}``

Where ``\\text{MAD}(x)`` is median absolute deviation of `x`.

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> BiweightScaleRMS()(data)
0.0

julia> BiweightScaleRMS(c=3.0)(data, dims=1)
1×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
```
"""
Base.@kwdef struct BiweightScaleRMS <: RMSEstimator
    c::Number = 9.0
    M::Union{Nothing,Number} = nothing
end

function biweight_scale(x::AbstractArray{T}, c, M) where T
    length(x) == 1 && return zero(T)
    M = M === nothing ? median(x) : M
    _mad = mad(x, normalize = false)
    # avoid divide by zero error
    _mad = _mad ≈ 0 ? 1 : _mad
    u = @. (x - M) / (c * _mad)

    num = den = zero(eltype(u))
        for ui in u
        abs(ui) < 1 || continue
        num += (c * _mad * ui)^2 * (1 - ui^2)^4
        den += (1 - ui^2) * (1 - 5ui^2)
    end

    return sqrt(length(x) * num) / abs(den)
end

_biweight_scale(x::AbstractArray, c, M, ::Colon) = biweight_scale(x, c, M)
_biweight_scale(x::AbstractArray, c, M, dims) = mapslices(x->biweight_scale(x, c, M), x, dims=dims)
(alg::BiweightScaleRMS)(data; dims = :) = _biweight_scale(data, alg.c, alg.M, dims)
