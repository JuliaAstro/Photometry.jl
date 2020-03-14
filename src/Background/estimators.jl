using Statistics
using StatsBase

########################################################################
# Location estimators

"""
    Mean

This estimator returns the mean of the input.

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> Mean()(data)
1.0

julia> Mean()(data, dims=1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```
"""
struct Mean <: BackgroundEstimator end

(::Mean)(data; dims = :) = mean(data, dims = dims)

"""
    Median

This estimator returns the median of the input.

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> Median()(data)
1.0

julia> Median()(data, dims=1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```
"""
struct Median <: BackgroundEstimator end

(::Median)(data; dims = :) = median(data, dims = dims)

"""
    Mode

This estimator returns the mode of the input.

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> Mode()(data)
1.0

julia> Mode()(data, dims=1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```
"""
struct Mode <: BackgroundEstimator end

(::Mode)(data; dims = :) = dims isa Colon ? mode(data) : mapslices(mode, data; dims = dims)

"""
    SourceExtractor

This estimator returns the background of the input using the SourceExtractor algorithm.

The background is calculated using a mode estimator of the form `(2.5 * median) - (1.5 * mean)`.

If `(mean - median) / std > 0.3` then the median is used and if `std = 0` then the mean is used.

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> SourceExtractor()(data)
1.0

julia> SourceExtractor()(data, dims=1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```
"""
struct SourceExtractor <: BackgroundEstimator end

""" Utility function for SourceExtractor algorithm"""
function validate_SE(background::Number, _mean::Number, _median::Number, _std::Number)
    _std ≈ 0 && return _mean
    abs(_mean - _median) / _std > 0.3 && return _median
    return background
end

function (::SourceExtractor)(data; dims = :)
    _mean = mean(data, dims = dims)
    _median = median(data, dims = dims)
    _std = std(data, dims = dims)

    background = @. 2.5 * _median - 1.5 * _mean
    return validate_SE.(background, _mean, _median, _std)
end

"""
    MMM(median_factor=3, mean_factor=2)

Estimate the background using a mode estimator of the form `median_factor * median - mean_factor * mean`.
This algorithm is based on the `MMM` routine originally implemented in DAOPHOT. `MMM` uses factors of `median_factor=3` and `mean_factor=2` by default.
This estimator assumes that contaminated sky pixel values overwhelmingly display positive departures from the true value.


# Examples
```jldoctest
julia> x = ones(3, 5);

julia> MMM()(x)
1.0

julia> MMM(4, 3)(x, dims = 1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```

# See Also
[`SourceExtractor`](@ref)
"""
struct MMM{T <: Number} <: BackgroundEstimator
    median_factor::T
    mean_factor::T
end

MMM() = MMM(3, 2)

(alg::MMM)(data; dims = :) = alg.median_factor * median(data, dims = dims) - alg.mean_factor * mean(data, dims = dims)

"""
    BiweightLocation(c = 6.0, M = nothing)

Estimate the background using the robust biweight location statistic.

``\\xi_{biloc}=M + \\frac{\\sum_{|u_i|<1}{(x_i - M)(1 - u_i^2)^2}}{\\sum_{|u_i|<1}{(1-u_i^2)^2}}``

``u_i = \\frac{(x_i - M)}{c\\cdot\\text{MAD}(x)}``

Where ``\\text{MAD}(x)`` is median absolute deviation of `x`.

# Examples
```jldoctest
julia> x = ones(3,5);

julia> BiweightLocation()(x)
1.0

julia> BiweightLocation(5.5)(x; dims = 1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```
"""
struct BiweightLocation{T <: Number} <: BackgroundEstimator
    c::T
    M::Union{Nothing,T}
end

BiweightLocation(c = 6.0) = BiweightLocation(c, nothing)

function biweight_location(data::AbstractArray, c = 6.0, M = median(data))

    M = M === nothing ? median(data) : M

    MAD = mad(data, normalize = false)

    MAD ≈ 0 && return M

    u = @. (data - M) / (c * MAD)

    num = den = zero(eltype(u))

    for ui in u
        if abs(ui) < 1
            num += ui * (1 - ui^2)^2
            den += (1 - ui^2)^2
        end
    end

    den ≈ 0 && return M
    return M + (c * MAD * num) / den
end

(alg::BiweightLocation)(data; dims = :) = dims isa Colon ? biweight_location(data, alg.c, alg.M) :
                                          mapslices(X->biweight_location(X, alg.c, alg.M), data, dims = dims)


########################################################################
# RMS Estimators

"""
    Std

Uses the standard deviation statistic for background RMS estimation.

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> Std()(data)
0.0

julia> Std()(data, dims=1)
1×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
```
"""
struct Std <: BackgroundRMSEstimator end

(::Std)(data; dims = :) = std(data; corrected = false, dims = dims)

"""
    MADStd

Uses the standard median absolute deviation (MAD) statistic for background RMS estimation.

This is typically given as

``\\sigma \\approx 1.4826 \\cdot \\text{MAD}``

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> MADStd()(data)
0.0

julia> MADStd()(data, dims=1)
1×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
```
"""
struct MADStd <: BackgroundRMSEstimator end

(::MADStd)(data; dims = :) = dims isa Colon ? mad(data, normalize = true) :
                                mapslices(x->mad(x, normalize = true), data; dims = dims)


"""
    BiweightScale(c=9.0, M=nothing)

Uses the robust biweight scale statistic for background RMS estimation.

The biweight scale is the square root of the biweight midvariance. The biweight midvariance uses a tuning constant, `c`, and an optional initial guess of the central value `M`.

``\\zeta^2_{biscl}= \\frac{n\\sum_{|u_i|<1}{(x_i - M)^2(1 - u_i^2)^4}}{\\left[\\sum_{|u_i|<1}{(1-u_i^2)(1-5u_i^2)}\\right]^2}``

``u_i = \\frac{(x_i - M)}{c\\cdot\\text{MAD}(x)}``

Where ``\\text{MAD}(x)`` is median absolute deviation of `x`.

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> BiweightScale()(data)
0.0

julia> BiweightScale(3.0)(data, dims=1)
1×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
```
"""
struct BiweightScale <: BackgroundRMSEstimator
    c::Number
M::Union{Nothing,Number}
end

BiweightScale(c = 9.0) = BiweightScale(c, nothing)

function biweight_scale(x::AbstractArray{T}, c = 9.0, M = median(x)) where {T}
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

(alg::BiweightScale)(data; dims = :) = dims isa Colon ? biweight_scale(data, alg.c, alg.M) :
                                          mapslices(x->biweight_scale(x, alg.c, alg.M), data, dims = dims)
