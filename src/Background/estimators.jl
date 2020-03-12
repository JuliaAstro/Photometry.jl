using Statistics
using StatsBase

"""
    Mean

This estimator returns the mean of the input.

# Example
```jldoctest
julia> data = ones(3, 5);

julia> estimate_background(Mean, data)
1.0

julia> estimate_background(Mean, data, dims=1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```
"""
struct Mean <: BackgroundEstimator end

estimate_background(::Mean, data; dims = :) = mean(data, dims = dims)

"""
    Median

This estimator returns the median of the input.

# Example
```jldoctest
julia> data = ones(3, 5);

julia> estimate_background(Median, data)
1.0

julia> estimate_background(Median, data, dims=1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```
"""
struct Median <: BackgroundEstimator end

estimate_background(::Median, data; dims = :) = median(data, dims = dims)

"""
    Mode

This estimator returns the mode of the input.

# Example
```jldoctest
julia> data = ones(3, 5);

julia> estimate_background(Mode, data)
1.0

julia> estimate_background(Mode, data, dims=1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```
"""
struct Mode <: BackgroundEstimator end

estimate_background(::Mode, data; dims = :) = dims isa Colon ? mode(data) : mapslices(mode, data; dims = dims)

"""
    SourceExtractor

This estimator returns the background of the input using the SourceExtractor algorithm.

The background is calculated using a mode estimator of the form `(2.5 * median) - (1.5 * mean)`.

If `(mean - median) / std > 0.3` then the median is used and if `std = 0` then the mean is used.

# Example
```jldoctest
julia> data = ones(3, 5);

julia> estimate_background(SourceExtractor, data)
1.0

julia> estimate_background(SourceExtractor, data, dims=1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```
"""
struct SourceExtractor <: BackgroundEstimator end

""" Utility function for SourceExtractor algorithm"""
function validate_SE(background::Number, _mean::Number, _median::Number, _std::Number)
    _std ≈ 0 && return _mean
    abs(_mean - _median) / std > 0.3 && return _median
    return background
end

function estimate_background(::SourceExtractor, data; dims = :)
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


# Example
```jldoctest
julia> x = ones(3, 5);

julia> estimate_background(MMM, x)
1.0

julia> estimate_background(MMM(4, 3), x, dims = 1)
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
estimate_background(alg::MMM, data; dims = :) = alg.median_factor * median(data, dims = dims) - alg.mean_factor * mean(data, dims = dims)

"""
    BiweightLocation(c = 6.0, M = nothing)

Estimate the background using the robust biweight location statistic.

``\\xi_{biloc}=M + \\frac{\\sum_{|u_i|<1}{(x_i - M)(1 - u_i^2)^2}}{\\sum_{|u_i|<1}{(1-u_i^2)^2}}``
with
``u_i = \\frac{(x_i - M)}{c\\cdot\\text{MAD}(x)}``
here ``\\text{MAD}(x)`` is median absolute deviation of `x`.

# Example
```jldoctest
julia> x = ones(3,5);

julia> estimate_background(BiweightLocation, x)
1.0

julia> estimate_background(BiweightLocation(5.5), x; dims = 1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
 ```
"""
struct BiweightLocation{T <: Number} <: BackgroundEstimator
    c::T
    M::Union{Nothing, T}
end

BiweightLocation(c) = BiweightLocation(c, nothing)
BiweightLocation() = BiweightLocation(6.0)

function biweight_location(data::AbstractArray, c, M=median(data))

    M = M === nothing ? median(data) : M

    MAD = mad(data, normalize = false)

    MAD ≈ 0 && return M

    u = @. (data - M) / (c * MAD)

    num = zero(eltype(u))
    den = zero(eltype(u))

    for ui in u
        if abs(ui) < 1
            num += ui * (1 - ui^2)^2
            den += (1 - ui^2)^2
        end
    end

    den ≈ 0 && return M
    return M + (c * MAD * num)/den
end


estimate_background(alg::BiweightLocation, data; dims = :) = dims isa Colon ? biweight_location(data, alg.c, alg.M) : mapslices(X -> biweight_location(X, alg.c, alg.M), data, dims=dims)
