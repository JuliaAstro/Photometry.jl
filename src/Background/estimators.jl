using Statistics
using StatsBase
using Compat

"""
    Mean

This estimator returns the mean of the input.

# Example
```jldoctest
julia> data = ones(5, 5);

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
julia> data = ones(5 ,5);

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
julia> data = ones(5, 5);

julia> estimate_background(Mode, data)
1.0

julia> estimate_background(Mode, data, dims=1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```
"""
struct Mode <: BackgroundEstimator end

estimate_background(::Mode, data; dims = :) = dims isa Colon ? mode(data) : @compat mode.(eachslice(data, dims = dims))

"""
    SourceExtractor

This estimator returns the background of the input using the SourceExtractor algorithm.

The background is calculated using a mode estimator of the form `(2.5 * median) - (1.5 * mean)`.

If `(mean - median) / std > 0.3` then the median is used and if `std = 0` then the mean is used.

# Example
```jldoctest
julia> data = ones(5 ,5);

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
    MMMBackground(median_factor=3, mean_factor=2)

Estimate the background using a mode estimator of the form `median_factor * median - mean_factor * mean`.
This algorithm is based on the `MMM` routine originally implemented in DAOPHOT. `MMM` uses factors of `median_factor=3` and `mean_factor=2` by default.
This estimator assumes that contaminated sky pixel values overwhelmingly display positive departures from the true value.


# Example
```jldoctest
julia> x = ones(5,5);

julia> estimate_background(MMMBackground, x)
1.0

julia> estimate_background(MMMBackground(4,3), x, dims = 1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```

# See Also
[`SourceExtractor`](@ref)
"""
struct MMMBackground{T <: Number} <: BackgroundEstimator
    median_factor::T
    mean_factor::T
end

MMMBackground() = MMMBackground(3, 2)
estimate_background(alg::MMMBackground, data; dims = :) = alg.median_factor * median(data, dims = dims) - alg.mean_factor * mean(data, dims = dims)
