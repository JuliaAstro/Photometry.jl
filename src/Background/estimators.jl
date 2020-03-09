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
    Biweight_location(M=nothing, c=6.0)

Estimate the background using robust statistic biweight location.

# Example
```jldoctest
julia> x = ones(3,5);

julia> estimate_background(Biweight_location, x)
1.0

julia> estimate_background(Biweight_location, x, dims = 1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
 ```
"""
struct Biweight_location{T <: Number} <: BackgroundEstimator
    M::Union{Nothing, T}
    c::T
end

Biweight_location() = Biweight_location(nothing, 6.0)

function Biweight_util(data::Array{T}, M::Union{Nothing, T}, c::T) where T <: Number

    if M == nothing
        M = median(data)
    end

    MAD = mad(data, normalize = false)

    if MAD == 0
        return M
    end

    u = float(data)
    u .-= M
    u ./= (c*MAD)

    num = 0
    den = 0

    for i in eachindex(u)
        if abs(u[i]) < 1
            num += u[i] * (1 - u[i]^2)^2
            den += (1 - u[i]^2)^2
        end
    end
    if den == 0
        return M
    end
    return M + (c * MAD * num)/den
end


estimate_background(alg::Biweight_location, data; dims = :) = dims isa Colon ? Biweight_util(data, alg.M, alg.c) : mapslices(X -> Biweight_util(X, alg.M, alg.c), data, dims=dims)
