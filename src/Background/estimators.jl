using Statistics
using StatsBase

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
!!! note
    `Mode` does not supports the dims keyword.

# Example
```jldoctest
julia> data = ones(5, 5);

julia> estimate_background(Mode, data)
1.0
```
"""
struct Mode <: BackgroundEstimator end

estimate_background(::Mode, data; dims = :) = mode(data)

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

function estimate_background(::SourceExtractor, data; dims = :)
    _mean = mean(data, dims = dims)
    _median = median(data, dims = dims)
    _std = std(data, dims = dims)

    background = 2.5 * _median - 1.5 * _mean
    if typeof(background) != Float64
        for i in length(background)
            if _std[i] == 0
                background[i] = _mean[i]
            elseif abs(_mean[i] - _median[i])/_std[i] > 0.3
                background[i] = _median[i]
            end
        end
    else
        if _std == 0
            background = _mean
        elseif abs(_mean - _median)/_std > 0.3
            background = _median
        end
    end
    return background
end
