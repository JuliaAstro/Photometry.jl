using Statistics, StatsBase

"""
    Mean <: BackgroundEstimator

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
    Median <: BackgroundEstimator

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
    Mode <: BackgroundEstimator

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
