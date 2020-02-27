using Statistics

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
