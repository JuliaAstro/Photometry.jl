using Statistics

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

This estimator returns the median of the input, after sigma clipping. The clipping on lower side is `sigma_low` and on higher side is `sigma_high`.

#Example
```jldoctest
julia> data = ones(5 ,5)

julia> estimate_background(Median(1, 1), data)
1.0

julia> estimate_background(Median(1,1), data, dims=1)
1×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
```
"""
struct Median{T <: Number} <: BackgroundEstimator
    sigma_high::T
    sigma_low::T

    function Median(sigma_low::T, sigma_high::T) where T <: Number
        sigma_low < 0 && error("Invalid sigma sigma_low=$sigma_low. sigma_low must be greater than or equal to 0")
        sigma_high < 0 && error("Invalid sigma sigma_high=$sigma_high. sigma_high must be greater than or equal to 0")
        new{T}(sigma_low, sigma_high)
    end
end

Median(sigma) = Median(sigma, sigma)
estimate_background(center::Median, data; dims = :, iterations=5) = median(sigma_clip(data, center.sigma_low, center.sigma_high, iterations = iterations), dims = dims)
