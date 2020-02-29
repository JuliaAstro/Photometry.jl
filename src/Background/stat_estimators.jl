using Statistics, StatsBase

"""
    sigma_clip(data, sigma_low, sigma_high; center=median, std=std)
    sigma_clip(data, sigma; center=median, std=std)

This function returns sigma clipped values of the input `data`. The funtion `sigma_clip!` is an inplace implementation and mutates the `data`.
`sigma_high` and `sigma_low` are for un-symmetrical clipping, when `sigma_low = sigma_high` then they can be passed as `sigma`.
`center` and `std` are optional parameters which are functions for finding central element and standard deviation.

# Example
```jldoctest
julia> data = [1, 2, 3];

julia> sigma_clip(data, 1)
3-element Array{Float64,1}:
 1.0
 2.0
 3.0

julia> sigma_clip(data, 1, 1)
3-element Array{Float64,1}:
 1.0
 2.0
 3.0
```
"""

function sigma_clip!(data::AbstractArray, sigma_low::Real, sigma_high::Real=sigma_low; center=median, std=std)
    mean = center(data)
    deviation = std(data)
    clamp!(data, mean - sigma_low * deviation, mean + sigma_high * deviation)
    return data
end

sigma_clip(data, rest...; kwargs...) = sigma_clip!(float(copy(data)), rest...; kwargs...)

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
