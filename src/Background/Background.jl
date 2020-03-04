module Background

using Statistics

export estimate_background,
       Mean,
       Median,
       Mode,
       sigma_clip,
       sigma_clip!


# Abstract types
"""
    Background.BackgroundEstimator

This abstract type embodies the possible background estimation algorithms for dispatch with [`estimate_background`](@ref).
"""
abstract type BackgroundEstimator end


"""
    estimate_background(::BackgroundEstimator, data; dims=:)

Perform 2D background estimation using the given estimator.

The value returned will be an two arrays corresponding to the estimated background, whose dimensionality will depend on the `dims` keyword and the estimator used.

If the background estimator has no parameters (like [`Mean`](@ref)), you can just specify the type without construction.

# See Also
[Background Estimators](@ref)
"""
estimate_background(::BackgroundEstimator, ::AbstractArray; dims = :)
estimate_background(T::Type{<:BackgroundEstimator}, d::AbstractArray; dims = :) = estimate_background(T(), d; dims = dims)

"""
    estimate_background(:BackgroundEstimator, data, box_size, kernel_size; dims=:)

Perform 2D background estimation using the given estimator using meshes and kernels.

This function will estimate backgrounds in meshes of size `box_size`, using a filter kernel of size `kernel_size`. These correspond to the dimension, so for 2D data you could specify (20,) or (20,20) as the box/kernel size, matching with dims=1 for the scalar variant.

If either size is an integer, the implicit shape will be square (eg. `box_size=4` is equivalent to `box_size=(4,4)`). Contrast this to a single dimension size, like `box_size=(4,)`.

If the background estimator has no parameters (like [`Mean`](@ref)), you can just specify the type without construction.

# See Also
[Background Estimators](@ref)
"""
estimate_background(::BackgroundEstimator, ::AbstractArray, ::Tuple, ::Tuple; dims = :) = error("Not implemented!")
estimate_background(T::Type{<:BackgroundEstimator}, d::AbstractArray, b::Tuple, k::Tuple; dims = :) = estimate_background(T(), d, b, k; dims = dims)
estimate_background(alg::BackgroundEstimator, data::AbstractArray, box_size::Integer, kernel_size; dims = :) = estimate_background(alg, data, (box_size, box_size), kernel_size; dims = dims)
estimate_background(alg::BackgroundEstimator, data::AbstractArray, box_size, kernel_size::Integer; dims = :) = estimate_background(alg, data, box_size, (kernel_size, kernel_size); dims = dims)


"""
    sigma_clip!(data, sigma; center=median, std=std)
    sigma_clip!(data, sigma_low, sigma_high; center=median, std=std)

In-place version of [`sigma_clip`](@ref)

!!! warning
    `sigma_clip!` mutates the element in place and mutation cannot lead to change in type.
    User should be careful about using the data-types. E.g.- `x = [1,2,3]`, calling `clamp!(x, 0.5, 0.5)`
    would lead to error because the value of 1 and 3 should have become `Float` from `Int`, but mutation of type is not
    permissible.
"""

function sigma_clip!(data::AbstractArray, sigma_low::Real, sigma_high::Real = sigma_low; center = median, std = std)
    mean = center(data)
    deviation = std(data)
    return clamp!(data, mean - sigma_low * deviation, mean + sigma_high * deviation)
end

"""
    sigma_clip(data, sigma; center=median, std=std)
    sigma_clip(data, sigma_low, sigma_high; center=median, std=std)

This function returns sigma-clipped values of the input `data`.

Specify the upper and lower bounds with `sigma_low` and `sigma_high`, otherwise assume they are equal. `center` and `std` are optional keyword arguments which are functions for finding central element and standard deviation.

This will replace values in `data` lower than `center - sigma_low * std` with that value, and values higher than `center + sigma_high * std` with that value.

# Example
```jldoctest
julia> x = randn(100000);

julia> extrema(x)
(-4.387579729097121, 4.518192547139076)

julia> x_clip = sigma_clip(x,1);

julia> extrema(x_clip) # should be close to (-1, 1)
(-1.0021043865183705, 1.0011542162690115)
```
"""
sigma_clip(data::AbstractArray, sigma_low::Real, sigma_high::Real = sigma_low; center = median, std = std) = sigma_clip!(float(data), sigma_low, sigma_high; center = center, std = std)

# Estimators
include("stat_estimators.jl")

end # Background
