module Background

export estimate_background,
       Mean,
       Median,
       Mode


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


# Estimators
include("stat_estimators.jl")

end # Background
