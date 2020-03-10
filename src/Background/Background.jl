module Background

using Statistics

export estimate_background,
       sigma_clip,
       sigma_clip!,
       # Estimators
       Mean,
       Median,
       Mode,
       MMM,
       SourceExtractor,
       BiweightLocation,
       # RMS Estimators
       StdRMS,
       MADStdRMS,
       BiweightScaleRMS


# Abstract types
"""
    Background.BackgroundEstimator

This abstract type embodies the possible background estimation algorithms for dispatch with [`estimate_background`](@ref).

To implement a new estimator, you must define the struct and define a method like `(::MyEstimator)(data::AbstractArray; dims=:)`.
"""
abstract type BackgroundEstimator end

"""
    Background.BackgroundRMSEstimator

This abstract type embodies the possible background RMS estimation algorithms for dispatch with [`estimate_background`](@ref).

To implement a new estimator, you must define the struct and define a method like `(::MyRMSEstimator)(data::AbstractArray; dims=:)`.
"""
abstract type BackgroundRMSEstimator end

# Estimators
include("estimators.jl")

###############################################################################

"""
    estimate_background(data, ::BackgroundEstimator=SourceExtractor, ::BackgroundRMSEstimator=StdRMS; dims=:)

Perform scalar background estimation using the given estimators.

The value returned will be two values corresponding to the estimated background and the estimated background RMS. The dimensionality will depend on the `dims` keyword.

If the background estimator has no parameters (like [`Mean`](@ref)), you can just specify the type without construction.

# See Also
[Background Estimators](@ref)
[Background RMS Estimators](@ref)

# Examples
"""
estimate_background(::AbstractArray, ::BackgroundEstimator = SourceExtractor(), ::BackgroundRMSEstimator = StdRMS(); dims = :)
estimate_background(d::AbstractArray, T::Type{<:BackgroundEstimator} = SourceExtractor, R::Type{<:BackgroundRMSEstimator} = StdRMS; dims = :) = estimate_background(d, T(), R(); dims = dims)

"""
    estimate_background(::BackgroundEstimator, data, mesh_size; edge=:pad, dims=:)

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
    sigma_clip!(x, sigma; center=median(x), std=std(x))
    sigma_clip!(x, sigma_low, sigma_high; center=median(x), std=std(x))

In-place version of [`sigma_clip`](@ref)

!!! warning
    `sigma_clip!` mutates the element in place and mutation cannot lead to change in type.
    Please be considerate of your input type, because if you are using `Int64` and we try to clip it to `0.5` an `InexactError` will be thrown.

    To avoid this, we recommend converting to float before clipping, or using [`sigma_clip`](@ref) which does this internally.
"""
sigma_clip!(x::AbstractArray, sigma_low::Real, sigma_high::Real = sigma_low; center = median(x), std = std(x)) = clamp!(x, center - sigma_low * std, center + sigma_high * std)

"""
    sigma_clip(x, sigma; center=median(x), std=std(x))
    sigma_clip(x, sigma_low, sigma_high; center=median(x), std=std(x))

This function returns sigma-clipped values of the input `x`.

Specify the upper and lower bounds with `sigma_low` and `sigma_high`, otherwise assume they are equal. `center` and `std` are optional keyword arguments which are functions for finding central element and standard deviation.

This will replace values in `x` lower than `center - sigma_low * std` with that value, and values higher than `center + sigma_high * std` with that value.

# Examples
```jldoctest
julia> x = randn(100_000);

julia> extrema(x)
(-4.387579729097121, 4.518192547139076)

julia> x_clip = sigma_clip(x,1);

julia> extrema(x_clip) # should be close to (-1, 1)
(-1.0021043865183705, 1.0011542162690115)
```
"""
sigma_clip(x::AbstractArray, sigma_low::Real, sigma_high::Real = sigma_low; center = median(x), std = std(x)) = sigma_clip!(float(x), sigma_low, sigma_high; center = center, std = std)


end # Background
