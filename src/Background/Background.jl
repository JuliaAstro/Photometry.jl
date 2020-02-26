module Background

# Abstract types
"""
    Background.BackgroundEstimator

This abstract type embodies the possible background estimation algorithms for dispatch with [`estimate_background`](@ref).
"""
abstract type BackgroundEstimator end


"""
    estimate_background(::Type{<:BackgroundEstimator}, data, args...; dims=nothing, kwargs...)
    estimate_background(::BackgroundEstimator, data; dims=nothing)

Perform 2D background estimation using the given estimator. 

The BackgroundEstimator may be constructed before, or if a type is given, any extra arguments and keyword arguments will be passed to the constructor. The value returned will be an array corresponding to the estimated background, whose shape will depend on the `dims` keyword. By default, it will return a background with the same shape as the input `data`.

# See Also
[Background Estimators](@ref)
"""
estimate_background(ALG::Type{<:BackgroundEstimator}, data::AbstractMatrix, args...; dims = nothing, kwargs...) = estimate_background(ALG(args...; kwargs...), data; dims = dims)
estimate_background(::BackgroundEstimator, ::AbstractMatrix; dims = nothing)

end # Background
