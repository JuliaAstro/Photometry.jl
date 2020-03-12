module Background

using Statistics
using Images: padarray, Fill

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
```jldoctest
julia> data = ones(3, 5);

julia> bkg, bkg_rms = estimate_background(data)
1.0, 0.0

julia> bkg, bkg_rms = estimate_background(data, MAD, MADStdRMS)
1.0, 0.0
```
"""
function estimate_background(data, bkg::BackgroundEstimator = SourceExtractor(), bkg_rms::BackgroundRMSEstimator = StdRMS(); dims = :)
    return bkg(data, dims = dims), bkg_rms(data, dims = dims)
end
estimate_background(d::AbstractArray, T::Type{<:BackgroundEstimator}, R::Type{<:BackgroundRMSEstimator}; dims = :) = estimate_background(d, T(), R(); dims = dims)

"""
    estimate_background(data, mesh_size, ::BackgroundEstimator=SourceExtractor, ::BackgroundRMSEstimator=StdRMS; edge_method=:pad)

Perform 2D background estimation using the given estimators using meshes.

This function will estimate backgrounds in meshes of size `mesh_size`. When `size(data)` is not an integer multiple of the mesh size, there are two edge methods: `:pad` and `:crop`. The default is to pad (and is recommend to avoid losing image data).

If either size is an integer, the implicit shape will be square (eg. `box_size=4` is equivalent to `box_size=(4,4)`). Contrast this to a single dimension size, like `box_size=(4,)`.

If the background estimator has no parameters (like [`Mean`](@ref)), you can just specify the type without construction.

# See Also
[Background Estimators](@ref)
[Background RMS Estimators](@ref)
"""
function estimate_background(data, 
    mesh_size::NTuple{2,<:Integer}, 
    BKG::BackgroundEstimator = SourceExtractor(), 
    BKG_RMS::BackgroundRMSEstimator = StdRMS(); 
    edge_method = :pad)

    if edge_method === :pad
        nextra = size(data) .% mesh_size
        npad = mesh_size .- nextra
        X = padarray(data, Fill(NaN, (0, 0), npad))
        nmesh = size(X) .÷ mesh_size
    elseif edge_method === :crop
        nmesh = size(data) .÷ mesh_size
        maxidx = nmesh .* mesh_size
        idxs = Base.OneTo.(maxidx)
        X = data[idxs...]
    else
        error("Invalid edge method: $edge_method")
    end

    bkg = Array{eltype(data)}(undef, nmesh)
    bkg_rms = Array{eltype(data)}(undef, nmesh)
    @inbounds for i in 1:nmesh[1], j in 1:nmesh[2]
        rows = i * nmesh[1]:i * nmesh[1] + mesh_size[1]
        cols = j * nmesh[2]:j * nmesh[2] + mesh_size[2]
        bkg[i, j] = BKG(X[rows, cols])
        bkg_rms[i, j] = BKG_RMS(X[rows, cols])
    end

    return bkg, bkg_rms
end

estimate_background(data, mesh_size::Int, bkg::BackgroundEstimator = SourceExtractor(), bkg_rms::BackgroundRMSEstimator = StdRMS(); edge_method = :pad) = estimate_background(data, (mesh_size, mesh_size), bkg, bkg_rms; edge_method = edge_method)

estimate_background(data, mesh_size, T::Type{<:BackgroundEstimator}, R::Type{<:BackgroundRMSEstimator}; edge_method = :pad) = estimate_background(data, T(), R(); edge_method = edge_method)


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
