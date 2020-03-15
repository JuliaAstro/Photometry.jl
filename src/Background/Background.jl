module Background

using Statistics
using ImageFiltering: padarray, Fill

export estimate_background,
       sigma_clip,
       sigma_clip!,
       # Estimators
       MeanBackground,
       MedianBackground,
       ModeBackground,
       MMMBackground,
       SourceExtractorBackground,
       BiweightLocationBackground,
       # RMS Estimators
       StdRMS,
       MADStdRMS,
       BiweightScaleRMS,
       # Interpolators
       ZoomInterpolator


###############################################################################
# Abstract types

# Estimators
"""
    Background.BackgroundEstimator

This abstract type embodies the possible background estimation algorithms for dispatch with [`estimate_background`](@ref).

To implement a new estimator, you must define the struct and define a method like `(::MyEstimator)(data::AbstractArray; dims=:)`.

# See Also
[Location Estimators](@ref)
"""
abstract type BackgroundEstimator end

"""
    Background.BackgroundRMSEstimator

This abstract type embodies the possible background RMS estimation algorithms for dispatch with [`estimate_background`](@ref).

To implement a new estimator, you must define the struct and define a method like `(::MyRMSEstimator)(data::AbstractArray; dims=:)`.

# See Also
[RMS Estimators](@ref)
"""
abstract type BackgroundRMSEstimator end


include("estimators.jl")

# Interpolators

"""
    Background.BackgroundInterpolator

This abstract type embodies the different ways of converting a low-resolution mesh into a high-resolution image, especially for dispatch with [`estimate_background`](@ref)

To implement a new interpolation scheme, you must define the struct and define a method like `(::MyInterpolator)(mesh)`

# See Also
[Interpolators](@ref)
"""
abstract type BackgroundInterpolator end

include("interpolators.jl")


###############################################################################
# Core functions

"""
    estimate_background(data,
        ::BackgroundEstimator=SourceExtractorBackground,
        ::BackgroundRMSEstimator=StdRMS;
        dims=:)

Perform scalar background estimation using the given estimators.

The value returned will be two values corresponding to the estimated background and the estimated background RMS. The dimensionality will depend on the `dims` keyword.

If the background estimator has no parameters (like [`MeanBackground`](@ref)), you can just specify the type without construction.

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> bkg, bkg_rms = estimate_background(data)
(1.0, 0.0)

julia> bkg, bkg_rms = estimate_background(data, BiweightLocationBackground, BiweightScaleRMS)
(1.0, 0.0)
```

# See Also
* [Location Estimators](@ref)
* [RMS Estimators](@ref)
"""
function estimate_background(data,
        bkg::BackgroundEstimator = SourceExtractorBackground(),
        bkg_rms::BackgroundRMSEstimator = StdRMS();
        dims = :)
    return bkg(data, dims = dims), bkg_rms(data, dims = dims)
end
estimate_background(d::AbstractArray, T::Type{<:BackgroundEstimator}, R::Type{<:BackgroundRMSEstimator}; dims = :) = estimate_background(d, T(), R(); dims = dims)
estimate_background(d::AbstractArray, b::BackgroundEstimator, R::Type{<:BackgroundRMSEstimator}; dims = :) = estimate_background(d, b, R(); dims = dims)
estimate_background(d::AbstractArray, T::Type{<:BackgroundEstimator}, r::BackgroundRMSEstimator; dims = :) = estimate_background(d, T(), r; dims = dims)

"""
    estimate_background(data,
        mesh_size,
        ::BackgroundEstimator=SourceExtractorBackground(),
        ::BackgroundRMSEstimator=StdRMS(),
        ::BackgroundInterpolator=ZoomInterpolator(mesh_size);
        edge_method=:pad)

Perform 2D background estimation using the given estimators using meshes.

This function will estimate backgrounds in meshes of size `mesh_size`. When `size(data)` is not an integer multiple of the mesh size, there are two edge methods: `:pad` and `:crop`. The default is to pad (and is recommend to avoid losing image data).

If either size is an integer, the implicit shape will be square (eg. `box_size=4` is equivalent to `box_size=(4,4)`). Contrast this to a single dimension size, like `box_size=(4,)`.

Once the meshes are created they will be passed to the `BackgroundInterpolator` to recreate a low-order estimate of the background at the same resolution as the input.

# See Also
* [Location Estimators](@ref)
* [RMS Estimators](@ref)
* [Interpolators](@ref)
"""
function estimate_background(data::AbstractArray{T},
        mesh_size::NTuple{2,<:Integer},
        BKG::BackgroundEstimator = SourceExtractorBackground(),
        BKG_RMS::BackgroundRMSEstimator = StdRMS(),
        ITP::BackgroundInterpolator = ZoomInterpolator(mesh_size);
        edge_method = :pad) where T

    if edge_method === :pad
        nextra = size(data) .% mesh_size
        # have to make sure to avoid padding an extra mesh_size if nextra is 0
        npad = Tuple([n == 0 ? 0 : sz - n for (sz, n) in zip(mesh_size, nextra)])
        X = padarray(float(data), Fill(NaN, (0, 0), npad))
        nmesh = size(X) .÷ mesh_size
    elseif edge_method === :crop
        nmesh = size(data) .÷ mesh_size
        maxidx = nmesh .* mesh_size
        idxs = Base.OneTo.(maxidx)
        X = data[idxs...]
    else
        error("Invalid edge method: $edge_method")
    end

    bkg = zeros(float(T), nmesh)
    bkg_rms = zeros(float(T), nmesh)
    @inbounds for i in 1:nmesh[1], j in 1:nmesh[2]
        rows = (i - 1) * mesh_size[1] + 1:i * mesh_size[1]
        cols = (j - 1) * mesh_size[2] + 1:j * mesh_size[2]
        d = @view X[rows, cols]
        d_ = @view d[d .!== NaN]
        length(d_) == 0 && continue
        bkg[i, j] = BKG(d_)
        bkg_rms[i, j] = BKG_RMS(d_)
    end

    # Now interpolate back to original size
    bkg = ITP(bkg)
    bkg_rms = ITP(bkg_rms)

    return bkg, bkg_rms
end

estimate_background(data,
    mesh_size::Int,
    bkg::BackgroundEstimator = SourceExtractorBackground(),
    bkg_rms::BackgroundRMSEstimator = StdRMS(),
    itp::BackgroundInterpolator = ZoomInterpolator(mesh_size);
    edge_method = :pad) = estimate_background(data, (mesh_size, mesh_size), bkg, bkg_rms, itp; edge_method = edge_method)

"""
    sigma_clip!(x, sigma; fill=:clamp, center=median(x), std=std(x))
    sigma_clip!(x, sigma_low, sigma_high; fill=:clamp, center=median(x), std=std(x))

In-place version of [`sigma_clip`](@ref)

!!! warning
    `sigma_clip!` mutates the element in place and mutation cannot lead to change in type.
    Please be considerate of your input type, because if you are using `Int64` and we try to clip it to `0.5` an `InexactError` will be thrown.

    To avoid this, we recommend converting to float before clipping, or using [`sigma_clip`](@ref) which does this internally.
"""
function sigma_clip!(x::AbstractArray{T},
    sigma_low::Real,
    sigma_high::Real = sigma_low;
    fill = :clamp,
    center = median(x),
    std = std(x, corrected = false)) where T
    # clamp
    fill === :clamp && return clamp!(x, center - sigma_low * std, center + sigma_high * std)
    # fill
    mask = center - sigma_low * std .< x .< center + sigma_high * std
    x[.!mask] .= T(fill)
    return x
end

"""
    sigma_clip(x, sigma; fill=:clamp, center=median(x), std=std(x, corrected=false))
    sigma_clip(x, sigma_low, sigma_high; fill=:clamp, center=median(x), std=std(x, corrected=false))

This function returns sigma-clipped values of the input `x`.

Specify the upper and lower bounds with `sigma_low` and `sigma_high`, otherwise assume they are equal. `center` and `std` are optional keyword arguments which are functions for finding central element and standard deviation.

If `fill === :clamp`, this will clamp values in `x` lower than `center - sigma_low * std` and values higher than `center + sigma_high * std`. Otherwise, they will be replaced with `fill`.

# Examples
```jldoctest
julia> x = randn(100_000);

julia> extrema(x)
(-4.387579729097121, 4.518192547139076)

julia> x_clip = sigma_clip(x, 1);

julia> extrema(x_clip) # should be close to (-1, 1)
(-1.0021043865183705, 1.0011542162690115)
```
"""
sigma_clip(x::AbstractArray{T}, sigma_low::Real, sigma_high::Real = sigma_low; fill = :clamp, center = median(x), std = std(x)) where T = sigma_clip!(float(x), sigma_low, sigma_high; fill = fill, center = center, std = std)


end # Background
