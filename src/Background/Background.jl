module Background

using Statistics
using ImageFiltering: padarray, Fill, mapwindow

export estimate_background,
       sigma_clip,
       sigma_clip!,
       # Estimators
       MMMBackground,
       SourceExtractorBackground,
       BiweightLocationBackground,
       # RMS Estimators
       StdRMS,
       MADStdRMS,
       BiweightScaleRMS,
       # Interpolators
       ZoomInterpolator,
       IDWInterpolator


###############################################################################
# Abstract types

# Estimators
"""
    Background.LocationEstimator

This abstract type embodies the possible background estimation algorithms for dispatch with [`estimate_background`](@ref).

To implement a new estimator, you must define the struct and define a method like `(::MyEstimator)(data::AbstractArray; dims=:)`.

# See Also
* [Location Estimators](@ref)
"""
abstract type LocationEstimator end

"""
    Background.RMSEstimator

This abstract type embodies the possible background RMS estimation algorithms for dispatch with [`estimate_background`](@ref).

To implement a new estimator, you must define the struct and define a method like `(::MyRMSEstimator)(data::AbstractArray; dims=:)`.

# See Also
* [RMS Estimators](@ref)
"""
abstract type RMSEstimator end


include("estimators.jl")

# Interpolators

"""
    Background.BackgroundInterpolator

This abstract type embodies the different ways of converting a low-resolution mesh into a high-resolution image, especially for dispatch with [`estimate_background`](@ref)

To implement a new interpolation scheme, you must define the struct and define a method like `(::MyInterpolator)(mesh)`

# See Also
* [Interpolators](@ref)
"""
abstract type BackgroundInterpolator end

include("interpolators.jl")

###############################################################################
# Core functions

"""
    estimate_background(data;
        location=SourceExtractorBackground(),
        rms=StdRMS(),
        dims=:)

Perform scalar background estimation using the given estimators.

The value returned will be two values corresponding to the estimated background and the estimated background RMS. The dimensionality will depend on the `dims` keyword.

`location` and `rms` can be anything that is callable, for example `median`, or one of the estimators we provide in [Background Estimators](@ref).

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> bkg, bkg_rms = estimate_background(data)
(1.0, 0.0)

julia> using Statistics: median

julia> bkg, bkg_rms = estimate_background(data; location=median, rms=MADStdRMS())
(1.0, 0.0)
```

# See Also
* [Location Estimators](@ref), [RMS Estimators](@ref)
"""
function estimate_background(data;
        location = SourceExtractorBackground(),
        rms = StdRMS(),
        dims = :)
    return location(data, dims = dims), rms(data, dims = dims)
end

"""
    estimate_background(data, mesh_size;
        location=SourceExtractorBackground(),
        rms=StdRMS(),
        itp=ZoomInterpolator(mesh_size),
        edge_method=:pad,
        [filter_size])

Perform 2D background estimation using the given estimators using meshes.

This function will estimate backgrounds in meshes of size `mesh_size`. When `size(data)` is not an integer multiple of the mesh size, there are two edge methods: `:pad` and `:crop`. The default is to pad (and is recommend to avoid losing image data). If `mesh_size` is an integer, the implicit shape will be square (eg. `mesh_size=4` is equivalent to `mesh_size=(4,4)`).

For evaluating the meshes, each mesh will be passed into `location` to estimate the background and then into `rms` to estimate the background root-mean-square value. These can be anything that is callable, like `median` or one of our [Background Estimators](@ref).

Once the meshes are created they will be median filtered if `filter_size` is given. `filter_size` can be either an integer or a tuple, with the integer being converted to a tuple the same way `mesh_size` is. Filtering is done via `ImageFiltering.MapWindow.mapwindow`. `filter_size` must be odd.

After filtering (if applicable), the meshes are passed to the `itp` to recreate a low-order estimate of the background at the same resolution as the input.

!!! note
    If your `mesh_size` is not an integer multiple of the input size, the output background and rms arrays will not have the same size.

# See Also
* [Location Estimators](@ref), [RMS Estimators](@ref), [Interpolators](@ref)
"""
function estimate_background(data::AbstractArray{T},
        mesh_size::NTuple{2,<:Integer};
        location = SourceExtractorBackground(),
        rms = StdRMS(),
        itp = ZoomInterpolator(mesh_size),
        edge_method = :pad,
        kwargs...) where T

    # handle border effects
    X, nmesh = _craft_array(Val(edge_method), data, mesh_size)

    bkg = zeros(float(T), nmesh)
    bkg_rms = zeros(float(T), nmesh)
    @inbounds for i in 1:nmesh[1], j in 1:nmesh[2]
        # get view of data where mesh is and filter out NaN
        rows = (i - 1) * mesh_size[1] + 1:i * mesh_size[1]
        cols = (j - 1) * mesh_size[2] + 1:j * mesh_size[2]
        d = @view X[rows, cols]
        d_ = @view d[d .!== NaN]
        # skip if only NaN
        length(d_) == 0 && continue
        # calculate background and rms via estimators
        bkg[i, j] = location(d_)
        bkg_rms[i, j] = rms(d_)
    end

    # filtering
    if haskey(kwargs, :filter_size)
        bkg, bkg_rms = _filter(bkg, bkg_rms, kwargs[:filter_size])
    end

    # Now interpolate back to original size
    bkg = itp(bkg)
    bkg_rms = itp(bkg_rms)

    return bkg, bkg_rms
end

estimate_background(data,
    mesh_size::Int;
    location = SourceExtractorBackground(),
    rms = StdRMS(),
    itp = ZoomInterpolator(mesh_size),
    edge_method = :pad, kwargs...) = estimate_background(data, (mesh_size, mesh_size); location = location, rms = rms, itp = itp, edge_method = edge_method, kwargs...)

# pad array
function _craft_array(::Val{:pad}, data, mesh_size)
    nextra = size(data) .% mesh_size
    # have to make sure to avoid padding an extra mesh_size if nextra is 0
    npad = Tuple([n == 0 ? 0 : sz - n for (sz, n) in zip(mesh_size, nextra)])
    X = padarray(float(data), Fill(NaN, (0, 0), npad))
    nmesh = size(X) .÷ mesh_size
    return X, nmesh
end

# crop array
function _craft_array(::Val{:crop}, data, mesh_size)
    nmesh = size(data) .÷ mesh_size
    maxidx = nmesh .* mesh_size
    idxs = Base.OneTo.(maxidx)
    X = data[idxs...]
    return X, nmesh
end

# filter meshes
function _filter(bkg, bkg_rms, filter_size::NTuple{2,<:Integer})
    # skip trivial
    filter_size == (1, 1) && return bkg, bkg_rms

    # perform median filter
    bkg = mapwindow(median!, bkg, filter_size, border = "replicate")
    bkg_rms = mapwindow(median!, bkg_rms, filter_size, border = "replicate")
    return bkg, bkg_rms
end

_filter(bkg, bkg_rms, f::Integer) = _filter(bkg, bkg_rms, (f, f))

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
