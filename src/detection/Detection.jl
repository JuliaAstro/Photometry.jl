module Detection

using Parameters
using ImageFiltering
using TypedTables
using FillArrays: Zeros

export PeakMesh, extract_sources

"""
    Detection.SourceFinder

Abstract super type for source detection algorithms used with [`extract_sources`](@ref).
"""
abstract type SourceFinder end


"""
    extract_sources(::SourceFinder, data, [error]; sort=true)

Find and extract point-like sources in `data` using the given
[`Detection.SourceFinder`](@ref) algorithm.

Returns a `TypedTables.Table` with the positions of the sources and any
information specific to the algorithm. For instance, [`PeakMesh`](@ref) adds a
`value` column with the peak values.

The returned `x`/`y` positions index the first/second axis of `data`
respectively, following the [Pixel Convention](@ref) shared with the aperture
types. This allows detected sources to be passed directly to the aperture
constructors, e.g., `CircularAperture.(sources.x, sources.y, r)`.

`data` is assumed to be background-subtracted. If `error` is provided it will be
propagated into the detection algorithm. If `sort` is `true` the sources will
be sorted by their amplitude, otherwise they are returned in array (column-major) order.

`error` should be `nothing` or an `AbstractArray` with the same axes as `data`
defining the expected error in each pixel. If `nothing` is provided, any local
maximum is returned, including negative values. The default is a lazy array of
zeros (`FillArrays.Zeros`), which means only positive pixels are returned.

# See Also
* [Source Detection Algorithms](@ref)

# Example
```jldoctest
julia> data = zeros(20, 30); data[4, 27] = 5; data[10, 3] = 3;

julia> sources = extract_sources(PeakMesh(), data)
Table with 3 columns and 2 rows:
     x   y   value
   ┌──────────────
 1 │ 4   27  5.0
 2 │ 10  3   3.0

julia> photometry(CircularAperture.(sources.x, sources.y, 2.0), data).aperture_sum
2-element Vector{Float64}:
 5.0
 3.0
```
"""
extract_sources


"""
    PeakMesh(box_size=(3, 3), nsigma=3.0)

Detect sources as local peaks above a threshold.

A pixel is extracted when it is a strict local maximum within the box of size
`box_size` centered on it, and its value is greater than the pixel-wise
threshold `error * nsigma` computed in [`extract_sources`](@ref).

`box_size` is the box size along `x` and `y` (the first and second array axes);
a single integer gives a square box. Sizes must be positive, and an even size
is rounded up to the next odd size so that the box stays centered on the pixel.

!!! note
    Only strict local maxima are extracted. A plateau of equal-valued pixels,
    e.g. a saturated or clipped star core, has no strict maximum and is skipped,
    as is any pixel with a `NaN` neighbour.

# Example
```jldoctest
julia> pm = PeakMesh((7, 7), 3.0)
PeakMesh
  box_size: Tuple{Int64, Int64}
  nsigma: Float64 3.0
```
"""
@with_kw struct PeakMesh <: SourceFinder
    box_size::NTuple{2, Int} = (3, 3)
    nsigma::Float64 = 3
    function PeakMesh(box_size::Tuple{Integer, Integer}, nsigma)
        all(>=(1), box_size) || throw(ArgumentError("box_size must be positive, got $box_size"))
        return new(Int.(box_size), nsigma)
    end
    PeakMesh(box_size::Integer, nsigma) = PeakMesh((box_size, box_size), nsigma)
end

function extract_sources(alg::PeakMesh, data::AbstractMatrix{T}, error = Zeros{T}(axes(data)); sort = true) where {T}
    peaks = findlocalmaxima(data; window = alg.box_size)
    if !isnothing(error)
        axes(error) == axes(data) || throw(DimensionMismatch("`error` must have the same axes as `data`, got $(axes(error)) and $(axes(data))"))
        filter!(ci -> data[ci] > alg.nsigma * error[ci], peaks)
    end
    # x indexes the first array axis and y the second (see the Pixel Convention).
    # Keep this the only place that maps between array axes and x/y.
    sm = Table(x = map(ci -> ci[1], peaks), y = map(ci -> ci[2], peaks), value = data[peaks])
    sort && sort!(sm, by = row -> row.value, rev = true)
    return sm
end

end # module Detection
