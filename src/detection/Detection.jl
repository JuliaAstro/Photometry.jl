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

Uses `method` to find and extract point-like sources.

Returns a `TypedTables.Table` with positions and information related to the
`method`. For instance, using `PeakMesh` returns a table column for the peak
values.

The returned `x`/`y` positions index the first/second axis of `data`
respectively, following the [Pixel Convention](@ref) shared with the aperture
types. This allows detected sources to be passed directly to the aperture
constructors, e.g., `CircularAperture.(sources.x, sources.y, r)`.

`data` is assumed to be background-subtracted. If `error` is provided it will be
propagated into the detection algorithm. If `sort` is `true` the sources will
be sorted by their amplitude, otherwise they are returned in array (column-major) order.

`error` should be `nothing` or an `AbstractArray` defining the expected error in each pixel.
If `nothing` is provided, any local maximum is returned, including negative values.
The default is `zeros(data)`, which means only positive pixels are returned.

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

Detect sources by finding peaks above a threshold in grids across the image.

This creates a pixel-wise threshold for sources by calculating `error * nsigma`
when used with [`extract_sources`](@ref).
The peaks are found by searching the image in boxes of size `box_size`. If the
maximum value in that box is greater than the threshold set above, the point is
extracted.

# Example
```jldoctest
julia> pm = PeakMesh((7, 7), 3.0)
PeakMesh
  box_size: Tuple{Int64, Int64}
  nsigma: Float64 3.0
```
"""
@with_kw struct PeakMesh <: SourceFinder
    box_size::NTuple{2, <:Integer} = (3, 3)
    nsigma::Float64 = 3
    PeakMesh(box_size::NTuple{2, <:Integer}, nsigma) = new(box_size, nsigma)
    PeakMesh(box_size::Integer, nsigma) = new((box_size, box_size), nsigma)
end

function extract_sources(alg::PeakMesh, data::AbstractMatrix{T}, error = Zeros(data); sort = true) where {T}
    peaks = findlocalmaxima(data; window = alg.box_size)
    if !isnothing(error)
        filter!(ci -> data[ci] > alg.nsigma * error[ci], peaks)
    end
    # x indexes the first array axis and y the second (see the Pixel Convention).
    # Keep this the only place that maps between array axes and x/y.
    sm = Table(map(ci -> (x = ci[1], y = ci[2], value = data[ci]), peaks))
    sort && sort!(sm, by = row -> row.value, rev = true)
    return sm
end

end # module Detection
