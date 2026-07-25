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
`method`. For instance, using `PeakMesh` returns a table a column for the peak
values. The returned `x`/`y` positions index the first/second axis of `data` respectively, following the same coordinate convention as the aperture types. This allows for detected sources to be passed directly to the aperture constructors, e.g., `CircularAperture.(sources.x, sources.y, r)`.

`data` is assumed to be background-subtracted. If `error` is provided it will be
propagated into the detection algorithm. If `sort` is `true` the sources will
be sorted by their amplitude.

`error` should be `nothing` or an `AbstractArray` defining the expected error in each pixel.
If `nothing` is provided, any local maximum is returned, including negative values.
The default is `zeros(data)`, which means only positive pixels are returned.

# See Also
* [Source Detection Algorithms](@ref)

# Example
```jldoctest
julia> data = rand(2048, 2048);

julia> pm = PeakMesh((7, 7), 3.0)
PeakMesh
  box_size: Tuple{Int64, Int64}
  nsigma: Float64 3.0

julia> sources = extract_sources(pm, data)
Table with 3 columns and 86004 rows:
      x     y     value
    ┌─────────────────────
 1  │ 1285  1120  1.0
 2  │ 845   1751  1.0
 3  │ 506   1670  1.0
 4  │ 666   1792  1.0
 5  │ 1456  314   0.999999
 6  │ 432   1723  0.999999
 7  │ 322   209   0.999999
 8  │ 334   1872  0.999999
 9  │ 1269  940   0.999999
 10 │ 493   1624  0.999998
 11 │ 1202  436   0.999998
 12 │ 107   363   0.999998
 13 │ 617   1355  0.999998
 14 │ 179   1355  0.999998
 15 │ 165   1916  0.999997
 16 │ 1963  931   0.999997
 17 │ 215   1246  0.999996
 ⋮  │  ⋮     ⋮       ⋮
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
    sm = findlocalmaxima(data; window = alg.box_size)
    to_nt(ci) = (x = ci[1], y = ci[2], value = data[ci])
    sm = Table(map(to_nt, sm))
    if !(isnothing(error))
        threshold = (error .* alg.nsigma)
        sm = filter(row -> row.value > threshold[row.x, row.y], sm)
    end
    sort && sort!(sm, by = row -> row.value, rev = true)
    return sm
end

end # module Detection
