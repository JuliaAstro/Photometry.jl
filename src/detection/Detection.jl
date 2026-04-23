module Detection

using Parameters
using ImageFiltering
using TypedTables

export PeakMesh, extract_sources

"""
    Detection.SourceFinder

Abstract super type for source detection algorithms used with [`extract_sources`](@ref).
"""
abstract type SourceFinder end


"""
    extract_sources(::SourceFinder, data, [error]; sorted=true)

Uses `method` to find and extract point-like sources.

Returns a `TypedTables.Table` with positions and information related to the
`method`. For instance, using `PeakMesh` returns a table a column for the peak
values.

`data` is assumed to be background-subtracted. If `error` is provided it will be
propagated into the detection algorithm. If `sorted` is `true` the sources will
be sorted by their amplitude.

`error` can be a scalar or an AbstractArray, defining the expected error in each pixel overall or individually.
If `nothing` is provided, any local maximum is returned. The default is `zero`, which means only positive pixels are returned.

# See Also
* [Source Detection Algorithms](@ref)

# Example
```jldoctest
julia> using Random: Xoshiro

julia> data = rand(Xoshiro(1234), 2048, 2048);

julia> pm = PeakMesh((7, 7), 3.0)
PeakMesh
  box_size: Tuple{Int64, Int64}
  nsigma: Float64 3.0


julia> sources = extract_sources(pm, data)
Table with 3 columns and 85668 rows:
      x     y     value
    ┌─────────────────────
 1  │ 1200  1203  1.0
 2  │ 704   827   1.0
 3  │ 1345  1285  1.0
 4  │ 875   1762  0.999999
 5  │ 298   1017  0.999999
 6  │ 580   508   0.999999
 7  │ 654   236   0.999999
 8  │ 237   1287  0.999999
 9  │ 380   725   0.999998
 10 │ 1113  337   0.999998
 11 │ 596   1397  0.999998
 12 │ 1351  744   0.999997
 13 │ 643   1471  0.999997
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
    box_size::NTuple{2,<:Integer} = (3, 3)
    nsigma::Float64 = 3
    PeakMesh(box_size::NTuple{2,<:Integer}, nsigma) = new(box_size, nsigma)
    PeakMesh(box_size::Integer, nsigma) = new((box_size, box_size), nsigma)
end

function extract_sources(alg::PeakMesh, data::AbstractMatrix{T}, error = Zeros(data), sort = true) where T
    sm = findlocalmaxima(data; window = alg.box_size)
    to_nt(ci) = (x=ci[2], y=ci[1], value=data[ci])
    sm = Table(map(to_nt, sm))
    if !(isnothing(error))
        threshold = (error .* alg.nsigma)
        sm = filter(row -> row.value > threshold[row.y, row.x], sm)
    end
    sort && sort!(sm, by = row -> row.value, rev = true)
    return sm
end

end # module Detection
