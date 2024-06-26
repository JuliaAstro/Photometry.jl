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

# See Also
* [Source Detection Algorithms](@ref)
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
"""
@with_kw struct PeakMesh <: SourceFinder
    box_size::NTuple{2,<:Integer} = (3, 3)
    nsigma::Float64 = 3
    PeakMesh(box_size::NTuple{2,<:Integer}, nsigma) = new(box_size, nsigma)
    PeakMesh(box_size::Integer, nsigma) = new((box_size, box_size), nsigma)
end

function extract_sources(alg::PeakMesh, data::AbstractMatrix{T}, error = zeros(T, size(data)), sort = true) where T
    threshold = @. error * alg.nsigma
    data_max = mapwindow(maximum, data, alg.box_size, border = Fill(zero(T)))

    rows = NamedTuple{(:x, :y, :value),Tuple{Int,Int,T}}[]
    @inbounds for idx in CartesianIndices(data)
        if data[idx] > threshold[idx] && data[idx] == data_max[idx]
            push!(rows, (x = idx.I[2], y = idx.I[1], value = data[idx]))
        end
    end
    sort && sort!(rows, by = row->row.value, rev = true)
    return Table(rows)
end

end # module Detection
