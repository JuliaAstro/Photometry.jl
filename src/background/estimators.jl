#=
Part of this work is derived from astropy/photutils and astropy/astropy. The relevant derivations
are considered under a BSD 3-clause license. =#

using Statistics
using StatsBase
using Parameters

########################################################################
# Location estimators

"""
    SourceExtractorBackground()

This estimator returns the background of the input using the
SourceExtractorBackground algorithm.

The background is calculated using a mode estimator of the form
`(2.5 * median) - (1.5 * mean)`.

If `(mean - median) / std > 0.3` then the median is used and
if `std = 0` then the mean is used.

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> SourceExtractorBackground()(data)
1.0

julia> SourceExtractorBackground()(data, dims=1)
1×5 Matrix{Float64}:
 1.0  1.0  1.0  1.0  1.0
```
"""
struct SourceExtractorBackground <: LocationEstimator end

""" Utility function for SourceExtractorBackground algorithm"""
function validate_SE(background, _mean, _median, _std)
    _std ≈ 0 && return _mean
    abs(_mean - _median) / _std > 0.3 && return _median
    return background
end

function (::SourceExtractorBackground)(data; dims = :)
    _mean = mean(data, dims = dims)
    _median = median(data, dims = dims)
    _std = std(data, dims = dims)

    background = @. 2.5 * _median - 1.5 * _mean
    return validate_SE.(background, _mean, _median, _std)
end

"""
    MMMBackground(median_factor=3, mean_factor=2)

Estimate the background using a mode estimator of the form
`median_factor * median - mean_factor * mean`.
This algorithm is based on the `MMMBackground` routine originally implemented in
DAOPHOT.
`MMMBackground` uses factors of `median_factor=3` and `mean_factor=2` by default.
This estimator assumes that contaminated sky pixel values overwhelmingly display
positive departures from the true value.

# Examples
```jldoctest
julia> x = ones(3, 5);

julia> MMMBackground()(x)
1.0

julia> MMMBackground(median_factor=4, mean_factor=3)(x, dims = 1)
1×5 Matrix{Float64}:
 1.0  1.0  1.0  1.0  1.0
```

# See Also
* [`SourceExtractorBackground`](@ref)
"""
@with_kw struct MMMBackground{T <: Number} <: LocationEstimator
    median_factor::T = 3
    mean_factor::T = 2
end

(alg::MMMBackground)(data; dims = :) =
    alg.median_factor .* median(data, dims = dims) .- alg.mean_factor .* mean(data, dims = dims)

"""
    BiweightLocationBackground(c = 6.0, M = nothing)

Estimate the background using the robust biweight location statistic.

```math
ξ_{biloc} = M + \\frac{∑_{|uᵢ|<1}{(xᵢ - M)(1 - uᵢ²)²}}{∑_{|uᵢ|<1}{(1-uᵢ²)²}}
```

```math
u_i = \\frac{(x_i - M)}{c⋅\\mathrm{MAD}(x)}
```

Where ``\\mathrm{MAD}(x)`` is median absolute deviation of `x`.

# Examples
```jldoctest
julia> x = ones(3,5);

julia> BiweightLocationBackground()(x)
1.0

julia> BiweightLocationBackground(c=5.5)(x; dims = 1)
1×5 Matrix{Float64}:
 1.0  1.0  1.0  1.0  1.0
```
"""
@with_kw struct BiweightLocationBackground{T <: Number} <: LocationEstimator
    c::T = 6.0
    M::Union{Nothing,T} = nothing
end

# Consider PR in StatsBase.jl for biweight statistics
function biweight_location(data::AbstractArray, c , M)
    M = M === nothing ? median(data) : M
    MAD = mad(data, normalize = false)
    MAD ≈ 0 && return M
    u = @. (data - M) / (c * MAD)

    num = den = zero(eltype(u))
    for ui in u
        abs(ui) < 1 || continue
        num += ui * (1 - ui^2)^2
        den += (1 - ui^2)^2
    end

    den ≈ 0 && return M
    return M + (c * MAD * num) / den
end

_biweight_location(data::AbstractArray, c, M, ::Colon) =
    biweight_location(data, c, M)
_biweight_location(data::AbstractArray, c, M, dims) =
    mapslices(X->biweight_location(X, c, M), data, dims = dims)

(alg::BiweightLocationBackground)(data; dims = :) = _biweight_location(data, alg.c, alg.M, dims)


########################################################################
# RMS Estimators

"""
    StdRMS()

Uses the standard deviation statistic for background RMS estimation.

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> StdRMS()(data)
0.0

julia> StdRMS()(data, dims=1)
1×5 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0
```
"""
struct StdRMS <: RMSEstimator end

(::StdRMS)(data; dims = :) = std(data; corrected = false, dims = dims)

"""
    MADStdRMS()

Uses the standard median absolute deviation (MAD) statistic for background RMS estimation.

This is typically given as

``σ ≈ 1.4826 ⋅ \\mathrm{MAD}``

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> MADStdRMS()(data)
0.0

julia> MADStdRMS()(data, dims=1)
1×5 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0
```
"""
struct MADStdRMS <: RMSEstimator end

_mad(data, ::Colon) = mad(data, normalize=true)
_mad(data, dims) = mapslices(x->mad(x, normalize = true), data; dims = dims)
(::MADStdRMS)(data; dims = :) = _mad(data, dims)


"""
    BiweightScaleRMS(c=9.0, M=nothing)

Uses the robust biweight scale statistic for background RMS estimation.

The biweight scale is the square root of the biweight midvariance. The biweight
midvariance uses a tuning constant, `c`, and an optional initial guess of the
central value `M`.

```math
ζ²_{biscl} = \\frac{n ∑_{|uᵢ|<1}{(xᵢ - M)²(1 - uᵢ²)⁴}}{\\left[∑_{|uᵢ|<1}{(1-uᵢ²)(1-5uᵢ²)}\\right]²}
```

```math
uᵢ = \\frac{(xᵢ - M)}{c⋅\\mathrm{MAD}(x)}
```

Where ``\\mathrm{MAD}(x)`` is median absolute deviation of `x`.

# Examples
```jldoctest
julia> data = ones(3, 5);

julia> BiweightScaleRMS()(data)
0.0

julia> BiweightScaleRMS(c=3.0)(data, dims=1)
1×5 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0
```
"""
@with_kw struct BiweightScaleRMS <: RMSEstimator
    c::Number = 9.0
    M::Union{Nothing,Number} = nothing
end

function biweight_scale(x::AbstractArray{T}, c, M) where T
    length(x) == 1 && return zero(T)
    M = M === nothing ? median(x) : M
    _mad = mad(x, normalize = false)
    # avoid divide by zero error
    _mad = _mad ≈ 0 ? 1 : _mad
    u = @. (x - M) / (c * _mad)

    num = den = zero(eltype(u))
    for ui in u
        abs(ui) < 1 || continue
        num += (c * _mad * ui)^2 * (1 - ui^2)^4
        den += (1 - ui^2) * (1 - 5ui^2)
    end

    return sqrt(length(x) * num) / abs(den)
end

_biweight_scale(x::AbstractArray, c, M, ::Colon) = biweight_scale(x, c, M)
_biweight_scale(x::AbstractArray, c, M, dims) = mapslices(x->biweight_scale(x, c, M), x, dims=dims)
(alg::BiweightScaleRMS)(data; dims = :) = _biweight_scale(data, alg.c, alg.M, dims)

@testsnippet estimators begin
    using Photometry.Background: MMMBackground, SourceExtractorBackground, BiweightLocationBackground,
                                 StdRMS, MADStdRMS, BiweightScaleRMS
    using StatsBase: mad, median, std
end

@testitem "background/estimators: Trivial estimator" setup=[estimators] begin
    @testset "Trivial $E"  for E in [MMMBackground, SourceExtractorBackground, BiweightLocationBackground]
        estimator = E()

        @test estimator(ones(1)) == 1

        data = ones(10, 10)

        @test estimator(data) ≈ 1.0
        @test estimator(data, dims = 1) ≈ ones(1, 10)
        @test estimator(data, dims = 2) ≈ ones(10)

        data = zeros(10, 10)

        @test estimator(data) ≈ 0.0
        @test estimator(data, dims = 1) ≈ zeros(1, 10)
        @test estimator(data, dims = 2) ≈ zeros(10)

        data = randn(100, 100)
    end
end

@testitem "background/estimators: SourceExtractorBackground" setup=[estimators] begin
    # test skewed distribution
    data = float(collect(1:100))
    data[71:end] .= 1e7

    @test SourceExtractorBackground()(data) ≈ median(data)
end

@testitem "background/estimators: BiweightLocationBackground" setup=[estimators] begin
    b = BiweightLocationBackground()
    @test b([1, 3, 5, 500, 2]) ≈ 2.745 atol = 1e-3
end

###############################################################################
# RMS Estimators

@testitem "background/estimators: Trivial RMS estimator" setup=[estimators] begin
    @testset "Trivial $E"  for E in [StdRMS, MADStdRMS, BiweightScaleRMS]
        estimator = E()

        @test estimator(ones(1)) == 0

        data = ones(10, 10)
        @test estimator(data) ≈ 0.0
        @test estimator(data, dims = 1) ≈ zeros(1, 10)
        @test estimator(data, dims = 2) ≈ zeros(10)


        data = randn(10000, 10000)
        @test estimator(data) ≈ 1 atol = 3e-2
    end
end

@testitem "background/estimators: StdRMS" setup=[estimators] begin
    s = StdRMS()
    data = randn(100)
    @test s(data) == std(data, corrected = false)
end

@testitem "background/estimators: MADStdRMS" setup=[estimators] begin
    s = MADStdRMS()
    data = randn(100)
    @test s(data) == mad(data, normalize = true)
end

@testitem "background/estimators: BiweightScaleRMS" setup=[estimators] begin
    s = BiweightScaleRMS()
    @test s([1, 3, 5, 500, 2]) ≈ 1.70992562072
end
