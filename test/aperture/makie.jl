using Photometry.Aperture:
    Subpixel,
    CircularAperture,
    CircularAnnulus,
    EllipticalAperture,
    EllipticalAnnulus,
    RectangularAperture,
    RectangularAnnulus

using Makie: Makie, Point2d

const MAKIE_APERTURES = [
    CircularAperture(3, 3, 3),
    CircularAnnulus(3, 3, 2, 4),
    EllipticalAperture(3, 3, 4, 2, 45),
    EllipticalAnnulus(3, 3, 3, 4, 2, -26),
    RectangularAperture(3, 3, 3, 4, 15),
    RectangularAnnulus(3, 3, 3, 4, 2, -5),
]

isannulus(ap) = ap isa Union{CircularAnnulus, EllipticalAnnulus, RectangularAnnulus}
nanpoint(p) = any(isnan, p)

# Split a NaN-separated point list into per-ring index ranges
function split_at_nan(points)
    ranges = UnitRange{Int}[]
    start = 1
    for (i, p) in enumerate(points)
        if nanpoint(p)
            push!(ranges, start:(i - 1))
            start = i + 1
        end
    end
    push!(ranges, start:length(points))
    return ranges
end

@testset "aperture/makie: Aperture Makie conversions" begin
    @testset "outline - $(typeof(ap))" for ap in MAKIE_APERTURES
        @test Makie.plottype(ap) == Makie.Lines
        (points,) = Makie.convert_arguments(Makie.PointBased(), ap)
        @test points isa Vector{Point2d}

        rings = [points[r] for r in split_at_nan(points) if !isempty(r)]

        # Unlike the Plots recipes (which shift by +0.5 for Plots' cell-edge
        # convention), outlines center on the aperture's own (x, y).
        finite = filter(!nanpoint, points)
        N = length(finite)
        x_mean = sum(p[1] for p in finite) / N
        y_mean = sum(p[2] for p in finite) / N
        @test x_mean ≈ 3 atol = 3 / sqrt(N)
        @test y_mean ≈ 3 atol = 3 / sqrt(N)

        # Annuli carry both rings separated by a NaN break; every ring closes.
        @test count(nanpoint, points) == (isannulus(ap) ? 1 : 0)
        for ring in rings
            @test first(ring) ≈ last(ring) atol = 1.0e-8
        end
    end

    @testset "npoints argument" begin
        ap = CircularAperture(3, 3, 3)
        (points,) = Makie.convert_arguments(Makie.PointBased(), ap, 25)
        @test length(points) == 25
    end

    @testset "vector of apertures" begin
        @test Makie.plottype(MAKIE_APERTURES) == Makie.Lines
        (points,) = Makie.convert_arguments(Makie.PointBased(), MAKIE_APERTURES)
        @test points isa Vector{Point2d}
        # one break between consecutive apertures plus one inside each annulus
        nannuli = count(isannulus, MAKIE_APERTURES)
        @test count(nanpoint, points) == length(MAKIE_APERTURES) - 1 + nannuli
        @test Makie.convert_arguments(Makie.PointBased(), CircularAperture{Float64}[]) == (Point2d[],)

        # npoints threads through the vector method too
        (points,) = Makie.convert_arguments(Makie.PointBased(), MAKIE_APERTURES, 25)
        nrings = length(MAKIE_APERTURES) + nannuli
        @test count(!nanpoint, points) == 25 * nrings
    end

    @testset "poly conversion - $(typeof(ap))" for ap in MAKIE_APERTURES
        (polygon,) = Makie.convert_arguments(Makie.Poly, ap)
        @test polygon isa Makie.Polygon
        # annuli carry the inner ring as an interior, cutting a hole in the fill
        @test length(polygon.interiors) == (isannulus(ap) ? 1 : 0)

        (polygon,) = Makie.convert_arguments(Makie.Poly, ap, 25)
        @test length(polygon.exterior) == 25
    end

    @testset "poly conversion - vector of apertures" begin
        (polygons,) = Makie.convert_arguments(Makie.Poly, MAKIE_APERTURES)
        @test polygons isa Vector{<:Makie.Polygon}
        @test length(polygons) == length(MAKIE_APERTURES)

        (polygons,) = Makie.convert_arguments(Makie.Poly, MAKIE_APERTURES, 25)
        @test all(p -> length(p.exterior) == 25, polygons)
    end

    @testset "Subpixel unwraps" begin
        ap = CircularAperture(3, 3, 3)
        @test Makie.convert_arguments(Makie.PointBased(), Subpixel(ap)) ==
            Makie.convert_arguments(Makie.PointBased(), ap)
        @test Makie.convert_arguments(Makie.Poly, Subpixel(ap)) ==
            Makie.convert_arguments(Makie.Poly, ap)
    end
end
