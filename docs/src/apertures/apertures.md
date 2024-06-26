```@meta
DocTestSetup = :(using Photometry)
```

# Apertures

All apertures will rely on a position and the shape parameters.
```julia
aperture = Aperture(x0, y0, shape_params...)
```
The position can be pixels or sky coordinates. The sky coordinate positions utilize [SkyCoords.jl](https://juliaastro.github.io/SkyCoords.jl/stable) and [WCS.jl](https://juliaastro.github.io/WCS.jl/stable) for conversion.

!!! warning
    Sky coordinates are not supported yet.

!!! note
    See [Pixel Convention](@ref) - The origin is the bottom-left with `(1, 1)` being the center of the pixel.


## API/Reference

```@docs
Aperture.AbstractAperture
Subpixel
Aperture.area_arc
Aperture.circular_overlap_core
Aperture.circular_overlap_single_exact
Aperture.inside_ellipse
Base.size(::Aperture.AbstractAperture)
Aperture.area_triangle
Aperture.inside_rectangle
Aperture.bounds
```

### Circular Apertures

These apertures are parameterized by radius.

```@docs
CircularAperture
CircularAnnulus
```

### Elliptical Apertures

These apertures are parameterized by the semi-major axis `a`, semi-minor axis `b` and position angle in degrees counter-clockwise from the positive x-axis `θ`


```@docs
EllipticalAperture
EllipticalAnnulus
```


### Rectangular Apertures

These apertures are parameterized by width `w`, height `h`, and position angle in degrees counter-clockwise from the positive x-axis `θ`.

```@docs
RectangularAperture
RectangularAnnulus
```
