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
```

### Circular Apertures

These apertures are parametrized by radius.

```@docs
CircularAperture
CircularAnnulus
```

### Elliptical Apertures

These apertures are parametrized by the semi-major axis `a`, semi-minor axis `b` and position angle in degrees counter-clockwise from the positive x-axis `θ`


```@docs
EllipticalAperture
EllipticalAnnulus
```


### Rectangular Apertures

These apertures are parametrized by width `w`, height `h`, and position angle in degrees counter-clockwise from the positive x-axis `θ`.

```@docs
RectangularAperture
RectangularAnnulus
```
