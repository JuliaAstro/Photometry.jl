# Apertures

All apertures will rely on a position and the shape parameters.
```julia
aperture = Aperture(x0, y0, shape_params...)
```
The position can be pixels or sky coordinates. The sky coordinate positions utilize [SkyCoords.jl](https://juliaastro.github.io/SkyCoords.jl/stable) and [WCS.jl](https://juliaastro.github.io/WCS.jl/stable) for conversion. 

## Circular Apertures

These apertures are parametrized by radius.

```@docs
CircularAperture
CircularAnnulus
```

## Elliptical Apertures

These apertures are parametrized by the semi-major axis `a` and semi-minor axis `b`.
