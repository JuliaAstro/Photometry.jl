# Examples

TODO

## Plotting

!!! warning
    This is still a work-in-progress

    
We have recipes for all our aperture types, so you can easily create overlays on your images.

```@example plot
using Photometry
using Plots

plot(CircularAperture(2, 3, 4), c=1, xlims=(0, 9), ylims=(0, 9))
plot!(CircularAnnulus(5, 5, 2.1, 3), c=2)
plot!(EllipticalAperture(0, 0, 10, 1, -32), c=3)

savefig("apertures.svg"); nothing # hide
```

![](apertures.svg)
