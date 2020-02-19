# Examples

## Plotting
    
We have recipes for all our aperture types, so you can easily create overlays on your images.

```@example plot
using Photometry
using Plots

plot(CircularAperture(2, 3, 4), c=1, xlims=(0, 9), ylims=(0, 9))
plot!(CircularAnnulus(5, 5, 2.1, 3), c=2)
plot!(EllipticalAperture(0, 0, 10, 1, 32), c=3)

savefig("apertures.svg"); nothing # hide
```

![](apertures.svg)

## Simple Stars

Here is an example where we will find aperture fluxes for stars from M67. The dataset is provided as part of the [astropy/photutils-datasets](https://github.com/astropy/photutils-datasets) repository.

Let's start by downloading and showing our image

```@example stars
using Photometry
using Plots
using FITSIO

hdu = FITS(download("https://github.com/astropy/photutils-datasets/raw/master/data/M6707HH.fits"))
image = read(hdu[1])'
cutout = @view image[71:150, 81:155]

heatmap(cutout, aspect_ratio=1, xlims=(1, size(cutout, 2)), ylims=(1, size(cutout, 1)), c=:inferno)
savefig("m67.svg"); nothing # hide
```

![](m67.svg)

Now let's add some apertures!

```@example stars
positions = [
    [48, 68],
    [30, 63],
    [24, 49],
    [18, 30],
    [14, 11],
    [66, 14.5]
]

radii = [3, 3, 2.7, 2, 2.7, 3]

aps = CircularAperture.(positions, radii)
```

now let's plot them up

```@example stars
heatmap(cutout, aspect_ratio=1, xlims=(1, size(cutout, 2)), ylims=(1, size(cutout, 1)), c=:inferno)
plot!.(aps, c=:white)
savefig("m67_aps.svg"); nothing # hide
```

![](m67_aps.svg)

