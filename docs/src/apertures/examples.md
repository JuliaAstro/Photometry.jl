# Examples

## Plotting
    
We have recipes for all our aperture types, so you can easily create overlays on your images.

```@example plot
using Photometry
using Plots

plot(CircularAperture(2, 3, 4), c=1, xlims=(0, 9), ylims=(0, 9))
plot!(CircularAnnulus(5, 5, 2.1, 3), c=2)
plot!(EllipticalAperture(0, 0, 10, 1, 32), c=3)

savefig("apertures.png"); nothing # hide
```

![](apertures.png)

## Simple Stars

Here is an example where we will find aperture fluxes for stars from M67. The dataset is provided as part of the [astropy/photutils-datasets](https://github.com/astropy/photutils-datasets) repository.

Let's start by downloading and showing our image

```@example stars
using Photometry
using Plots
using FITSIO

hdu = FITS(download("https://github.com/astropy/photutils-datasets/raw/master/data/M6707HH.fits"))
image = read(hdu[1])'
chunk = @view image[71:150, 81:155]

heatmap(chunk, aspect_ratio=1, c=:inferno,
    xlims=(1, size(chunk, 2)), ylims=(1, size(chunk, 1)))
savefig("m67.png"); nothing # hide
```

![](m67.png)

Now let's add some apertures!

```@example stars
positions = [
    [47.5, 67.5],
    [29.5, 62.5],
    [23.5, 48.5],
    [17.5, 29.5],
    [13.5, 10.5],
    [65.5, 14.0]
]

radii = [3, 3, 2.7, 2, 2.7, 3]

aps = CircularAperture.(positions, radii)
```

now let's plot them up

```@example stars
heatmap(chunk, aspect_ratio=1, c=:inferno,
    xlims=(1, size(chunk, 2)), ylims=(1, size(chunk, 1)))
plot!.(aps, c=:white)
savefig("m67_aps.png"); nothing # hide
```

![](m67_aps.png)

and finally let's get our output table for the photometry

```@example stars
table = aperture_photometry(aps, chunk)
```
