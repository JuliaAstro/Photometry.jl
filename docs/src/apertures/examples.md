# Examples

## Plotting
We have recipes for all our aperture types, so you can easily create overlays on your images.

```@example plot
using Photometry
using Plots

plot(CircularAperture(2, 3, 4), c=1, xlims=(-1, 12), ylims=(0, 9))
plot!(CircularAnnulus(5, 5, 2.1, 3), c=2)
plot!(EllipticalAperture(0, 0, 10, 1, 32), c=3)
plot!(EllipticalAnnulus(5, 5, 4, 5, 2, -32), c=4)
plot!(RectangularAperture(0, 0, 4, 4, 4), c=5)
plot!(RectangularAnnulus(5, 1, 3, 4, 4, 4), c=6)
savefig("plot-1.png"); nothing #hide
```

![](plot-1.png)

## Simple Stars

Here is an example where we will find aperture fluxes for stars from M67. The dataset is provided as part of the [astropy/photutils-datasets](https://github.com/astropy/photutils-datasets) repository.

Let's start by downloading and showing our image

```@example stars
using Photometry
using Plots
using FITSIO

# Load data in
hdu = FITS(download("https://github.com/astropy/photutils-datasets/raw/master/data/M6707HH.fits"))
image = read(hdu[1])'
chunk = image[71:150, 81:155]

# Plot
default(aspect_ratio=1, xlims=(1, size(chunk, 2)), ylims=(1, size(chunk, 1)))

heatmap(chunk)
savefig("plot-2.png"); nothing # hide
```

![](plot-2.png)

Now let's add some apertures!

```@example stars
positions = [
    [47.5 , 67.5],
    [29.5 , 62.5],
    [23.5 , 48.5],
    [17.5 , 29.5],
    [13.25, 10.5],
    [65.5 , 14.0]
]

radii = [3, 3, 2.7, 2, 2.7, 3]

aps = CircularAperture.(positions, radii)
```

now let's plot them up

```@example stars
p = heatmap(chunk)
plot!.(aps, c=:white)
p
savefig("plot-3.png"); nothing # hide
```

![](plot-3.png)

and finally let's get our output table for the photometry

```@example stars
table = aperture_photometry(aps, chunk)
```

## Stars with Spatial Background Subtraction

This example will be the same as [Simple Stars](@ref) but will add background estimation using the tools in [Background Estimation](@ref)

```@example stars
clipped = sigma_clip(chunk, 1, fill=NaN)
# Estimate 2D spatial background using boxes of size (5, 5)
bkg, bkg_rms = estimate_background(clipped, 5)

plot(layout=(2, 2), size=(600, 600), ticks=false)
heatmap!(chunk, title="Original", subplot=1)
heatmap!(clipped, title="Sigma-Clipped", subplot=2)
heatmap!(bkg, title="Background", subplot=3)
heatmap!(bkg_rms, title="Background RMS", subplot=4)
savefig("plot-4.png"); nothing # hide
```

![](plot-4.png)

Now, using the same apertures, let's find the output using the background-subtracted image

```@example stars
p = plot(layout=(1, 2),
    clims=(minimum(chunk .- bkg),
    maximum(chunk)),
    size=(600, 260),
    ticks=false)
heatmap!(chunk, title="Original", colorbar=false, subplot=1)
heatmap!(chunk .- bkg, title="Subtracted", subplot=2)
plot!.(aps, c=:white, subplot=1)
plot!.(aps, c=:white, subplot=2)
p
savefig("plot-5.png"); nothing # hide
```

![](plot-5.png)

```@example stars
table = aperture_photometry(aps, chunk .- bkg, bkg_rms)
```
