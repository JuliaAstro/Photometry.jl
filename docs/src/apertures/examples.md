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
```

## Simple Stars

Here is an example where we will find aperture fluxes for stars from M67. The dataset is provided as part of the [astropy/photutils-datasets](https://github.com/astropy/photutils-datasets) repository.

Let's start by downloading and showing our image

```@example stars
using Photometry
using Plots
using FITSIO

# Load data in
url = "https://rawcdn.githack.com/astropy/photutils-datasets/8c97b4fa3a6c9e6ea072faeed2d49a20585658ba/data/M6707HH.fits"
hdu = FITS(download(url))
chunk = read(hdu[1], 81:155, 71:150)

# Plot
function imshow(image; kwargs...)
    xs, ys = axes(image)
    data = transpose(image)
    heatmap(xs, ys, data; aspect_ratio=1, xlim=extrema(xs), ylim=extrema(ys), kwargs...)
end

imshow(chunk)
```

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
imshow(chunk)
plot!(aps, c=:white)
```

and finally let's get our output table for the photometry

```@example stars
table = photometry(aps, chunk)
```

## Stars with Spatial Background Subtraction

This example will be the same as [Simple Stars](@ref) but will add background estimation using the tools in [Background Estimation](@ref)

```@example stars
clipped = sigma_clip(chunk, 1, fill=NaN)
# Estimate 2D spatial background using boxes of size (5, 5)
bkg, bkg_rms = estimate_background(clipped, 5)

plot(
    imshow(chunk, title="Original"),
    imshow(clipped, title="Sigma-Clipped"),
    imshow(bkg, title="Background"),
    imshow(bkg_rms, title="Background RMS");
    layout=(2, 2), size=(600, 600), ticks=false
)
```

Now, using the same apertures, let's find the output using the background-subtracted image

```@example stars
plot(
    imshow(chunk, title="Original"),
    imshow(chunk .- bkg, title="Subtracted");
    layout=2, size=(600, 260), ticks=false, colorbar=false
)
plot!(aps, c=:white, subplot=1)
plot!(aps, c=:white, subplot=2)
```

```@example stars
table = photometry(aps, chunk .- bkg, bkg_rms)
```
