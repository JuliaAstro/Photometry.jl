# Background Estimation

The module provides tools and algorithms for estimating the background of astronomical data.

## Usage

Estimating backgrounds is an important step in performing photometry. Ideally, we could perfectly describe the background with a scalar value or with some distribution. Unfortunately, it's impossible for us to precisely separate the background and foreground signals. Here, we use mixture of robust statistical estimators and meshing to let us get the spatially varying background from an astronomical photo.

Let's show an example

```@example bkg
using Photometry
using FITSIO
using Plots

# Download our image, courtesy of astropy
hdu = FITS(download("https://github.com/astropy/photutils-datasets/raw/master/data/M6707HH.fits"))
image = read(hdu[1])'

default(aspect_ratio=1, xlims=(1, size(image, 2)), ylims=(1, size(image, 1)))

heatmap(image)
savefig("m67_full.png"); nothing # hide
```

![](m67_full.png)

Now let's try and estimate the background using [`estimate_background`](@ref). First, we'll si
gma-clip to try and remove the signals from the stars. Then, the background is broken down into meshes, in this case of size `(50, 50)`. Within each mesh, the given statistical estimators get the background value and RMS. By default, we use [`SourceExtractorBackground`](@ref) and [`StdRMS`](@ref). This creates a low-resolution image, which we then need to resize. We can accomplish this using an interpolator, by default a cubic-spline interpolator via [`ZoomInterpolator`](@ref). The end result is a smooth estimate of the spatially varying background and background RMS.

```@example bkg
# sigma-clip
clipped = sigma_clip(image, 1, fill=NaN)

# get background and background rms with mesh-size (50, 50)
bkg, bkg_rms = estimate_background(clipped, 50)

# plot
plot(layout=(2, 2), size=(800, 800), link=:all, ticks=false)
heatmap!(image, title="Original", subplot=1)
heatmap!(clipped, title="Sigma-Clipped", subplot=2)
heatmap!(bkg, title="Background", subplot=3)
heatmap!(bkg_rms, title="Background RMS", subplot=4)

savefig("bkg.png"); nothing # hide
```

![](bkg.png)

We could apply a median filter, too, by specifying `filter_size`

```@example bkg
# get background and background rms with mesh-size (50, 50) and filter_size (5, 5)
bkg_f, bkg_rms_f = estimate_background(clipped, 50, filter_size=5)

# plot
plot(layout=(2, 2), size=(800, 800), link=:all, ticks=false)
heatmap!(bkg, title="Unfiltered", ylabel="Background", subplot=1)
heatmap!(bkg_f, title="Filtered", subplot=2)
heatmap!(bkg_rms, ylabel="RMS", subplot=3)
heatmap!(bkg_rms_f, subplot=4)

savefig("bkg_filtered.png"); nothing # hide
```

![](bkg_filtered.png)

Now we can see our image after subtracting the filtered background and ready for [Aperture Photometry](@ref)!

```@example bkg
subt = image .- bkg_f[axes(image)...]
plot(layout=(1, 2),
    link=:all,
    size=(800, 350),
    xlims=(400, 800),
    ylims=(400, 800),
    clims=(minimum(subt), maximum(image)),
    ticks=false)
heatmap!(image, title="Original", colorbar=false, subplot=1)
heatmap!(subt, title="Subtracted", subplot=2)

savefig("bkg_final.png"); nothing # hide
```

![](bkg_final.png)

### IDW Interpolator

Here is a quick example using the [`IDWInterpolator`](@ref)

```@example bkg
b1, r1 = estimate_background(clipped, 50, filter_size=5)
b2, r2 = estimate_background(clipped, 50, SourceExtractorBackground(), StdRMS(), IDWInterpolator(50), filter_size=5)

plot(layout=(2, 2), size=(800, 800), link=:all, ticks=false)
heatmap!(b1, title="ZoomInterpolator", ylabel="Background", subplot=1)
heatmap!(b2, title="IDWInterpolator", subplot=2)
heatmap!(r1, ylabel="RMS", subplot=3)
heatmap!(r2, subplot=4)

savefig("bkg_idw_comp.png"); nothing # hide
```

![](bkg_idw_comp.png)

## API/Reference

```@docs
estimate_background
sigma_clip
sigma_clip!
```
