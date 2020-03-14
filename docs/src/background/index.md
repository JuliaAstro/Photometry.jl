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
gma-clip to try and remove the signals from the stars. Then, the background is broken down into meshes, in this case of size `(50, 50)`. Within each mesh, the given statistical estimators get the background value and RMS. By default, we use [`SourceExtractor`](@ref) and [`StdRMS`](@ref). This creates a low-resolution image, which we then need to resize. We can accomplish this using an interpolator, by default a cubic-spline interpolator via [`ZoomInterpolator`](@ref). The end result is a smooth estimate of the spatially varying background and background RMS.

```@example bkg
# sigma-clip
clipped = sigma_clip(image, 1, fill=NaN)

# get background and background rms with mesh-size (50, 50)
bkg, bkg_rms = estimate_background(clipped, 50)

# plot
plot(layout=(2, 2), size=(800, 800), link=:all)
heatmap!(image,   title="Original",       subplot=1)
heatmap!(clipped, title="Sigma-Clipped",  subplot=2)
heatmap!(bkg,     title="Background",     subplot=3)
heatmap!(bkg_rms, title="Background RMS", subplot=4)

savefig("bkg.png"); nothing # hide
```

![](bkg.png)

and now we can see our image with the background subtracted and ready for [Aperture Photometry](@ref)!

```@example bkg
subt = image .- bkg[:1059, :1059]
plot(layout=(1, 2), link=:all, size=(800, 350), xlims=(400, 800), ylims=(400, 800))
heatmap!(image, title="Original",   subplot=1)
heatmap!(subt,  title="Subtracted", subplot=2)

savefig("bkg_final.png"); nothing # hide
```

![](bkg_final.png)


## API/Reference

```@docs
estimate_background
sigma_clip
sigma_clip!
```
