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
url = "https://rawcdn.githack.com/astropy/photutils-datasets/8c97b4fa3a6c9e6ea072faeed2d49a20585658ba/data/M6707HH.fits"
hdu = FITS(download(url))
image = read(hdu[1])

# Plot
function imshow(image; kwargs...)
    xs, ys = axes(image)
    data = transpose(image)
    heatmap(xs, ys, data;
            aspect_ratio=1,
            xlim=extrema(xs), ylim=extrema(ys),
            kwargs...)
end

imshow(image)
```

Now let's try and estimate the background using [`estimate_background`](@ref). First, we'll sigma-clip to try and remove the signals from the stars. Then, the background is broken down into boxes, in this case of size `(50, 50)`. Within each box, the given statistical estimators get the background value and RMS. By default, we use [`SourceExtractorBackground`](@ref) and [`StdRMS`](@ref). This creates a low-resolution image, which we then need to resize. We can accomplish this using an interpolator, by default a cubic-spline interpolator via [`ZoomInterpolator`](@ref). The end result is a smooth estimate of the spatially varying background and background RMS.

```@example bkg
# sigma-clip
clipped = sigma_clip(image, 1, fill=NaN)

# get background and background rms with box-size (50, 50)
bkg, bkg_rms = estimate_background(clipped, 50)

# plot
plot(
    imshow(image, title="Original"),
    imshow(clipped, title="Sigma-Clipped"),
    imshow(bkg, title="Background"),
    imshow(bkg_rms, title="Background RMS"),
    layout=(2, 2), ticks=false,
)
```

We could apply a median filter, too, by specifying `filter_size`

```@example bkg
# get background and background rms with box-size (50, 50) and filter_size (5, 5)
bkg_f, bkg_rms_f = estimate_background(clipped, 50, filter_size=5)

# plot
plot(
    imshow(bkg, title="Unfiltered", ylabel="Background"),
    imshow(bkg_f, title="Filtered"),
    imshow(bkg_rms, ylabel="RMS"),
    imshow(bkg_rms_f);
    layout=(2, 2), ticks=false,
)
```

Now we can see our image after subtracting the filtered background and ready for [Aperture Photometry](@ref)!

```@example bkg
subt = image .- bkg_f[axes(image)...]
plot(
    imshow(image, title="Original", colorbar=false),
    imshow(subt, title="Subtracted");
    layout=(1, 2), size=(600, 260),
    xlims=(400, 800), ylims=(400, 800),
    clims=(minimum(subt), maximum(image)),
    ticks=false, aspect_ratio=1,
)
```

### IDW Interpolator

Here is a quick example using the [`IDWInterpolator`](@ref)

```@example bkg
b1, r1 = estimate_background(clipped, 50, filter_size=5)
b2, r2 = estimate_background(clipped, 50, itp=IDWInterpolator(50), filter_size=5)

plot(
    imshow(b1, title="ZoomInterpolator", ylabel="Background"),
    imshow(b2, title="IDWInterpolator"),
    imshow(r1, ylabel="RMS"),
    imshow(r2);
    layout=(2, 2), ticks=false,
)
```

## API/Reference

```@docs
estimate_background
sigma_clip
sigma_clip!
Background.validate_SE
```
