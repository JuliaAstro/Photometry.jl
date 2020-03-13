# Background Estimation

The module provides tools and algorithms for estimating the background of astronomical data.

## Usage

## Interpolators

Background interpolators provide a method for converting low-resolution meshes into low-order high-resolution images.

```@docs
Background.BackgroundInterpolator
ZoomInterpolator
```

## API/Reference

```@docs
estimate_background
sigma_clip
sigma_clip!
```
