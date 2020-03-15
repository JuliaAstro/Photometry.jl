# Background Estimators

All of these estimators are subtypes of [`Background.BackgroundEstimator`](@ref) or [`Background.BackgroundRMSEstimator`](@ref) and are derived using various statistical and image processing methods.

## Location Estimators

These estimators are used for estimating the background using some form of a central statistic.

```@docs
Background.BackgroundEstimator
MeanBackground
MedianBackground
ModeBackground
MMMBackground
SourceExtractorBackground
BiweightLocationBackground
```

## RMS Estimators

These estimators are used for estimating the root-mean-square (RMS) of the background using some form of a deviation statistic.

```@docs
Background.BackgroundRMSEstimator
StdRMS
MADStdRMS
BiweightScaleRMS
```
