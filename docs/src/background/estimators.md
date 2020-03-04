```@meta
DocTestSetup = quote
    using Random
    Random.seed!(123456)
end
```

# Background Estimators

All of these estimators are subtypes of [`Background.BackgroundEstimator`](@ref) and are derived using various statistical and image processing methods.

## Simple Statistics

```@docs
Mean
Median
Mode
```

## API/Reference

```@docs
Background.BackgroundEstimator
```
