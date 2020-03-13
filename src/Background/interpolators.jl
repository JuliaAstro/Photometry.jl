using Interpolations: CubicSplineInterpolation

"""
    ZoomInterpolator(factors)

Use a cubic-spline interpolation scheme to increase resolution of a mesh.

`factors` represents the level of "zoom", so an input mesh of size `(10, 10)` with factors `(2, 2)` will have an output size of `(20, 20)`. If only an integer is provided, it will be used as the factor for every axis. 

# Examples
```jldoctest
julia> ZoomInterpolator(2)([1 0; 0 1])
4×4 Array{Float64,2}:
  1.0          0.666667  0.333333  -2.77556e-17
  0.666667     0.555556  0.444444   0.333333   
  0.333333     0.444444  0.555556   0.666667   
 -5.55112e-17  0.333333  0.666667   1.0

julia> ZoomInterpolator(2, 3)([1 0; 0 1])
4×6 Array{Float64,2}:
  1.0          0.8  0.6       0.4       0.2  -2.77556e-17
  0.666667     0.6  0.533333  0.466667  0.4   0.333333   
  0.333333     0.4  0.466667  0.533333  0.6   0.666667   
 -5.55112e-17  0.2  0.4       0.6       0.8   1.0

```
"""
struct ZoomInterpolator <: BackgroundInterpolator
factors::NTuple{2,<:Integer}
end

ZoomInterpolator(factor::Integer) = ZoomInterpolator((factor, factor))
ZoomInterpolator(f1::Integer, args...) = ZoomInterpolator((f1, args...))

function (z::ZoomInterpolator)(mesh)
out_axes = [range(1, n, length = n * f) for (n, f) in zip(size(mesh), z.factors)]
itp = CubicSplineInterpolation(axes(mesh), mesh)
return [itp(i, j) for i in out_axes[1], j in out_axes[2]]
end
