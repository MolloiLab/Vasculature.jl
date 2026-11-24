# Why the package's X table misses X = A − kY at small k: its recursion in Float32 against the same in Float64, and with the
# slope weight ∫(1 − s/dt) e^{−k(dt−s)} ds = a − b. One made-up input, no XCAT.
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
import Vasculature as V
using Printf
const F = V.Flow
inp = F.gamma_input(; t0 = 12.0f0, tpeak = 19.0f0, peak_blood = 15.0f0, alpha = 3.0f0, hd = 0.45f0, dt = 0.02f0, tmax = 120.0f0)
c = Float64.(inp.c); dt = Float64(inp.dt); nt = length(c)
A = zeros(nt); for i in 2:nt; A[i] = A[i - 1] + 0.5dt * (c[i - 1] + c[i]); end
function table(k, T; slope_weight = :package)
    k = T(k); h = T(dt); e = exp(-k * h); a = (one(T) - e) / k; b = (h - a) / (k * h)
    w = slope_weight === :package ? b : a - b
    X = zeros(T, nt); Y = zeros(T, nt)
    for i in 2:nt
        X[i] = e * X[i - 1] + T(c[i]) * a - (T(c[i]) - T(c[i - 1])) * w
        Y[i] = Y[i - 1] + T(0.5) * h * (X[i - 1] + X[i])
    end
    return Float64.(X), Float64.(Y), Float64(b), Float64(a - b)
end
for k in (1.0e-4, 1.52e-3, 7.17e-3, 0.1, 1.0, 10.0)
    gap(X, Y) = maximum(abs.(X .- A .+ k .* Y)) / maximum(X)
    X32, Y32, b32, _ = table(k, Float32)
    X64, Y64, b64, ab64 = table(k, Float64)
    Xe, Ye, _, _ = table(k, Float64; slope_weight = :exact)
    @printf("k %.2e (k·dt %.1e): b in Float32 %.4g s, in Float64 %.4g s (a − b %.4g); |X − A + kY| / max X: package Float32 %.1e, Float64 %.1e, Float64 with a − b %.1e\n",
        k, k * dt, b32, b64, ab64, gap(X32, Y32), gap(X64, Y64), gap(Xe, Ye))
end
