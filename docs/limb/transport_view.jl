# The leg's contrast transport at rest as one page (transport_view_page.jl): the main vessels (pass 1) in front view,
# coloured by their blood iodine every half second, every patch's tissue as a dot, and the iodine in the leg's arteries,
# tissue and veins over time, the two in step.
#
#     julia -t auto transport_view.jl          from transport_rest (run transport_leg.jl first): the data, then the page
#     julia -t auto transport_view.jl page     the page again from the saved data (view_data.jls)
#
# Every bed's iodine comes from its saved arrival times (beds_<k>.jls) and the final arterial input, in double precision:
#   - a patch's arteries hold its inflow less what has reached its beds;
#   - a bed holds what reached it less what left into its venule, the exchange through X_k(t) = ∫ c e^{−k(t−s)} ds;
#   - a patch's veins hold what left its beds less what reached the venous patches' roots.
# The main arteries are plug flow from the leg's inlet, so each carries the input delayed by its own arrival time (the
# page applies the delay). Each main vein carries the outflows of the venous patches it drains, each delayed by its path.
#
# Numbers only. The page is made from XCAT: it stays in LIMB_CACHE_DIR\transport_rest (lab storage) and is never published.

using Printf, Serialization, Statistics, Base64

include(joinpath(@__DIR__, "transport_leg.jl"))
include(joinpath(@__DIR__, "transport_view_page.jl"))

const VIEW_DT = 0.5
const VT = collect(0.0:VIEW_DT:Float64(TMAX))
const NVF = length(VT)
const KG = collect(range(log(1.0e-6), log(10.0), length = 200)) # back-flux rates (1/s), wider than InputCurve's: no bed clamped
const DKG = KG[2] - KG[1]
const VIEW_DATA = joinpath(T_DIR, "view_data.jls")
const VIEW_PAGE = joinpath(T_DIR, "transport_view.html")

"X_k(t) = ∫ c(s) e^{−k(t−s)} ds on the grid for every k of KG, in double precision, exact for piecewise-linear c; (k, time), so neighbouring rates sit together."
function xtable(c)
    X = zeros(length(KG), NG)
    Threads.@threads for m in eachindex(KG)
        k = exp(KG[m])
        x = k * FINE
        e = exp(-x)
        a = -expm1(-x) / k
        d = x < 1.0e-3 ? x^2 / 2 - x^3 / 6 + x^4 / 24 : x + expm1(-x)
        b = d / (k * x) # the later sample's weight; the earlier one's is a − b
        v = 0.0
        for i in 2:NG
            v = e * v + b * c[i] + (a - b) * c[i - 1]
            X[m, i] = v
        end
    end
    return X
end

@inline function lerp(v, u)
    u <= 0.0 && return 0.0
    x = u / FINE
    i = unsafe_trunc(Int, x)
    i >= NG - 1 && return @inbounds v[NG]
    f = x - i
    return @inbounds v[i + 1] + f * (v[i + 2] - v[i + 1])
end

@inline function xlook(X, m, h, u)
    u <= 0.0 && return 0.0
    x = u / FINE
    i = unsafe_trunc(Int, x)
    i >= NG - 1 && return @inbounds (1 - h) * X[m, NG] + h * X[m + 1, NG]
    f = x - i
    @inbounds lo = X[m, i + 1] + f * (X[m, i + 2] - X[m, i + 1])
    @inbounds hi = X[m + 1, i + 1] + f * (X[m + 1, i + 2] - X[m + 1, i + 1])
    return (1 - h) * lo + h * hi
end

"Every bed's iodine (mg) every VIEW_DT, added to its patch's arteries, beds and veins and to its venous patch's outflow (a function of its own, so the loop compiles for the arrays' types)."
function bed_loop!(art, bed, ven, out, fpm, fps, exs, kxs, tins, tas, tvs, tcs, svs, cols, A, X)
    for i in eachindex(fps)
        fp = Float64(fps[i])
        E = Float64(exs[i])
        y = clamp((log(Float64(kxs[i])) - KG[1]) / DKG, 0.0, length(KG) - 1 - 1.0e-9)
        m = unsafe_trunc(Int, y)
        h = y - m
        m += 1
        tin = Float64(tins[i])
        t4 = tin + tas[i] + tvs[i] + tcs[i]
        t5 = t4 + svs[i]
        c = cols[i]
        fpm[c] += fp
        @inbounds for f in searchsortedfirst(VT, tin):NVF
            T = VT[f]
            ain = fp * lerp(A, T - tin)
            g4 = fp * (lerp(A, T - t4) - E * xlook(X, m, h, T - t4))
            g5 = fp * (lerp(A, T - t5) - E * xlook(X, m, h, T - t5))
            art[f] -= ain
            bed[f] += ain - g4
            ven[f] += g4 - g5
            out[f, c] += g5
        end
    end
    return nothing
end

"Patch `k`'s iodine (mg) every VIEW_DT in its arteries, beds and veins, and what it has sent to each venous patch's root, with that root's plasma flow."
function patch_iodine(k, W, A, X)
    b = deserialize(joinpath(T_DIR, "beds_$k.jls"))
    vleaf = deserialize(joinpath(FLOW_DIR, "vleaf_$k.jls"))
    vf = forest(k, W)
    mi = vf.segm[vleaf]
    length(mi) == length(b.fp) || error("patch $k: $(length(mi)) venules for $(length(b.fp)) beds")
    ms = sort!(unique(mi))
    col = Dict(m => i for (i, m) in enumerate(ms))
    art, bed, ven = zeros(NVF), zeros(NVF), zeros(NVF)
    out, fpm = zeros(NVF, length(ms)), zeros(length(ms))
    bed_loop!(art, bed, ven, out, fpm, b.fp, b.ex, b.kx, b.tin, b.ta, b.tv, b.tc, b.sv, Int[col[m] for m in mi], A, X)
    fpk = sum(fpm)
    for f in 1:NVF
        art[f] += fpk * lerp(A, VT[f] - b.D)
    end
    return (; k, region = b.region, D = b.D, fp = fpk, art, bed, ven, ms, out, fpm)
end

"Linear interpolation of samples `y` on VT at time `s`, 0 before the start."
@inline vlerp(y, s) = s <= 0 ? 0.0 : s >= VT[end] ? y[end] :
    (x = s / VIEW_DT; i = unsafe_trunc(Int, x); y[i + 1] + (x - i) * (y[i + 2] - y[i + 1]))

"The front-view silhouette of the leg's soft tissue: rows of runs on cells `2h` mm."
function silhouette()
    soft, g = LEG.soft, LEG.grid
    nx, nz = size(soft, 1), size(soft, 3)
    proj = falses(cld(nx, 2), cld(nz, 2))
    for kz in 1:nz, jy in axes(soft, 2), ix in 1:nx
        soft[ix, jy, kz] && (proj[cld(ix, 2), cld(kz, 2)] = true)
    end
    hs = 2 * g.h
    edge(i, d) = g.origin[d] + (2i - 2.5) * g.h # the lower edge of cell i along axis d
    rows = Vector{Float64}[]
    for K in axes(proj, 2)
        row = Float64[edge(K, 3)]
        I = 1
        while I <= size(proj, 1)
            if proj[I, K]
                J = I
                while J < size(proj, 1) && proj[J + 1, K]
                    J += 1
                end
                push!(row, edge(I, 1), edge(J, 1))
                I = J + 1
            else
                I += 1
            end
        end
        length(row) > 1 && push!(rows, row)
    end
    return (; h = hs, rows)
end

function view_data()
    isfile(VIEW_DATA) && return deserialize(VIEW_DATA)
    t0 = time()
    R = deserialize(joinpath(T_DIR, "result.jls"))
    W = venule_map()
    ks = [s.k for s in R.patches]
    c = R.c
    A = cumin(c)
    X = xtable(c)
    tsay(@sprintf("view: input tables in %.0f s", time() - t0))
    P = Vector{Any}(undef, length(PATCHES.voxels))
    failed = run_each(ks, "view") do k
        P[k] = patch_iodine(k, W, A, X)
    end
    isempty(failed) || error("view failed for patches $(failed)")
    # every venous patch's outflow (mg, cumulative) and plasma flow
    M = length(VEIN_PATCHES.voxels)
    G5, Fm = zeros(NVF, M), zeros(M)
    for k in ks, (i, m) in enumerate(P[k].ms)
        G5[:, m] .+= P[k].out[:, i]
        Fm[m] += P[k].fpm[i]
    end
    q = zeros(length(PATCHES.voxels))
    for k in ks
        q[k] = deserialize(joinpath(FLOW_DIR, "final_$k.jls")).rest.q
    end
    PA = pass1_arteries(q)
    PV = pass1_veins(Fm)
    fp = sum(Fm)
    # the leg's compartments
    G1 = [fp * lerp(A, T) for T in VT]
    G2 = [sum(P[k].fp * lerp(A, T - P[k].D) for k in ks) for T in VT]
    G6 = [sum(vlerp(view(G5, :, m), T - PV.DV[m]) for m in 1:M if Fm[m] > 0) for T in VT]
    art = (G1 .- G2) .+ sum(P[k].art for k in ks)
    tis = sum(P[k].bed for k in ks)
    ven = sum(P[k].ven for k in ks) .+ (vec(sum(G5; dims = 2)) .- G6)
    total = G1 .- G6
    for tq in (30, 60, 120)
        i, j = findfirst(==(tq), VT), findfirst(f -> abs(f - tq) < 1.0e-3, R.frames)
        tsay(@sprintf(
            "view at %d s: total %.1f mg (transport %.1f), tissue %.1f (%.1f), arteries %.1f (%.1f), veins %.1f (%.1f)",
            tq, total[i], R.total[j], tis[i], sum(r -> R.beds[r][j], 1:4), art[i], R.art1[j] + sum(r -> R.art2[r][j], 1:4),
            ven[i], R.ven1[j] + sum(r -> R.ven2[r][j], 1:4),
        ))
    end
    # the main veins: each venous patch's outflow concentration (plasma), carried down its path, flow-weighted
    cm = zeros(NVF, M)
    for m in 1:M
        Fm[m] > 0 || continue
        for f in 1:NVF
            lo, hi = max(f - 1, 1), min(f + 1, NVF)
            cm[f, m] = (G5[hi, m] - G5[lo, m]) / ((hi - lo) * VIEW_DT) / Fm[m]
        end
    end
    tv = pass1.veins.tree
    nv = length(tv.parent)
    num, den = zeros(NVF, nv), zeros(nv)
    for (m, j0) in enumerate(VEIN_PATCHES.segment)
        Fm[m] > 0 || continue
        num[:, j0] .+= Fm[m] .* view(cm, :, m) # its own root segment carries its outflow
        den[j0] += Fm[m]
        p1 = tv.parent[j0]
        j = p1
        while j > 0
            d = PV.sv[p1] - PV.sv[j] + PV.tau[j] / 2
            for f in 1:NVF
                num[f, j] += Fm[m] * vlerp(view(cm, :, m), VT[f] - d)
            end
            den[j] += Fm[m]
            j = tv.parent[j]
        end
    end
    vconc = PF .* num ./ max.(den', 1.0e-30) # blood, mg/mL; (frames, segments)
    # the main arteries: each segment's arrival at its middle (a pass-1 terminal: at its start)
    ta = pass1.arteries.tree
    tbeg = [ta.parent[j] == 0 ? 0.0 : PA.tend[ta.parent[j]] for j in eachindex(ta.parent)]
    inner = falses(length(ta.parent))
    foreach(p -> p > 0 && (inner[p] = true), ta.parent)
    tmid = [inner[j] ? (tbeg[j] + PA.tend[j]) / 2 : tbeg[j] for j in eachindex(ta.parent)]
    # the patches' tissue: their vessels and beds per mL
    ml = [patch_ml(k) for k in ks]
    pconc = reduce(hcat, [(P[k].art .+ P[k].bed .+ P[k].ven) ./ patch_ml(k) for k in ks]) # (frames, patches)
    pos = pass1.arteries.tree.b[:, PATCHES.segment[ks]]
    sil = silhouette()
    D = (;
        vt = VT, aif = PF .* c[1:STRIDE:end], aif_dt = Float64(DT_OUT),
        art = (x1 = ta.a[1, :], z1 = ta.a[3, :], x2 = ta.b[1, :], z2 = ta.b[3, :], r = Float64.(PASS1_MURRAY.radius), tau = tmid),
        vein = (x1 = tv.a[1, :], z1 = tv.a[3, :], x2 = tv.b[1, :], z2 = tv.b[3, :], r = Float64.(PASS1_MURRAY_VEINS.radius), conc = vconc),
        tis = (x = pos[1, :], z = pos[3, :], ml, conc = pconc),
        sil, series = (art = art, tis = tis, ven = ven, total = total),
        injected = Float64(R.body.injected), leg_ml_min = 60 * R.body.share * R.body.co, npatch = length(ks),
    )
    save(VIEW_DATA, D)
    tsay(@sprintf("view data in %.0f s: %d main arteries, %d main veins, %d patches", time() - t0, length(D.art.x1), nv, length(ks)))
    return D
end

"The page's data from the view data."
function page_data(D)
    framemajor(v) = vec(permutedims(v)) # (frames, items) -> frame by frame
    tobytes(v, top) = UInt8.(round.(255 .* clamp.(v ./ top, 0, 1)))
    vmax = maximum(D.vein.conc)
    tmax = maximum(D.tis.conc)
    xs = vcat(D.art.x1, D.art.x2, D.vein.x1, D.vein.x2, [r[i] for r in D.sil.rows for i in 2:length(r)])
    zs = vcat(D.art.z1, D.art.z2, D.vein.z1, D.vein.z2, [r[1] for r in D.sil.rows])
    pad = 8.0
    return (;
        t = D.vt, aif = D.aif, aif_dt = D.aif_dt, aif_max = maximum(D.aif),
        art = D.art, vein = (; D.vein.x1, D.vein.z1, D.vein.x2, D.vein.z2, D.vein.r, q = tobytes(framemajor(D.vein.conc), vmax), max = vmax),
        tis = (; D.tis.x, D.tis.z, D.tis.ml, q = tobytes(framemajor(D.tis.conc), tmax), max = tmax),
        sil = D.sil, bounds = (x0 = minimum(xs) - pad, x1 = maximum(xs) + D.sil.h + pad, z0 = minimum(zs) - pad, z1 = maximum(zs) + D.sil.h + pad),
        series = (keys = ["artery", "tissue", "vein", "total"], names = ["Arteries", "Tissue", "Veins", "Total"],
                  values = [D.series.art, D.series.tis, D.series.ven, D.series.total]),
        meta = (
            title = "Contrast through the leg at rest",
            subtitle = "Simulated at rest: 35 mL of contrast (350 mg iodine per mL) at 5 mL/s into the right arm, then 35 mL of a 30:70 contrast:saline chaser. Press Play, drag the time, or click the chart.",
            legnote = @sprintf("Front view of the right leg. Lines: the main arteries and veins (the first pass, with XCAT's imaged femoral vessels), drawn at their radius. Dots: the %d tissue patches, by size. An artery carries the input delayed by its own arrival time; a vein carries the outflows of the patches it drains. Colours on a square-root scale.", D.npatch),
            chartnote = @sprintf("Arteries: the main and the patches' arteries. Tissue: arterioles, capillaries, the space between the cells, venules. Veins: the patches' and the main veins. Total: the whole leg, %.0f mg of the %.0f mg injected by %.0f s. The leg takes %.0f mL/min of blood at rest.", D.series.total[end], D.injected, D.vt[end], D.leg_ml_min),
        ),
    )
end

function main_view(redo_page_only = false)
    D = redo_page_only ? deserialize(VIEW_DATA) : view_data()
    path = write_view(VIEW_PAGE, page_data(D))
    tsay(@sprintf("page: %s (%.1f MB)", path, filesize(path) / 1.0e6))
    return path
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_view(get(ARGS, 1, "") == "page")
end
