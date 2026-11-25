# A narrowing in the leg's superficial femoral artery: the flow in every part of the leg, down to the foot, at rest and at
# maximal dilation. Built on flow_leg.jl's solved leg; the package is not changed.
#
#     julia -t auto stenosis_leg.jl test      the site, step 1 on a few patches, and step 2 on passive curves
#     julia -t auto stenosis_leg.jl           every step; patches already saved are skipped
#
# The narrowing: LENGTH_MM of XCAT's imaged femoral artery from its node nearest ABOVE_KNEE_MM above the knee joint line
# (halfway between XCAT's femur and tibia), below the deep femoral artery: the distal adductor canal, where superficial
# femoral stenoses are most common. Its loss is the package circuit's own (`circuit`'s `stenosis`): ΔP = f·Q + s·Q²
# (Young & Tsai 1973), viscous f = 8πμL_s/A_s², separation s = ρ·1.52/(2A₀²)·(A₀/A_s − 1)², counted once at the
# narrowing's exit. Its throat keeps its size at maximal dilation; the rest of its segments dilate as pass 1's do.
#
#   1. Curves, per patch: its flow against its inlet pressure, at its healthy outlet pressure. At rest from its healthy
#      circuit as the package's `reference`: the healthy tone constant is kept, so its arterioles answer a lower pressure
#      by dilating, and its flow falls only when they run out of dilation. At maximal dilation at its κ. At twelve inlet
#      pressures between its venous pressure and the femoral artery's (closer where the resting flow bends), and at its
#      healthy one. Saved as stenosis/curves/curve_<k>.jls.
#   2. The network, per stenosis (0–95% of the diameter) and state: pass 1 with the narrowing, every patch on its curve.
#      Solved exactly on the tree: each segment's flow as a function of its inlet pressure from the patches up, then the
#      pressures from the femoral artery down. The venous pressures are iterated, each patch's curve shifted by its outlet
#      pressure's change.
#   3. Exact, at 50, 70 and 90%: every patch whose pressures moved is solved again at them, the points join its curve and
#      the network is solved again, until no patch's pressure moves more than MOVED.
#   4. The results: each part of the leg's flow, its oxygen and its arterioles' tone, the trunk's pressures and flows.
#      Saved as stenosis/result.jls and summary.txt.
#
# Numbers only. Everything is made from XCAT, so it stays in LIMB_CACHE_DIR\stenosis (lab storage).

using Printf, Serialization, Statistics

include(joinpath(@__DIR__, "flow_leg.jl")) # the notebook and the flow's steps; included, it runs nothing

const S_DIR = mkpath(joinpath(CACHE_DIR, "stenosis"))
const CURVES = mkpath(joinpath(S_DIR, "curves"))
const SLOG = open(joinpath(S_DIR, "progress.log"), "a")
function ssay(s)
    lock(LOCK) do
        line = @sprintf("[%6.0f s] %s", time() - STARTED, s)
        println(line)
        flush(stdout)
        println(SLOG, line)
        flush(SLOG)
    end
    return nothing
end

const LENGTH_MM = 20.0                       # the narrowing's length
const ABOVE_KNEE_MM = 140.0                  # where it starts, above the knee joint line
const XS = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9, 1.0] # a curve's inlet pressures: p_out + x (P_ROOT - p_out), close where resting flow bends
const SEVERITIES = collect(0.0:0.05:0.95)    # diameter stenosis
const DETAIL = [SEVERITIES[argmin(abs.(SEVERITIES .- x))] for x in (0.5, 0.7, 0.9)] # solved exactly
const MOVED = 0.25                           # mmHg: a patch is solved again when its pressures move more
const P_ROOT = Float64(LIMB_FLOW.p_in)       # the femoral artery's pressure
const P_VEIN = Float64(LIMB_FLOW.p_out)      # the femoral vein's
const DP = 0.01                              # mmHg: the network's pressure grid
const PG = collect(0.0:DP:P_ROOT)
const NPG = length(PG)
const CV_SIGMA = sqrt(log(1 + LIMB_FLOW.radius_cv^2))

# ---------------------------------------------------------------------------------------------------------------------
# the narrowing and the parts of the leg

"The knee joint line (z, mm): halfway between the lowest point of XCAT's femur and the highest of its tibia."
const KNEE_Z = let ext(name) = (X = reduce(hcat, [first(mesh(s)) for s in leg_surfaces if group(s.name) == name]); extrema(view(X, 3, :)))
    (ext("r_femur")[1] + ext("r_tibia")[2]) / 2
end

const ART = pass1.arteries.tree
const NSEG = length(ART.parent)
const KIDS = let c = [Int[] for _ in 1:NSEG]
    foreach(j -> ART.parent[j] > 0 && push!(c[ART.parent[j]], j), 1:NSEG)
    c
end
seg_mm(j) = sqrt(sum(abs2, view(ART.b, :, j) .- view(ART.a, :, j)))

"The femoral trunk: pass 1's path from the femoral artery to the lowest end of XCAT's imaged arteries (femoral, popliteal, on to the foot)."
const TRUNK = let base = findall(ART.isbase), path = Int[]
    j = base[argmin(ART.b[3, base])]
    while j > 0
        pushfirst!(path, j)
        j = ART.parent[j]
    end
    path
end

"The narrowing: `LENGTH_MM` of the trunk from its node nearest `ABOVE_KNEE_MM` above the knee joint line, as (segment, narrowed length in mm)."
const SITE = let z0 = KNEE_Z + ABOVE_KNEE_MM, out = Tuple{Int, Float64}[], left = LENGTH_MM
    i = argmin(i -> abs(ART.a[3, TRUNK[i]] - z0), eachindex(TRUNK))
    while left > 1.0e-9
        l = min(seg_mm(TRUNK[i]), left)
        push!(out, (TRUNK[i], l))
        left -= l
        i += 1
    end
    out
end

"Distance along the trunk (mm) from the femoral artery to where it first reaches height `z`."
function trunk_dist(z)
    d = 0.0
    for j in TRUNK
        za, zb, L = ART.a[3, j], ART.b[3, j], seg_mm(j)
        (za - z) * (zb - z) <= 0 && return d + (za == zb ? 0.0 : (za - z) / (za - zb)) * L
        d += L
    end
    return NaN
end

"The narrowing's start and end: distance along the trunk from the femoral artery (mm) and height (z, mm)."
const SITE_SPAN = let i = findfirst(==(SITE[1][1]), TRUNK), from = sum(seg_mm, TRUNK[1:(i - 1)]; init = 0.0)
    j2, l2 = SITE[end]
    (from = from, to = from + LENGTH_MM, z_from = ART.a[3, SITE[1][1]], z_to = ART.a[3, j2] + (ART.b[3, j2] - ART.a[3, j2]) * l2 / seg_mm(j2))
end

"Whether each patch is fed through the narrowing."
const FED = map(PATCHES.segment) do j
    s = SITE[1][1]
    while j > 0 && j != s
        j = ART.parent[j]
    end
    j == s
end

"Each patch's tissue centre (x, z in mm)."
const TISSUE_XZ = let lab = PATCHES.label, g = LEG.grid, K = length(PATCHES.voxels)
    sx, sz, n = zeros(K), zeros(K), zeros(Int, K)
    for kz in axes(lab, 3), jy in axes(lab, 2), ix in axes(lab, 1)
        l = lab[ix, jy, kz]
        l > 0 || continue
        sx[l] += g.origin[1] + (ix - 1) * g.h
        sz[l] += g.origin[3] + (kz - 1) * g.h
        n[l] += 1
    end
    (x = sx ./ max.(n, 1), z = sz ./ max.(n, 1))
end

const PART_NAMES = ["Thigh, fed above the narrowing", "Thigh, fed through it", "Knee", "Calf", "Foot"]
const PART_SHORT = ["Thigh above", "Thigh below", "Knee", "Calf", "Foot"]
"""
Each patch's part of the leg: the foot (most of its tissue below the ankle, as the notebook's flow table), the calf and the
knee (tissue centre within 5 cm of the joint line) by height, and the thigh by whether the narrowing feeds it.
"""
const PART = map(eachindex(PATCHES.segment)) do k
    PATCH_MAX_SHARE[k] < (1 + FOOT_MAX_SHARE) / 2 ? 5 : TISSUE_XZ.z[k] < KNEE_Z - 50 ? 4 : TISSUE_XZ.z[k] < KNEE_Z + 50 ? 3 : FED[k] ? 2 : 1
end

"Pass-1 artery `j`'s radius (mm) with the circuit's lognormal factor, as flow_leg.jl's `resistances` takes it (before any dilation)."
lumen_r(j) = Float32(PASS1_MURRAY.radius[j]) * Float32(exp(CV_SIGMA * V.Flow._hashnormal(j + (1 << 41)) - CV_SIGMA^2 / 2))

"One narrowed segment as the package's circuit takes it: segment `j`, the narrowed fraction of its length, f (mmHg·s/mL), s (mmHg·s²/mL²)."
struct Narrow
    j::Int
    frac::Float64
    f::Float64
    s::Float64
end

"The narrowing at diameter stenosis `ds`, with the package circuit's formulas (src/perfusion/flow.jl, `_build`)."
function narrowing(ds)
    ds > 0 || return Narrow[]
    f = LIMB_FLOW
    area = Float32((1 - ds)^2)
    segs = Set(first.(SITE))
    return map(SITE) do (j, Ls)
        r = lumen_r(j)
        A0 = Float32(π) * (r * 1.0f-3)^2
        As = area * A0
        mu = f.plasma_viscosity_mPas * 1.0f-3 * V.Flow.mu_rel(2000.0f0 * r, f.hd)
        Lj = Float32(seg_mm(j))
        L = min(Float32(Ls), Lj)
        fv = 8.0f0 * Float32(π) * mu * L * 1.0f-3 / As^2 * 1.0f-6 / 133.322f0
        exit = !any(in(segs), KIDS[j])
        s2 = exit ? 1060.0f0 * 1.52f0 / (2A0^2) * (1 / area - 1)^2 * 1.0f-12 / 133.322f0 : 0.0f0
        Narrow(j, Float64(L / Lj), Float64(fv), Float64(s2))
    end
end

# ---------------------------------------------------------------------------------------------------------------------
# step 1: each patch's flow against its inlet pressure

flowspec(pin, pout) = V.FlowSpec(LIMB_FLOW; p_in = F32(pin), p_out = F32(pout))

"Patch `k`'s trees and beds, as flow_leg.jl's final step builds its circuits."
function patch_inputs(k, W)
    P = deserialize(joinpath(FLOW_DIR, "pair_$k.jls"))
    vleaf = deserialize(joinpath(FLOW_DIR, "vleaf_$k.jls"))
    art = arterial(k)
    terms = leaves(art)
    terr = fill(F32(P.ml / length(terms)), length(terms))
    vf = venous_forest(k, W)
    vf.m == P.m || error("patch $k: its venous patches differ from the pairing's")
    dem = V.Territory.demand_field(art.b[:, terms], terr)
    return (; k, P, vleaf, art, terr, vf, dem, roots = findall(iszero, vf.tree.parent))
end

"Patch `k`'s healthy outlet pressure at rest, as flow_leg.jl's final step sets it."
rest_pout(I, Br) = sum(I.vf.n .* Br.pv[I.vf.m]) / sum(I.vf.n)

"Patch `I` solved at `state` with inlet and outlet pressures `pin`, `pout` (mmHg), at its network κ."
solve_patch(I, N, state, pin, pout; reference = nothing) = V.Flow.circuit(
    I.art, I.terr; state, veins = I.vf.tree, vleaf = I.vleaf, phys = flowspec(pin, pout), kappa = F32(N.kappa[I.k]),
    demand = I.dem, reference,
)

"What a solved circuit tells about its patch."
function summary_of(C, I, pin, pout)
    r = C.resist
    return (;
        pin = Float64(pin), pout = Float64(pout), q = Float64(sum(r.root_flow)), qroot = Float64.(C.qv[I.roots]),
        so2 = Float64(r.so2_sinus), starved = Float64(r.starved), act = Float64(r.mean_act), dilated = Int(r.n_dilated),
        max_tone = Int(r.n_max_tone), reversed = reversed(C), iters = Int(r.iters),
    )
end

"Step 1 for patch `k` (not saved): its healthy rest state, and its flow at rest and at maximal dilation against its inlet pressure."
function curve_core(k, W, N, Br)
    t0 = time()
    I = patch_inputs(k, W)
    po_r, po_s = rest_pout(I, Br), N.pout[k]
    Cr = solve_patch(I, N, :rest, Br.pin[k], po_r)
    healthy = summary_of(Cr, I, Br.pin[k], po_r)
    rest = [summary_of(solve_patch(I, N, :rest, p, po_r; reference = Cr), I, p, po_r) for p in sort!(unique([Br.pin[k]; po_r .+ XS .* (P_ROOT - po_r)]))]
    stress = [summary_of(solve_patch(I, N, :stress, p, po_s), I, p, po_s) for p in sort!(unique([N.pin[k]; po_s .+ XS .* (P_ROOT - po_s)]))]
    return (;
        k, healthy, rest, stress, pin_rest = Br.pin[k], pout_rest = po_r, pin_stress = N.pin[k], pout_stress = po_s,
        seconds = time() - t0,
    )
end

curve_patch(k, W, N, Br) = save(joinpath(CURVES, "curve_$k.jls"), curve_core(k, W, N, Br))

"Run `work(k)` for every `k` in `ks`, `CONCURRENT` at a time, with progress every 50 patches; failures are listed, not fatal."
function run_each(work, ks, label)
    isempty(ks) && return Int[]
    ch = Channel{Int}(length(ks))
    foreach(k -> put!(ch, k), ks)
    close(ch)
    done = Threads.Atomic{Int}(0)
    failed = Int[]
    t0 = time()
    @sync for _ in 1:min(CONCURRENT, length(ks))
        Threads.@spawn for k in ch
            try
                work(k)
            catch e
                lock(LOCK) do
                    push!(failed, k)
                end
                ssay("$label: patch $k failed: " * first(sprint(showerror, e), 300))
            end
            n = Threads.atomic_add!(done, 1) + 1
            if n % 50 == 0 || n == length(ks)
                dt = time() - t0
                ssay(@sprintf("%s: %d of %d patches in %.0f s, about %.0f min left", label, n, length(ks), dt, dt / n * (length(ks) - n) / 60))
            end
        end
    end
    return failed
end

# ---------------------------------------------------------------------------------------------------------------------
# step 2: the network

"A patch's flow (mL/s) against its inlet pressure (mmHg): from its outlet pressure, where it carries nothing, linear between points and beyond the last."
struct Curve
    p::Vector{Float64}
    q::Vector{Float64}
end
Curve(pout, points) = (p = [pout; [x.pin for x in points]]; q = [0.0; [x.q for x in points]]; o = sortperm(p); Curve(p[o], q[o]))
function (c::Curve)(p)
    p <= c.p[1] && return 0.0
    i = min(searchsortedlast(c.p, p), length(c.p) - 1)
    return max(c.q[i] + (c.q[i + 1] - c.q[i]) * (p - c.p[i]) / (c.p[i + 1] - c.p[i]), 0.0)
end
"Add the point (`p`, `q`) to curve `c`, replacing one within 1e-6 mmHg."
function add!(c::Curve, p, q)
    i = findfirst(x -> abs(x - p) < 1.0e-6, c.p)
    i === nothing || (c.q[i] = q; return c)
    i = searchsortedfirst(c.p, p)
    insert!(c.p, i, p)
    insert!(c.q, i, q)
    return c
end

"Pass 1's arterial losses per segment at `state` with the narrowing: ΔP = a Q + s Q²."
function losses(state, sten)
    a = copy(state === :rest ? PASS1.ra_rest : PASS1.ra_stress)
    s = zeros(NSEG)
    for x in sten
        a[x.j] = a[x.j] * (1 - x.frac) + x.f
        s[x.j] = x.s
    end
    return a, s
end

const PATCH_AT = Dict(j => k for (k, j) in enumerate(PATCHES.segment))

"Linear interpolation of `g` (on `PG`) at pressure `p`."
at_grid(g, p) = (u = clamp(p / DP, 0.0, NPG - 1.0); i = min(floor(Int, u), NPG - 2); g[i + 1] + (u - i) * (g[i + 2] - g[i + 1]))

"""
Pass 1's arteries with patch `k` taking `flow(k, p)` at inlet pressure p, and losses ΔP = a Q + s Q² per segment: every
segment's inlet pressure and flow. Exact on the tree: each segment's flow as a function of its inlet pressure (on `PG`),
from the patches up; then the pressures from the femoral artery (`P_ROOT`) down. `G` is a workspace (NPG × segments).
"""
function solve_arteries!(G, flow::F, a, s) where {F}
    H, env = zeros(NPG), zeros(NPG)
    for j in Iterators.reverse(PASS1.a_order)
        k = get(PATCH_AT, j, 0)
        if k > 0
            for i in 1:NPG
                H[i] = flow(k, PG[i])
            end
        else
            fill!(H, 0.0)
            for c in KIDS[j]
                H .+= view(G, :, c)
            end
        end
        Gj = view(G, :, j)
        if a[j] == 0 && s[j] == 0
            Gj .= H
            continue
        end
        # the segment's inlet pressure when its end is at PG[i]; inverted onto the grid
        m = -Inf
        for i in 1:NPG
            m = max(m, PG[i] + a[j] * H[i] + s[j] * H[i]^2)
            env[i] = m
        end
        i = 1
        for n in 1:NPG
            P = PG[n]
            while i < NPG && env[i + 1] <= P
                i += 1
            end
            if i == NPG
                Gj[n] = H[NPG]
            else
                d = env[i + 1] - env[i]
                Gj[n] = H[i] + (d > 0 ? clamp((P - env[i]) / d, 0.0, 1.0) : 0.0) * (H[i + 1] - H[i])
            end
        end
    end
    pin, q = zeros(NSEG), zeros(NSEG)
    for j in PASS1.a_order
        p = ART.parent[j]
        pin[j] = p == 0 ? P_ROOT : pin[p] - a[p] * q[p] - s[p] * q[p]^2
        q[j] = at_grid(view(G, :, j), pin[j])
    end
    return pin, q
end

"Each patch's outlet pressure, pass 1's veins carrying `split[k]` of patch `k`'s flow `Q[k]` into its venous patches `ms[k]`, weighted by `weight[k]`."
function venous_pout(Q, ms, split, weight)
    qv = zeros(length(VEIN_PATCHES.segment))
    for k in eachindex(Q), (i, m) in enumerate(ms[k])
        qv[m] += split[k][i] * Q[k]
    end
    Pv = pass_pressures(pass1.veins.tree, PASS1.v_order, PASS1.rv, VEIN_PATCHES.segment, qv, P_VEIN, +1)[VEIN_PATCHES.segment]
    return [isempty(ms[k]) ? NaN : sum(weight[k] .* Pv[ms[k]]) for k in eachindex(Q)]
end

"""
The network's patch data at `state`: each patch's curve (`curve(k)`), its healthy pressures, and how its outflow splits
among its venous patches and weighs their pressures, as flow_leg.jl set them (maximal: the pairing's weights for both;
rest: its healthy venous outflows, and its venules' counts). `PS` and `FS` hold the patches' pairing and final summaries.
"""
function network_data(state, ks, PS, FS, curve)
    K = length(PATCHES.segment)
    ok = falses(K)
    curves = [Curve([0.0, 1.0], [0.0, 0.0]) for _ in 1:K]
    pout_h, pin_h = zeros(K), zeros(K)
    ms, split, weight = [Int[] for _ in 1:K], [Float64[] for _ in 1:K], [Float64[] for _ in 1:K]
    for k in ks
        P, F = PS[k], FS[k]
        ms[k] = P.m
        if state === :rest
            split[k] = F.rest.qroot ./ sum(F.rest.qroot)
            weight[k] = P.n ./ sum(P.n)
            pout_h[k], pin_h[k] = F.pout_rest, F.pin_rest
        else
            split[k] = P.w
            weight[k] = P.w
            pout_h[k], pin_h[k] = F.pout, F.pin
        end
        curves[k] = curve(k)
        ok[k] = true
    end
    return (; state, ok, curves, pout_h, pin_h, ms, split, weight)
end

"""
Pass 1 at `D.state` with the narrowing `sten`, every patch on its curve (shifted by its outlet pressure's change): each patch's
inlet and outlet pressure and flow, and pass 1's arterial inlet pressures and flows.
"""
function solve_network(G, D, sten; tol = 1.0e-6)
    a, s = losses(D.state, sten)
    K = length(D.ok)
    shift = zeros(K)
    flow(k, p) = D.ok[k] ? D.curves[k](p - shift[k]) : 0.0
    pin_seg, q_seg = zeros(NSEG), zeros(NSEG)
    Q, pout, it, d = zeros(K), zeros(K), 0, Inf
    for outer in 1:100
        it = outer
        pin_seg, q_seg = solve_arteries!(G, flow, a, s)
        for k in 1:K
            Q[k] = flow(k, pin_seg[PATCHES.segment[k]])
        end
        pout = venous_pout(Q, D.ms, D.split, D.weight)
        d = 0.0
        for k in 1:K
            D.ok[k] || continue
            new = pout[k] - D.pout_h[k]
            d = max(d, abs(new - shift[k]))
            shift[k] = new
        end
        d < tol && break
    end
    # the patches at the final shifts
    pin_seg, q_seg = solve_arteries!(G, flow, a, s)
    for k in 1:K
        Q[k] = flow(k, pin_seg[PATCHES.segment[k]])
    end
    pin = pin_seg[PATCHES.segment]
    root = PASS1.a_order[1]
    return (; pin, pout, Q, pin_seg, q_seg, pend_seg = pin_seg .- a .* q_seg .- s .* q_seg .^ 2, iters = it, change = d, balance = q_seg[root] / sum(Q) - 1)
end

# ---------------------------------------------------------------------------------------------------------------------
# step 3: exact at the chosen stenoses

"Patch `k` solved again at each (state, ds, pin, pout) of `todo`: at rest from its healthy circuit as the reference."
function exact_core(k, W, N, Br, todo)
    I = patch_inputs(k, W)
    Cr = any(t -> t.state === :rest, todo) ? solve_patch(I, N, :rest, Br.pin[k], rest_pout(I, Br)) : nothing
    return [(; t.state, t.ds, x = summary_of(solve_patch(I, N, t.state, t.pin, t.pout; reference = t.state === :rest ? Cr : nothing), I, t.pin, t.pout)) for t in todo]
end

# ---------------------------------------------------------------------------------------------------------------------
# the trunk

"Each trunk segment's distance along the trunk from the femoral artery, at its end (mm)."
const TRUNK_MM = cumsum(seg_mm.(TRUNK))

"The trunk's pressure at height `z` (its first crossing, going down), from its segments' inlet and end pressures."
function trunk_at(z, pin_seg, pend_seg)
    for (i, j) in enumerate(TRUNK)
        za, zb = ART.a[3, j], ART.b[3, j]
        if zb <= z <= za || zb >= z >= za
            t = za == zb ? 0.0 : (z - za) / (zb - za)
            return pin_seg[j] + t * (pend_seg[j] - pin_seg[j])
        end
    end
    return NaN
end

"The front-view silhouette of the leg's soft tissue, rows of runs on cells `2h` mm (transport_view.jl's `silhouette`)."
function leg_silhouette()
    soft, g = LEG.soft, LEG.grid
    nx, nz = size(soft, 1), size(soft, 3)
    proj = falses(cld(nx, 2), cld(nz, 2))
    for kz in 1:nz, jy in axes(soft, 2), ix in 1:nx
        soft[ix, jy, kz] && (proj[cld(ix, 2), cld(kz, 2)] = true)
    end
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
    return (; h = 2 * g.h, rows)
end

"What the network says about one stenosis at one state."
function network_summary(S, ml)
    first_j, last_j = SITE[1][1], SITE[end][1]
    parts = [findall(==(p), PART) for p in eachindex(PART_NAMES)]
    return (;
        S.pin, S.pout, S.Q, iters = S.iters, change = S.change, balance = S.balance,
        trunk_pin = S.pin_seg[TRUNK], trunk_pend = S.pend_seg[TRUNK], trunk_q = S.q_seg[TRUNK],
        above = S.pin_seg[first_j], below = S.pend_seg[last_j], through = S.q_seg[first_j],
        knee = trunk_at(KNEE_Z, S.pin_seg, S.pend_seg), ankle = trunk_at(ANKLE.z, S.pin_seg, S.pend_seg),
        part_q = [sum(S.Q[ks]) for ks in parts], part_ml = [sum(ml[ks]) for ks in parts], leg_q = sum(S.Q),
    )
end

# ---------------------------------------------------------------------------------------------------------------------
# steps 3 and 4

"""
Step 3: the network at each stenosis of `DETAIL`, every patch whose pressures moved by more than `MOVED` (from its healthy
ones, or its last exact solve) solved again at them, the points added to its curves, until none moves. Returns the exact
solves by (patch, state, stenosis) and, per round, how far the curves were from them.
"""
function exact_rounds!(G, D, ks, W, N, Br)
    path = joinpath(S_DIR, "exact.jls")
    if isfile(path)
        X = deserialize(path)
        for (k, state, p, q) in X.added
            add!(D[state].curves[k], p, q)
        end
        return X.E, X.misses
    end
    E = Dict{Tuple{Int, Symbol, Float64}, Any}()
    misses = Any[]
    added = Tuple{Int, Symbol, Float64, Float64}[]
    for round in 1:4
        todo = Dict{Int, Vector{Any}}()
        for state in (:rest, :stress), ds in DETAIL
            S = solve_network(G, D[state], narrowing(ds))
            for k in ks
                last = get(E, (k, state, ds), nothing)
                p0, o0 = last === nothing ? (D[state].pin_h[k], D[state].pout_h[k]) : (last.pin, last.pout)
                if abs(S.pin[k] - p0) > MOVED || abs(S.pout[k] - o0) > MOVED
                    push!(get!(todo, k, Any[]), (; state, ds, pin = S.pin[k], pout = S.pout[k]))
                end
            end
        end
        n = sum(length, values(todo); init = 0)
        ssay(@sprintf("exact, round %d: %d solves in %d patches", round, n, length(todo)))
        if n == 0
            save(path, (; E, misses, added))
            return E, misses
        end
        out = Dict{Int, Any}()
        failed = run_each(collect(keys(todo)), "exact $round") do k
            r = exact_core(k, W, N, Br, todo[k])
            lock(LOCK) do
                out[k] = r
            end
        end
        isempty(failed) || error("exact solves failed for patches $(failed)")
        miss = Float64[]
        for (k, rs) in out, r in rs
            c = D[r.state].curves[k]
            p = r.x.pin - (r.x.pout - D[r.state].pout_h[k]) # on the curve's outlet pressure
            r.x.q > 0 && push!(miss, c(p) / r.x.q - 1)
            add!(c, p, r.x.q)
            push!(added, (k, r.state, p, r.x.q))
            E[(k, r.state, r.ds)] = r.x
        end
        push!(misses, miss)
        ssay(@sprintf("exact, round %d: the curves were off by median %.2g%%, 99%% %.2g%%, max %.2g%%", round, 100 * median(abs.(miss)), 100 * quantile(abs.(miss), 0.99), 100 * maximum(abs.(miss))))
    end
    ssay("exact: pressures still moving after 4 rounds")
    save(path, (; E, misses, added))
    return E, misses
end

"Patch `k`'s rest or maximal state at stenosis `ds`: its last exact solve, or its healthy state where its pressures did not move."
function patch_state(state, ds, k, CS, E)
    x = get(E, (k, state, ds), nothing)
    x === nothing || return x
    c = CS[k]
    return state === :rest ? c.healthy : c.stress[findfirst(p -> p.pin == c.pin_stress, c.stress)]
end

"Each part's resting oxygen and tone at stenosis `ds` (`ds` 0: healthy): venous SO₂ by flow, beds starved of O₂, mean activation by tissue."
function part_oxygen(ds, ks, CS, E, FS)
    xs = Dict(k => patch_state(:rest, ds, k, CS, E) for k in ks)
    return map(eachindex(PART_NAMES)) do p
        kp = [k for k in ks if PART[k] == p]
        q, ml, nt = [xs[k].q for k in kp], [FS[k].ml for k in kp], [FS[k].nterm for k in kp]
        (;
            so2 = sum(q .* [xs[k].so2 for k in kp]) / sum(q), starved = sum(nt .* [xs[k].starved for k in kp]) / sum(nt),
            act = sum(ml .* [xs[k].act for k in kp]) / sum(ml), q = sum(q),
        )
    end
end

"The site, the parts, and the narrowing's loss at the healthy flows through it."
function describe_site(ks, FS)
    ml = [k in ks ? FS[k].ml : 0.0 for k in eachindex(PATCHES.segment)]
    S = SITE_SPAN
    ssay(@sprintf(
        "knee joint line z %.0f mm, ankle z %.0f mm; the narrowing: %s, z %.0f to %.0f mm (its middle %.0f mm above the knee joint line), %.0f to %.0f mm along the trunk from the femoral artery (knee joint at %.0f, ankle at %.0f mm)",
        KNEE_Z, ANKLE.z, join([@sprintf("segment %d %.1f of %.1f mm, lumen %.2f mm across", j, l, seg_mm(j), 2lumen_r(j)) for (j, l) in SITE], "; "),
        S.z_from, S.z_to, (S.z_from + S.z_to) / 2 - KNEE_Z, S.from, S.to, trunk_dist(KNEE_Z), trunk_dist(ANKLE.z),
    ))
    ssay(@sprintf("fed through the narrowing: %d of %d patches, %.0f of %.0f mL (%.1f%%)", count(FED[ks]), length(ks), sum(ml[FED]), sum(ml), 100 * sum(ml[FED]) / sum(ml)))
    for p in eachindex(PART_NAMES)
        kp = [k for k in ks if PART[k] == p]
        ssay(@sprintf("  %-32s %4d patches, %5.0f mL, %5.1f%% fed through the narrowing", PART_NAMES[p], length(kp), sum(ml[kp]), 100 * sum(ml[kp] .* FED[kp]) / sum(ml[kp])))
    end
    qr = sum(FS[k].rest.q for k in ks if FED[k])
    qs = sum(FS[k].stress.q for k in ks if FED[k])
    for ds in DETAIL
        x = narrowing(ds)
        f, s = sum(n.f for n in x), sum(n.s for n in x)
        ssay(@sprintf("  %.0f%%: f %.3g mmHg·s/mL, s %.3g mmHg·s²/mL²; at the healthy flow through it, the loss would be %.1f mmHg at rest (%.2f mL/s), %.0f mmHg at maximal dilation (%.1f mL/s)", 100ds, f, s, f * qr + s * qr^2, qr, f * qs + s * qs^2, qs))
    end
    return nothing
end

"Load the patches' pairing and final summaries."
function patch_summaries(ks)
    PS = Dict(k => deserialize(joinpath(FLOW_DIR, "pair_$k.jls")) for k in ks)
    FS = Dict(k => deserialize(joinpath(FLOW_DIR, "final_$k.jls")) for k in ks)
    return PS, FS
end

solved_patches() = [k for k in eachindex(PATCHES.segment) if isfile(joinpath(FLOW_DIR, "final_$k.jls"))]

function run_stenosis()
    W = venule_map()
    N = deserialize(joinpath(FLOW_DIR, "network.jls"))
    Br = equal_share(:rest)
    ks = solved_patches()
    PS, FS = patch_summaries(ks)
    describe_site(ks, FS)
    # step 1
    todo = filter(k -> !isfile(joinpath(CURVES, "curve_$k.jls")), ks)
    failed = run_each(k -> curve_patch(k, W, N, Br), todo, "curves")
    isempty(failed) || error("curves failed for patches $(failed)")
    CS = Dict(k => deserialize(joinpath(CURVES, "curve_$k.jls")) for k in ks)
    rebuilt = [CS[k].healthy.q / FS[k].rest.q - 1 for k in ks]
    kept = [CS[k].rest[findfirst(p -> p.pin == CS[k].pin_rest, CS[k].rest)].q / CS[k].healthy.q - 1 for k in ks]
    stress0 = [CS[k].stress[findfirst(p -> p.pin == CS[k].pin_stress, CS[k].stress)].q / FS[k].stress.q - 1 for k in ks]
    ssay(@sprintf(
        "curves: healthy rest rebuilt max |%.2g|; rest with the healthy tone kept, at the healthy pressure, max |%.2g| (median %.2g); maximal rebuilt max |%.2g|",
        maximum(abs, rebuilt), maximum(abs, kept), median(abs.(kept)), maximum(abs, stress0),
    ))
    falls = count(k -> any(i -> CS[k].rest[i + 1].q < CS[k].rest[i].q, 1:(length(CS[k].rest) - 1)), ks)
    ssay(@sprintf("curves: %d patches with a resting flow that falls somewhere as the pressure rises", falls))
    # steps 2 and 3
    D = Dict(
        :rest => network_data(:rest, ks, PS, FS, k -> Curve(CS[k].pout_rest, CS[k].rest)),
        :stress => network_data(:stress, ks, PS, FS, k -> Curve(CS[k].pout_stress, CS[k].stress)),
    )
    G = Matrix{Float64}(undef, NPG, NSEG)
    E, misses = exact_rounds!(G, D, ks, W, N, Br)
    ml = [k in ks ? FS[k].ml : 0.0 for k in eachindex(PATCHES.segment)]
    nets = Dict{Tuple{Symbol, Float64}, Any}()
    for state in (:rest, :stress), ds in SEVERITIES
        nets[(state, ds)] = network_summary(solve_network(G, D[state], narrowing(ds)), ml)
    end
    oxygen = Dict(ds => part_oxygen(ds, ks, CS, E, FS) for ds in [0.0; DETAIL])
    report(nets, oxygen, D, E, misses, CS, FS, ks, ml)
    ssay(@sprintf("done in %.1f h", (time() - STARTED) / 3600))
end

"Step 4: result.jls and summary.txt."
function report(nets, oxygen, D, E, misses, CS, FS, ks, ml)
    lines = String[]
    out(s) = (push!(lines, s); ssay(s))
    h_r, h_s = nets[(:rest, 0.0)], nets[(:stress, 0.0)]
    out(@sprintf(
        "healthy, by the network on the curves: rest %.1f mL/min (flow_leg.jl %.1f), maximal %.0f (%.0f); inlet pressures moved by max %.3f mmHg at rest, %.3f at maximal",
        60 * h_r.leg_q, 60 * sum(FS[k].rest.q for k in ks), 60 * h_s.leg_q, 60 * sum(FS[k].stress.q for k in ks),
        maximum(k -> abs(h_r.pin[k] - D[:rest].pin_h[k]), ks), maximum(k -> abs(h_s.pin[k] - D[:stress].pin_h[k]), ks),
    ))
    out(@sprintf("networks: at most %d iterations, venous pressures settled to %.1g mmHg, flow balance within %.1g", maximum(n.iters for n in values(nets)), maximum(n.change for n in values(nets)), maximum(abs(n.balance) for n in values(nets))))
    out("")
    for (state, label) in ((:rest, "REST"), (:stress, "MAXIMAL DILATION"))
        out("$label: flow of each part, % of healthy (mL/min healthy in brackets)")
        out(@sprintf("%-10s %s %10s", "stenosis", join([@sprintf("%16s", n) for n in PART_SHORT]), "leg"))
        h = nets[(state, 0.0)]
        out(@sprintf("%-10s %s %10s", "healthy", join([@sprintf("%16s", @sprintf("(%.1f)", 60 * q)) for q in h.part_q]), @sprintf("(%.0f)", 60 * h.leg_q)))
        for ds in SEVERITIES[2:end]
            n = nets[(state, ds)]
            out(@sprintf("%-10s %s %9.1f%%", @sprintf("%.0f%%", 100ds), join([@sprintf("%15.1f%%", 100 * q / q0) for (q, q0) in zip(n.part_q, h.part_q)]), 100 * n.leg_q / h.leg_q))
        end
        out("$label: pressures (mmHg) and the flow through the narrowing")
        out(@sprintf("%-10s %10s %10s %10s %10s %10s %12s", "stenosis", "above", "below", "knee", "ankle", "ankle/90", "through mL/min"))
        for ds in SEVERITIES
            n = nets[(state, ds)]
            out(@sprintf("%-10s %10.1f %10.1f %10.1f %10.1f %10.2f %12.1f", @sprintf("%.0f%%", 100ds), n.above, n.below, n.knee, n.ankle, n.ankle / P_ROOT, 60 * n.through))
        end
        out("")
    end
    out("REST: oxygen and tone of each part (venous SO2 by flow; beds starved of O2; mean arteriolar activation, 0 = fully dilated)")
    for ds in [0.0; DETAIL]
        out(@sprintf("%-10s %s", ds == 0 ? "healthy" : @sprintf("%.0f%%", 100ds), join([@sprintf("  %s: SO2 %.3f, starved %.1f%%, act %.2f", PART_SHORT[p], o.so2, 100 * o.starved, o.act) for (p, o) in enumerate(oxygen[ds])], ";")))
    end
    out(@sprintf("exact solves: %d in %d rounds; curves off by max %s", length(E), length(misses), join([@sprintf("%.2g%%", 100 * maximum(abs.(m); init = 0.0)) for m in misses], ", then ")))
    open(joinpath(S_DIR, "summary.txt"), "w") do io
        foreach(l -> println(io, l), lines)
    end
    curves = Dict(state => [D[state].ok[k] ? (p = D[state].curves[k].p, q = D[state].curves[k].q, pout = D[state].pout_h[k]) : nothing for k in eachindex(PATCHES.segment)] for state in (:rest, :stress))
    detail = Dict((ds, k) => patch_state(:rest, ds, k, CS, E) for ds in DETAIL for k in ks)
    save(joinpath(S_DIR, "result.jls"), (;
        severities = SEVERITIES, detail_severities = DETAIL, part_names = PART_NAMES, part = PART, fed = FED, ks, ml,
        nterm = [k in ks ? FS[k].nterm : 0 for k in eachindex(PATCHES.segment)], tissue_x = TISSUE_XZ.x, tissue_z = TISSUE_XZ.z,
        pos = ART.b[:, PATCHES.segment], knee_z = KNEE_Z, ankle_z = ANKLE.z, p_root = P_ROOT, p_vein = P_VEIN,
        site = [(j = j, mm = l, za = ART.a[3, j], zb = ART.b[3, j], xa = ART.a[1, j], xb = ART.b[1, j], d = 2lumen_r(j)) for (j, l) in SITE],
        trunk = TRUNK, trunk_mm = TRUNK_MM, trunk_x = ART.b[1, TRUNK], trunk_z = ART.b[3, TRUNK], trunk_za = ART.a[3, TRUNK],
        part_short = PART_SHORT, site_from = SITE_SPAN.from, site_to = SITE_SPAN.to,
        site_above_knee = (SITE_SPAN.z_from + SITE_SPAN.z_to) / 2 - KNEE_Z, knee_mm = trunk_dist(KNEE_Z), ankle_mm = trunk_dist(ANKLE.z),
        nets, oxygen, curves, detail, healthy_rest = Dict(k => CS[k].healthy for k in ks), misses,
        sil = leg_silhouette(), art = (x1 = ART.a[1, :], z1 = ART.a[3, :], x2 = ART.b[1, :], z2 = ART.b[3, :], r = Float64.(PASS1_MURRAY.radius)),
    ))
    return nothing
end

# ---------------------------------------------------------------------------------------------------------------------

"The site; step 1 on one patch per part (one at a time, then all at once); step 2 on passive curves through each patch's healthy state."
function test_stenosis()
    W = venule_map()
    N = deserialize(joinpath(FLOW_DIR, "network.jls"))
    Br = equal_share(:rest)
    ks = solved_patches()
    PS, FS = patch_summaries(ks)
    describe_site(ks, FS)
    picks = map(eachindex(PART_NAMES)) do p
        kp = [k for k in ks if PART[k] == p]
        mid = median(FS[k].ml for k in kp)
        kp[argmin([abs(FS[k].ml - mid) for k in kp])]
    end
    ssay("test patches, one per part: $(picks)")
    t1 = @elapsed one = [curve_core(k, W, N, Br) for k in picks]
    for c in one
        k = c.k
        ssay(@sprintf(
            "patch %d (%s, %.1f mL, %d terminals), %.0f s: healthy rest %.4f mL/s rebuilt %+.2g; kept tone at the healthy pressure %+.2g; maximal at the healthy pressure %+.2g",
            k, PART_NAMES[PART[k]], FS[k].ml, FS[k].nterm, c.seconds, c.healthy.q, c.healthy.q / FS[k].rest.q - 1,
            c.rest[findfirst(p -> p.pin == c.pin_rest, c.rest)].q / c.healthy.q - 1,
            c.stress[findfirst(p -> p.pin == c.pin_stress, c.stress)].q / FS[k].stress.q - 1,
        ))
        ssay("   rest:    " * join([@sprintf("%.1f→%.3f (SO2 %.2f, act %.2f, %d it)", p.pin, p.q / c.healthy.q, p.so2, p.act, p.iters) for p in c.rest], ", "))
        ssay("   maximal: " * join([@sprintf("%.1f→%.3f", p.pin, p.q / FS[k].stress.q) for p in c.stress], ", "))
    end
    many = Vector{Any}(undef, length(picks))
    t2 = @elapsed @sync for (i, k) in enumerate(picks)
        Threads.@spawn many[i] = curve_core(k, W, N, Br)
    end
    same = all(i -> all(j -> one[i].rest[j].q == many[i].rest[j].q, eachindex(one[i].rest)) && all(j -> one[i].stress[j].q == many[i].stress[j].q, eachindex(one[i].stress)), eachindex(picks))
    ssay(@sprintf("step 1: %.0f s one at a time, %.0f s all at once (%.1f× faster); identical: %s", t1, t2, t1 / t2, same))
    # step 2 on passive curves
    G = Matrix{Float64}(undef, NPG, NSEG)
    ml = [k in ks ? FS[k].ml : 0.0 for k in eachindex(PATCHES.segment)]
    for state in (:rest, :stress)
        D = network_data(state, ks, PS, FS, k -> state === :rest ? Curve(FS[k].pout_rest, [(pin = FS[k].pin_rest, q = FS[k].rest.q)]) : Curve(FS[k].pout, [(pin = FS[k].pin, q = FS[k].stress.q)]))
        for ds in (0.0, 0.5, 0.7, 0.9)
            t = @elapsed S = solve_network(G, D, narrowing(ds))
            n = network_summary(S, ml)
            ssay(@sprintf(
                "%s %.0f%%: leg %.1f mL/min, parts %s mL/min; above %.1f below %.1f knee %.1f ankle %.1f mmHg; %d iterations (%.1g), balance %.1g; inlet pressures moved max %.3f mmHg from flow_leg.jl's; %.1f s",
                state, 100ds, 60 * n.leg_q, join([@sprintf("%.1f", 60q) for q in n.part_q], "/"), n.above, n.below, n.knee, n.ankle,
                n.iters, n.change, n.balance, maximum(k -> abs(S.pin[k] - D.pin_h[k]), ks), t,
            ))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    get(ARGS, 1, "") == "test" ? test_stenosis() : run_stenosis()
end
