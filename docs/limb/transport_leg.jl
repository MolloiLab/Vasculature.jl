# The leg's contrast transport at rest (limb_pipeline.jl, section 7, step 3): the
# injection (`INJECTION`) through the whole-body model, with the grown leg in place of the body's right leg.
#
#     julia -t auto transport_leg.jl          the whole leg (steps A-D); steps already saved are skipped
#     julia -t auto transport_leg.jl test     a few patches, checked against the package's own frame_masses and sinus_outflow
#
# The package is not changed. Its `transport` (src/perfusion/transport.jl) takes only the heart or the lungs as the organ
# (`ORGAN_SLOTS`). This is the same loop with the leg in the body's right-leg slot, calling the package's body model
# (`Body.human_circulation`, `Body.simulate` with `detailed`) and its circuit as they are. The leg is solved patch by
# patch, as flow_leg.jl solves its flow:
#   A. Venous flows. A patch's venous trees are the venous patches' trees cut down to its venules, so a vein that drains
#      several arterial patches is in each of their circuits, carrying only that patch's flow. Its transit is its plasma
#      volume over the plasma flow of every bed it drains, summed here over the patches. Saved as venous_flows.jls.
#   B. Per patch, the resting circuit exactly as flow_leg.jl's final step solves it, with those venous transits, and pass
#      1's transits at the resting flows: to the patch's inlet, and from each venous patch's root to the femoral vein.
#      Every bed's arrival times are saved as beds_<k>.jls.
#   C. The leg's step responses: every bed's arrival at the patch's inlet, its own inlet and outlet, its venous patch's
#      root and the femoral vein, binned on a 0.02 s grid (its back-flux rate on the package's InputCurve grid). The leg
#      is linear, so its response to any arterial input is a convolution with them.
#   D. The body: `human_circulation` with the right leg's share of the cardiac output set to the grown leg's resting flow
#      (the other regions scaled to keep the total), the leg fed from the abdominal aorta. The leg's femoral venous
#      outflow replaces the body's right-leg beds until it changes by at most `tol`, as `transport` iterates. Then the
#      iodine in every compartment over time, by region, and the checks. Saved as result.jls and summary.txt.
#
# Numbers only. Everything is made from XCAT, so it stays in LIMB_CACHE_DIR\transport_rest (lab storage).

using Printf, Serialization, Statistics

include(joinpath(@__DIR__, "flow_leg.jl")) # the notebook and the flow's steps; included, it runs nothing

const T_DIR = mkpath(joinpath(CACHE_DIR, "transport_rest"))
const TLOG = open(joinpath(T_DIR, "progress.log"), "a")
function tsay(s)
    lock(LOCK) do
        line = @sprintf("[%6.0f s] %s", time() - STARTED, s)
        println(line)
        flush(stdout)
        println(TLOG, line)
        flush(TLOG)
    end
    return nothing
end

const HD = LIMB_FLOW.hd
const PF = 1 - HD                 # plasma fraction of the blood in large vessels
const TMAX = 120.0f0              # s: at rest the leg's blood volume over its flow is over a minute
const FINE = 0.02f0               # s: the input grid, as `transport`'s `fine`
const DT_OUT = 0.1f0              # s: the venous outflow's grid, as `transport`'s `dt_out`
const GRID = collect(0.0f0:FINE:TMAX)
const NG = length(GRID)
const TC = collect(0.0f0:DT_OUT:TMAX)
const STRIDE = round(Int, DT_OUT / FINE)
const LOGK = collect(Float64, range(log(1.0e-4), log(10.0), length = 160)) # V.Flow.InputCurve's back-flux grid
const NK = length(LOGK)
const DLOGK = LOGK[2] - LOGK[1]
const SIG = sqrt(log(1 + LIMB_FLOW.radius_cv^2))
const REST = (co = 81.6f0, beat = 1.0f0, f_cor = 0.033f0) # HEART_TRANSPORT's resting values
const LEG_BEDS = [Symbol(t, :_leg_r) for t in (:muscle, :skin, :fat, :bone, :other)]
const REGIONS = ("foot", "calf", "knee", "thigh")
"Patch `k`'s region, as the notebook's flow table: the foot (most of its tissue below the ankle), then by height."
region(k) = PATCH_MAX_SHARE[k] < (1 + FOOT_MAX_SHARE) / 2 ? 1 : HEIGHT[k] < 0.4 ? 2 : HEIGHT[k] < 0.7 ? 3 : 4

"Run `work(k)` for every `k` in `ks`, `CONCURRENT` at a time (flow_leg.jl's `run_patches`, logged here)."
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
                tsay("$label: patch $k failed: " * first(sprint(showerror, e), 300))
            end
            n = Threads.atomic_add!(done, 1) + 1
            if n % 50 == 0 || n == length(ks)
                dt = time() - t0
                tsay(@sprintf("%s: %d of %d patches in %.0f s, about %.0f min left", label, n, length(ks), dt, dt / n * (length(ks) - n) / 60))
            end
        end
    end
    return failed
end

# ---------------------------------------------------------------------------------------------------------------------
# step A: the venous flows

"flow_leg.jl's `prune`: the segments it keeps (the leaves `keep_leaf` marks and every vessel above them), in order."
function kept(t, keep_leaf)
    keep = copy(keep_leaf)
    for j in findall(keep_leaf)
        p = t.parent[j]
        while p > 0 && !keep[p]
            keep[p] = true
            p = t.parent[p]
        end
    end
    return findall(keep)
end

"flow_leg.jl's `venous_forest`, with each forest segment's venous patch (`segm`) and index in that patch's tree (`segj`)."
function forest(k, W)
    pieces, ms, ns, nm = V.Flow.VTree[], Int[], Int[], Int[]
    segm, segj = Int32[], Int32[]
    for m in W.by_k[k]
        js = last(W.per_m[m][findfirst(p -> first(p) == k, W.per_m[m])])
        T = deserialize(joinpath(CACHE_DIR, "vein_patch_$m.jls")).tree
        keep = falses(length(T.parent))
        keep[js] .= true
        push!(pieces, prune(T, keep))
        kj = kept(T, keep)
        append!(segm, fill(Int32(m), length(kj)))
        append!(segj, Int32.(kj))
        push!(ms, m)
        push!(ns, length(js))
        push!(nm, length(T.parent))
    end
    tree = V.Flow.merge_vtrees(pieces...)
    count(iszero, tree.parent) == length(ms) || error("patch $k: a venous patch tree has more than one root")
    length(segm) == length(tree.parent) || error("patch $k: $(length(segm)) segments mapped, the forest has $(length(tree.parent))")
    return (tree = tree, m = ms, n = ns, nm = nm, segm = segm, segj = segj)
end

"The plasma each forest segment carries (mL/s): every terminal's plasma flow (blood flow `qt`) up its own venous path, as the circuit routes it."
function forest_flows(tree, vleaf, qt)
    nv = length(tree.parent)
    qr = zeros(nv)
    for i in eachindex(vleaf)
        qr[vleaf[i]] += max(qt[i], 1.0f-9) * PF
    end
    order = V.Flow._bfs(tree.parent)
    for h in nv:-1:1
        s = order[h]
        p = tree.parent[s]
        p > 0 && (qr[p] += qr[s])
    end
    return qr
end

"Every venous patch segment's plasma flow (mL/s) from the beds of all the arterial patches it drains, and each venous patch's outflow."
function venous_flows(W, ks)
    path = joinpath(T_DIR, "venous_flows.jls")
    isfile(path) && return deserialize(path)
    t0 = time()
    M = length(VEIN_PATCHES.voxels)
    Q = [Float32[] for _ in 1:M]
    qroot = zeros(M)
    failed = run_each(ks, "venous flows") do k
        vf = forest(k, W)
        vleaf = deserialize(joinpath(FLOW_DIR, "vleaf_$k.jls"))
        qt = deserialize(joinpath(FLOW_DIR, "flows_$k.jls")).rest
        qr = forest_flows(vf.tree, vleaf, qt)
        lock(LOCK) do
            for (i, m) in enumerate(vf.m)
                isempty(Q[m]) && (Q[m] = zeros(Float32, vf.nm[i]))
            end
            for s in eachindex(qr)
                Q[vf.segm[s]][vf.segj[s]] += qr[s]
                vf.tree.parent[s] == 0 && (qroot[vf.segm[s]] += qr[s])
            end
        end
    end
    isempty(failed) || error("venous flows failed for patches $(failed)")
    shared = sum(count(>(0), q) for q in Q)
    tsay(@sprintf("venous flows: %d venous patches, %.1f M veins carrying flow, outflow %.3f mL/s plasma, in %.0f s", count(!isempty, Q), shared / 1.0e6, sum(qroot), time() - t0))
    return save(path, (Q = Q, qroot = qroot))
end

# ---------------------------------------------------------------------------------------------------------------------
# pass 1 at rest

"Radius `r` (mm) times the circuit's fixed lognormal factor for `key` (`radius_cv`), as `circuit` and flow_leg.jl's `resistances` vary radii."
vary(r, key) = Float32(r) * Float32(exp(SIG * V.Flow._hashnormal(key) - SIG^2 / 2))

"Plasma volume (mL) of a segment of radius `r` and length `L` (mm), as the circuit counts it (Fåhræus)."
plasma_ml(r, L) = Float64(Float32(π) * r^2 * L * 1.0f-3 * (1 - V.Flow.fahraeus(2000.0f0 * r, HD) * HD))

seglen(t, j) = Float32(sqrt(sum(abs2, t.b[:, j] .- t.a[:, j])))

"Pass 1's arteries at rest, patch `k` taking `q[k]` mL/s of blood: each segment's plasma volume and transit, and each patch's arrival at its inlet (the start of its pass-1 segment, s)."
function pass1_arteries(q)
    t = pass1.arteries.tree
    n = length(t.parent)
    Qp = zeros(n)
    for (k, j) in enumerate(PATCHES.segment)
        Qp[j] += q[k] * PF
    end
    order = PASS1.a_order
    for j in Iterators.reverse(order)
        p = t.parent[j]
        p > 0 && (Qp[p] += Qp[j])
    end
    inner = falses(n)
    foreach(p -> p > 0 && (inner[p] = true), t.parent)
    vp, tau = zeros(n), zeros(n)
    for j in 1:n
        inner[j] || continue # a pass-1 terminal is its patch's root segment, in the patch's own circuit
        vp[j] = plasma_ml(vary(PASS1_MURRAY.radius[j], j + (1 << 41)), seglen(t, j))
        tau[j] = vp[j] / max(Qp[j], 1.0e-12)
    end
    tend = zeros(n)
    for j in order
        p = t.parent[j]
        tend[j] = (p == 0 ? 0.0 : tend[p]) + tau[j]
    end
    D = [t.parent[j] == 0 ? 0.0 : tend[t.parent[j]] for j in PATCHES.segment]
    return (; Qp, vp, tau, tend, D)
end

"Pass 1's veins at rest, venous patch `m` draining `qp[m]` mL/s of plasma: each segment's plasma volume and transit, and each venous patch's delay from its root to the femoral vein's (s)."
function pass1_veins(qp)
    t = pass1.veins.tree
    n = length(t.parent)
    Qp = zeros(n)
    for (m, j) in enumerate(VEIN_PATCHES.segment)
        Qp[j] += qp[m]
    end
    order = PASS1.v_order
    for j in Iterators.reverse(order)
        p = t.parent[j]
        p > 0 && (Qp[p] += Qp[j])
    end
    inner = falses(n)
    foreach(p -> p > 0 && (inner[p] = true), t.parent)
    vp, tau = zeros(n), zeros(n)
    for j in 1:n
        (inner[j] && Qp[j] > 0) || continue
        vp[j] = plasma_ml(vary(PASS1_MURRAY_VEINS.radius[j], j + (1 << 42)), seglen(t, j))
        tau[j] = vp[j] / Qp[j]
    end
    sv = zeros(n)
    for j in order
        p = t.parent[j]
        sv[j] = (p == 0 ? 0.0 : sv[p]) + tau[j]
    end
    DV = [t.parent[j] == 0 ? 0.0 : sv[t.parent[j]] for j in VEIN_PATCHES.segment]
    return (; Qp, vp, tau, sv, DV)
end

# ---------------------------------------------------------------------------------------------------------------------
# step B: per patch

"Patch `k`'s resting circuit, exactly as flow_leg.jl's `final_patch` solves it, and its venous forest."
function rest_circuit(k, W, N, Br)
    f = LIMB_FLOW
    P = deserialize(joinpath(FLOW_DIR, "pair_$k.jls"))
    vleaf = deserialize(joinpath(FLOW_DIR, "vleaf_$k.jls"))
    art = arterial(k)
    terms = leaves(art)
    terr = fill(F32(P.ml / length(terms)), length(terms))
    vf = forest(k, W)
    vf.m == P.m || error("patch $k: its venous patches differ from the pairing's")
    dem = V.Territory.demand_field(art.b[:, terms], terr)
    kappa = F32(N.kappa[k])
    pout_rest = sum(vf.n .* Br.pv[vf.m]) / sum(vf.n)
    spec_r = V.FlowSpec(f; p_in = F32(Br.pin[k]), p_out = F32(pout_rest))
    C = V.Flow.circuit(art, terr; state = :rest, veins = vf.tree, vleaf, phys = spec_r, kappa, demand = dem)
    return (; C, vf, vleaf, terms, ml = P.ml)
end

"""
Each forest segment's plasma volume (its radius varied by its own venous patch segment, so every patch that shares a vein
sees the same one), its transit at the plasma flow `Q[m][j]` of every bed it drains, and the delay from each segment's
inlet to its venous patch's root (s).
"""
function venous_transits(vf, Q)
    t = vf.tree
    nv = length(t.parent)
    vp, tau = zeros(nv), zeros(nv)
    for s in 1:nv
        m, j = vf.segm[s], vf.segj[s]
        vp[s] = plasma_ml(vary(t.r[s], (Int64(m) << 32) | Int64(j)), seglen(t, s))
        tau[s] = vp[s] / max(Float64(Q[m][j]), 1.0e-12)
    end
    order = V.Flow._bfs(t.parent)
    sv = zeros(nv)
    for h in 1:nv
        s = order[h]
        p = t.parent[s]
        sv[s] = (p == 0 ? 0.0 : sv[p]) + tau[s]
    end
    return (; vp, tau, sv)
end

"Every bed's arrival times and exchange: at its inlet (`tin`, from the leg's root), then the circuit's own transits; `D` is the patch's arrival, `DV` each venous patch's delay to the femoral vein."
function bed_data(k, C, vf, vt, D, DV)
    return (;
        k, D = Float64(D), region = region(k),
        tin = Float32.(D .+ C.t_end[C.terminals]), ta = C.ta, tc = C.tc, tv = C.tv,
        sv = Float32.(vt.sv[C.vleaf]), dv = Float32.(DV[vf.segm[C.vleaf]]),
        fp = C.fp, ex = C.ex, kx = C.k,
    )
end

"Patch `k` at rest: its beds' arrival times (saved) and its volumes for the checks."
function patch_beds(k, W, N, Br, Q, D, DV)
    t0 = time()
    R = rest_circuit(k, W, N, Br)
    C, vf = R.C, R.vf
    F = deserialize(joinpath(FLOW_DIR, "final_$k.jls"))
    q = Float64(sum(C.resist.root_flow))
    vt = venous_transits(vf, Q)
    beds = bed_data(k, C, vf, vt, D[k], DV)
    qr = forest_flows(vf.tree, C.vleaf, C.q[C.terminals])
    share = [qr[s] / max(Float64(Q[vf.segm[s]][vf.segj[s]]), 1.0e-12) for s in eachindex(qr)]
    summ = (;
        k, region = region(k), ml = R.ml, q, dq = q / F.rest.q - 1, fp = Float64(sum(C.fp)),
        vp_art = Float64(sum(C.vp_a)), vp_beds = Float64(sum(C.fp .* (C.ta .+ C.tc .+ C.tv))),
        vp_veins = sum(share .* vt.vp), vp_veins_own = Float64(sum(C.vp_v)), max_share = maximum(share),
        ex = Float64(median(C.ex)), seconds = time() - t0,
    )
    save(joinpath(T_DIR, "beds_$k.jls"), beds)
    return save(joinpath(T_DIR, "patch_$k.jls"), summ)
end

# ---------------------------------------------------------------------------------------------------------------------
# step C: the step responses

"Arrivals binned on the grid: `a`, the plasma flow (mL/s) arriving at each time; `b`, the part that leaves the plasma, by time and back-flux rate."
struct XHist
    a::Vector{Float64}
    b::Matrix{Float64} # (NG, NK)
end
XHist() = XHist(zeros(NG), zeros(NG, NK))

"The leg's boundaries, per region: patch inlets and bed inlets (no exchange), bed outlets and venous patch roots; and the femoral vein."
struct Hists
    pin::Vector{Vector{Float64}}
    bin::Vector{Vector{Float64}}
    bout::Vector{XHist}
    vroot::Vector{XHist}
    fem::XHist
end
Hists() = Hists([zeros(NG) for _ in REGIONS], [zeros(NG) for _ in REGIONS], [XHist() for _ in REGIONS], [XHist() for _ in REGIONS], XHist())

@inline function bin_a!(h::Vector{Float64}, t, w)
    x = t / FINE
    x >= NG - 1 && return nothing # arrives after TMAX
    n = floor(Int, x)
    g = x - n
    @inbounds h[n + 1] += (1 - g) * w
    @inbounds h[n + 2] += g * w
    return nothing
end

@inline function bin_x!(H::XHist, t, a, b, y)
    x = t / FINE
    x >= NG - 1 && return nothing
    n = floor(Int, x)
    g = x - n
    m = floor(Int, y)
    u = y - m
    @inbounds begin
        H.a[n + 1] += (1 - g) * a
        H.a[n + 2] += g * a
        H.b[n + 1, m + 1] += (1 - g) * (1 - u) * b
        H.b[n + 1, m + 2] += (1 - g) * u * b
        H.b[n + 2, m + 1] += g * (1 - u) * b
        H.b[n + 2, m + 2] += g * u * b
    end
    return nothing
end

"Bins one patch's beds into `H`; returns its plasma flow and Σ flow × arrival at the femoral vein (all beds, also those after TMAX)."
function bin_patch!(H::Hists, b)
    r = b.region
    w, wt = 0.0, 0.0
    for i in eachindex(b.fp)
        fp = Float64(b.fp[i])
        bx = fp * Float64(b.ex[i])
        y = clamp((log(Float64(b.kx[i])) - LOGK[1]) / DLOGK, 0.0, NK - 1 - 1.0e-9)
        tin = Float64(b.tin[i])
        bin_a!(H.bin[r], tin, fp)
        t4 = tin + b.ta[i] + b.tv[i] + b.tc[i]
        bin_x!(H.bout[r], t4, fp, bx, y)
        t5 = t4 + b.sv[i]
        bin_x!(H.vroot[r], t5, fp, bx, y)
        t6 = t5 + b.dv[i]
        bin_x!(H.fem, t6, fp, bx, y)
        w += fp
        wt += fp * t6
    end
    bin_a!(H.pin[r], b.D, w)
    return w, wt
end

"Step response (mL/s of plasma through the boundary for a unit plasma concentration from t = 0)."
step_a(h::Vector{Float64}) = cumsum(h)
function step_x(H::XHist)
    s = cumsum(H.a)
    for m in 1:NK
        e = exp(-exp(LOGK[m]) * FINE)
        R = 0.0
        @inbounds for n in 1:NG
            R = R * e + H.b[n, m]
            s[n] -= R
        end
    end
    return s
end

"The cumulative input A(t) = ∫ c (mg·s/mL), by the trapezoid rule on the grid, as the package's InputCurve."
cumin(c) = [0.0; cumsum(0.5 * FINE .* (c[1:(end - 1)] .+ c[2:end]))]

"""
The cumulative iodine (mg) through a boundary with step response `s` for the cumulative input `A`, both on the grid:
the arrivals (the step response's increments) times A at the time since each, so a bed arriving between two grid points
sees A interpolated linearly, as the package's `_lin` reads it.
"""
function cumflux(A, s)
    ds = diff([0.0; s])
    G = zeros(NG)
    Threads.@threads for n in 2:NG
        acc = 0.0
        @inbounds for j in 1:n
            acc += ds[j] * A[n - j + 1]
        end
        G[n] = acc
    end
    return G
end

"The whole leg's step responses, from every patch's beds."
function leg_hists(ks)
    H = Hists()
    tot = Threads.Atomic{Float64}(0.0)
    tott = Threads.Atomic{Float64}(0.0)
    t0 = time()
    failed = run_each(ks, "binning") do k
        b = deserialize(joinpath(T_DIR, "beds_$k.jls"))
        lock(LOCK) do
            w, wt = bin_patch!(H, b)
            Threads.atomic_add!(tot, w)
            Threads.atomic_add!(tott, wt)
        end
    end
    isempty(failed) || error("binning failed for patches $(failed)")
    tsay(@sprintf("binned %d patches in %.0f s", length(ks), time() - t0))
    return H, (fp = tot[], mean_arrival = tott[] / tot[])
end

# ---------------------------------------------------------------------------------------------------------------------
# step D: the body

"transport.jl's `_interp`: linear interpolation of samples y at times ts, held flat outside."
interp(t, ts, y) = t <= ts[1] ? y[1] : t >= ts[end] ? y[end] :
    (k = searchsortedlast(ts, t); w = (t - ts[k]) / (ts[k + 1] - ts[k]); y[k] + w * (y[k + 1] - y[k]))

"The whole-body model at rest with the right leg's share of the cardiac output at the grown leg's flow `q_leg` (mL/s); the other regions keep their proportions."
function leg_body(q_leg)
    kw = (; co = REST.co, beat = REST.beat, hct = HD, f_cor = REST.f_cor, coronary = V.LUNG_TRANSPORT.coronary)
    legshare(BC) = sum(BC.q[i] for i in eachindex(BC.segs) if BC.segs[i].name in LEG_BEDS)
    base = V.Body.human_circulation(; kw...)
    f0 = legshare(base)
    s = 0.15 * (q_leg / REST.co) / f0 # the default right leg is 0.15 of the peripheral tissues
    o = (1 - s) / (1 - 0.15)
    regions = (
        (:arm_r, :antecubital_r, F32(0.05o)), (:arm_l, :antecubital_l, F32(0.05o)), (:leg_r, :femoral_r, F32(s)),
        (:leg_l, :femoral_l, F32(0.15o)), (:upper, :svc_top, F32(0.3o)), (:lower, :ivc_low, F32(0.3o)),
    )
    BC = V.Body.human_circulation(; kw..., regions)
    f = legshare(BC)
    abs(f * REST.co / q_leg - 1) < 1.0e-3 || error("the body's right leg takes $(f * REST.co) mL/s, the grown leg $(q_leg)")
    return (; BC, share = Float64(f), lumped = Float64(f0 * REST.co))
end

simulate_body(B, det = Dict{Symbol, Any}()) = V.Body.simulate(B.BC, INJECTION; tmax = TMAX, dt = 0.002f0, every = 0.1f0, detailed = det)

"The leg's arterial input: the plasma concentration (mg/mL) leaving the abdominal aorta, on the grid."
aif(r) = [Float64(interp(t, r.t, r.out[:desc_aorta])) / PF for t in GRID]

"The leg's femoral venous plasma concentration (mg/mL) every DT_OUT, for the input `c`: the cumulative outflow's differences, as `transport` takes them."
function femoral_outflow(c, sfem, fp_leg)
    G = cumflux(cumin(c), sfem)[1:STRIDE:end]
    return [k == 1 ? 0.0 : (G[k] - G[k - 1]) / DT_OUT / fp_leg for k in eachindex(G)]
end

"`transport`'s loop: the leg's venous outflow fed back for the body's right-leg beds until it changes by at most `tol`."
function body_loop(B, sfem, fp_leg; passes = 10, tol = 1.0e-3)
    r = simulate_body(B)
    gap0 = Float64(r.injected - r.held - r.excreted) # the body's own bookkeeping, its lumped right leg in place
    c = aif(r)
    cs = femoral_outflow(c, sfem, fp_leg)
    npass, change = 0, 0.0
    while npass < passes
        fb = cs
        back = t -> PF * interp(Float32(t), TC, fb)
        r = simulate_body(B, Dict{Symbol, Any}(b => back for b in LEG_BEDS))
        c = aif(r)
        cs = femoral_outflow(c, sfem, fp_leg)
        npass += 1
        change = B.share * maximum(abs.(cs .- fb)) / max(maximum(abs.(cs)), floatmin(Float64))
        tsay(@sprintf("body pass %d: weighted change of the leg's outflow %.2e", npass, change))
        change <= tol && break
    end
    return (; r, c, cs, passes = npass, change, gap0)
end

"The index of time `tq` (s) in the times `ts`."
at_time(ts, tq) = findfirst(t -> abs(t - tq) < 1.0e-3, ts)

"First time a curve reaches `frac` of its peak, its peak and when (s)."
function timing(t, y; frac = 0.1)
    i = argmax(y)
    j = findfirst(>=(frac * y[i]), y)
    return (arrive = t[j], peak = y[i], tpeak = t[i])
end

function report(B, L, H, stats, PA, PV, S, ks, cmp)
    lines = String[]
    out(s) = (push!(lines, s); tsay(s))
    c = L.c
    fp = stats.fp
    # cumulative fluxes (mg) through every boundary
    A = cumin(c)
    G1 = fp .* A
    nr = length(REGIONS)
    G2 = [cumflux(A, step_a(H.pin[r])) for r in 1:nr]
    G3 = [cumflux(A, step_a(H.bin[r])) for r in 1:nr]
    G4 = [cumflux(A, step_x(H.bout[r])) for r in 1:nr]
    G5 = [cumflux(A, step_x(H.vroot[r])) for r in 1:nr]
    G6 = cumflux(A, step_x(H.fem))
    fr = 1:round(Int, 1 / FINE):NG # every second
    tf = GRID[fr]
    art1 = (G1 .- sum(G2))[fr]
    art2 = [(G2[r] .- G3[r])[fr] for r in 1:nr]
    beds = [(G3[r] .- G4[r])[fr] for r in 1:nr]
    ven2 = [(G4[r] .- G5[r])[fr] for r in 1:nr]
    ven1 = (sum(G5) .- G6)[fr]
    total = (G1 .- G6)[fr]
    ml = [sum(s.ml for s in S if s.region == r) for r in 1:nr]

    out(@sprintf(
        "body: the right leg takes %.0f mL/min of the cardiac output's %.0f (%.2f%%; the body's own right leg %.0f); %d feedback passes, last change %.1e",
        60 * B.share * REST.co, 60 * REST.co, 100 * B.share, 60 * B.lumped, L.passes, L.change,
    ))
    a = timing(GRID, PF .* c)
    out(@sprintf("arterial input (the abdominal aorta's outflow, blood): arrives (10%% of peak) at %.1f s, peaks at %.2f mg/mL at %.1f s", a.arrive, a.peak, a.tpeak))
    v = timing(TC, PF .* L.cs)
    out(@sprintf(
        "femoral venous outflow (blood): arrives at %.1f s, peaks at %.3f mg/mL at %.1f s (%.0f%% of the arterial peak); at 60 s %.3f, at 120 s %.3f",
        v.arrive, v.peak, v.tpeak, 100 * v.peak / a.peak, PF * L.cs[at_time(TC, 60)], PF * L.cs[end],
    ))
    vp = (
        art1 = sum(PA.vp), art2 = sum(s.vp_art for s in S), beds = sum(s.vp_beds for s in S), ven2 = sum(s.vp_veins for s in S),
        ven1 = sum(PV.vp), own = sum(s.vp_veins_own for s in S),
    )
    vtot = vp.art1 + vp.art2 + vp.beds + vp.ven2 + vp.ven1
    out(@sprintf(
        "plasma in the leg: pass-1 arteries %.1f mL, patch arteries %.1f, beds %.1f, patch veins %.1f (each patch's own circuit counts %.1f: the shared veins once per patch), pass-1 veins %.1f; total %.1f",
        vp.art1, vp.art2, vp.beds, vp.ven2, vp.own, vp.ven1, vtot,
    ))
    out(@sprintf(
        "mean transit at rest: plasma volume ÷ plasma flow %.1f s; the beds' flow-weighted arrival at the femoral vein %.1f s (the same by the central volume principle)",
        vtot / fp, stats.mean_arrival,
    ))
    injected, held, excreted = Float64(L.r.injected), Float64(L.r.held), Float64(L.r.excreted)
    out(@sprintf(
        "iodine at %.0f s: injected %.0f mg; in the grown leg %.0f mg; the body outside it holds %.0f and has excreted %.0f, leaving %.0f for the leg; with its own lumped leg the body leaves %.0f unaccounted",
        TMAX, injected, total[end], held, excreted, injected - held - excreted, L.gap0,
    ))
    for tq in (15, 20, 30, 60, 90, 120)
        i = at_time(tf, tq)
        out(@sprintf(
            "  %3d s: in the leg %6.1f mg: pass-1 arteries %5.1f, patch arteries %5.1f, beds %6.1f, patch veins %5.1f, pass-1 veins %5.1f",
            tq, total[i], art1[i], sum(r -> art2[r][i], 1:nr), sum(r -> beds[r][i], 1:nr), sum(r -> ven2[r][i], 1:nr), ven1[i],
        ))
    end
    for r in 1:nr
        conc = (art2[r] .+ beds[r] .+ ven2[r]) ./ ml[r]
        tb = timing(tf, conc)
        at(tq) = conc[at_time(tf, tq)]
        out(@sprintf(
            "%-5s %6.0f mL: the patches' iodine per mL of tissue peaks at %.3f mg/mL at %.0f s (beds %.0f%% of it); at 30 s %.3f, 60 s %.3f, 120 s %.3f",
            REGIONS[r], ml[r], tb.peak, tb.tpeak, 100 * beds[r][argmax(conc)] / (conc[argmax(conc)] * ml[r]), at(30), at(60), at(120),
        ))
    end
    dq = [s.dq for s in S]
    secs = [s.seconds for s in S]
    out(@sprintf(
        "patches: %d; resting flow reproduced to %.1e (max); first-pass extraction median %.2f; the most a patch takes of a shared vein's flow, median %.2f; %.1f s each",
        length(S), maximum(abs, dq), median([s.ex for s in S]), median([s.max_share for s in S]), median(secs),
    ))
    out(@sprintf("arrival at the patch inlets: median %.1f s after the leg's root (max %.1f); from the venous patches' roots to the femoral vein: median %.1f s (max %.1f)", median(PA.D[ks]), maximum(PA.D[ks]), median(filter(>(0), PV.DV)), maximum(PV.DV)))
    save(joinpath(T_DIR, "result.jls"), (;
        grid = GRID, c, tc = TC, cs = L.cs, frames = tf, art1, art2, beds, ven2, ven1, total, ml, regions = REGIONS,
        steps = (pin = step_a.(H.pin), bin = step_a.(H.bin), bout = step_x.(H.bout), vroot = step_x.(H.vroot), fem = step_x(H.fem)),
        fp, body = (share = B.share, lumped = B.lumped, co = REST.co, passes = L.passes, change = L.change, injected, held, excreted, gap0 = L.gap0, aorta = L.r.out[:desc_aorta], t = L.r.t),
        vp, mean_arrival = stats.mean_arrival, D = PA.D, DV = PV.DV, patches = S, compare = cmp,
    ))
    open(joinpath(T_DIR, "summary.txt"), "w") do io
        foreach(l -> println(io, l), lines)
    end
    return nothing
end

# ---------------------------------------------------------------------------------------------------------------------

function main()
    W = venule_map()
    Br = equal_share(:rest)
    N = deserialize(joinpath(FLOW_DIR, "network.jls"))
    ks = [k for k in eachindex(PATCHES.voxels) if isfile(joinpath(FLOW_DIR, "final_$k.jls"))]
    q = zeros(length(PATCHES.voxels))
    for k in ks
        q[k] = deserialize(joinpath(FLOW_DIR, "final_$k.jls")).rest.q
    end
    tsay(@sprintf("%d patches, resting flow %.1f mL/min", length(ks), 60 * sum(q)))
    VQ = venous_flows(W, ks)
    PA = pass1_arteries(q)
    PV = pass1_veins(VQ.qroot)
    todo = filter(k -> !isfile(joinpath(T_DIR, "beds_$k.jls")), ks)
    failed = run_each(k -> patch_beds(k, W, N, Br, VQ.Q, PA.D, PV.DV), todo, "patches")
    isempty(failed) || error("patches failed: $(failed)")
    S = [deserialize(joinpath(T_DIR, "patch_$k.jls")) for k in ks]
    H, stats = leg_hists(ks)
    B = leg_body(sum(q))
    L = body_loop(B, step_x(H.fem), stats.fp)
    cmp = compare_three(W, N, Br, L.c, ks)
    report(B, L, H, stats, PA, PV, S, ks, cmp)
    return tsay(@sprintf("done in %.1f h", (time() - STARTED) / 3600))
end

rel(x, y) = maximum(abs.(x .- y)) / max(maximum(abs.(y)), 1.0e-30)

"""
Patch `k` alone, from its own inlet, for the plasma input `c`: its compartments from the step responses against the
package's `frame_masses` and `sinus_outflow` on the same circuit (its veins at its own flows, as if no other patch shared
them), each as the largest difference relative to its peak; and the package's own compartments against its inflow less
its outflow (relative to the beds' peak).
"""
function compare_patch(k, W, N, Br, c)
    inp = V.Flow.InputCurve(Float32.(c), FINE)
    frames = collect(0.0f0:1.0f0:TMAX)
    fr = 1:round(Int, 1 / FINE):NG
    t = @elapsed R = rest_circuit(k, W, N, Br)
    C, vf = R.C, R.vf
    F = deserialize(joinpath(FLOW_DIR, "final_$k.jls"))
    ok = C.terminals == R.terms && C.vleaf == R.vleaf
    Qloc = [Float32[] for _ in 1:length(VEIN_PATCHES.voxels)]
    for (i, m) in enumerate(vf.m)
        Qloc[m] = zeros(Float32, vf.nm[i])
    end
    qr = forest_flows(vf.tree, C.vleaf, C.q[C.terminals])
    for s in eachindex(qr)
        Qloc[vf.segm[s]][vf.segj[s]] += qr[s]
    end
    vt = venous_transits(vf, Qloc)
    off2 = C.t_end[C.terminals] .+ C.ta .+ C.tv .+ Float32.(vt.sv[C.vleaf])
    swap = Dict(:tau_v => Float32.(vt.tau), :s_v => Float32.(vt.sv), :vp_v => Float32.(vt.vp), :offset => off2)
    C2 = V.Flow.Circuit((get(swap, f, getfield(C, f)) for f in fieldnames(V.Flow.Circuit))...)
    am, vm, tm = V.Flow.frame_masses(C2, inp, frames)
    cum = V.Flow.sinus_outflow(C2, inp, TC)
    cum_own = V.Flow.sinus_outflow(C, inp, TC)
    H = Hists()
    b = bed_data(k, C, vf, vt, 0.0, zeros(length(VEIN_PATCHES.voxels)))
    bin_patch!(H, b)
    rg = b.region
    A = cumin(c)
    G2 = cumflux(A, step_a(H.pin[rg]))
    G3 = cumflux(A, step_a(H.bin[rg]))
    G4 = cumflux(A, step_x(H.bout[rg]))
    G5 = cumflux(A, step_x(H.vroot[rg]))
    art, tis, ven = vec(sum(am; dims = 2)), vec(sum(tm; dims = 2)), vec(sum(vm; dims = 2))
    inflow = Float64(sum(C.fp)) .* A[fr]
    return (;
        k, region = REGIONS[rg], height = HEIGHT[k], beds = length(C.terminals), ml = R.ml, seconds = t, ok,
        dq = sum(C.resist.root_flow) / F.rest.q - 1,
        arteries = rel((G2 .- G3)[fr], art), tissue = rel((G3 .- G4)[fr], tis), veins = rel((G4 .- G5)[fr], ven),
        outflow = rel(G5[1:STRIDE:end], cum), key = rel(cum, cum_own),
        package_balance = maximum(abs.(art .+ tis .+ ven .- (inflow .- cum[1:10:end]))) / maximum(tis),
    )
end

"Three patches (the foot, the calf, the thigh) against the package, for the plasma input `c`."
function compare_three(W, N, Br, c, ks)
    near(h) = argmin(k -> PATCHES.voxels[k] > 0 ? abs(HEIGHT[k] - h) : Inf, eachindex(PATCHES.voxels))
    foot = first(k for k in ks if region(k) == 1)
    out = [compare_patch(k, W, N, Br, c) for k in unique([foot, near(0.25), near(0.8)])]
    for x in out
        tsay(@sprintf(
            "patch %4d (%s, height %.2f, %d beds): vs the package, arteries %.1e, beds %.1e, veins %.1e, venous outflow %.1e; the package's own compartments vs its inflow less outflow %.1e; the volume key moves the outflow %.1e; flow %.1e; orders %s",
            x.k, x.region, x.height, x.beds, x.arteries, x.tissue, x.veins, x.outflow, x.package_balance, x.key, x.dq, x.ok,
        ))
    end
    return out
end

"A few patches alone, from their own inlets, against the package's frame_masses and sinus_outflow on the same circuit."
function test()
    W = venule_map()
    Br = equal_share(:rest)
    N = deserialize(joinpath(FLOW_DIR, "network.jls"))
    ks = [k for k in eachindex(PATCHES.voxels) if isfile(joinpath(FLOW_DIR, "final_$k.jls"))]
    q_leg = sum(deserialize(joinpath(FLOW_DIR, "final_$k.jls")).rest.q for k in ks)
    B = leg_body(q_leg)
    r = simulate_body(B)
    c = aif(r)
    a = timing(GRID, PF .* c)
    tsay(@sprintf(
        "body: right leg %.0f mL/min (%.2f%% of %.0f; the body's own %.0f); input arrives %.1f s, peaks %.2f mg/mL blood at %.1f s; injected %.0f mg, held %.0f, excreted %.0f",
        60 * q_leg, 100 * B.share, 60 * REST.co, 60 * B.lumped, a.arrive, a.peak, a.tpeak, r.injected, r.held, r.excreted,
    ))
    compare_three(W, N, Br, c, ks)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    get(ARGS, 1, "") == "test" ? test() : main()
end
