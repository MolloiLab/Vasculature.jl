# The leg's blood flow (limb_pipeline.jl, section 7, step 2): V.flow's circuit on every
# arterial patch with the venules in its own tissue, the patches joined through pass 1.
#
#     julia -t auto flow_leg.jl            the whole leg (steps 0-4); steps already saved are skipped
#     julia -t auto flow_leg.jl test       steps 0 and 1 on a few patches, one at a time and then several at once
#
# Each tree has 342 million segments, too many for one circuit on a workstation, so every arterial patch is solved on
# its own:
#   0. Venules. Every venous terminal goes to the arterial patch whose tissue it ends in, so each patch's beds drain into
#      the venules of its own tissue, about one each, and no venule serves two patches. A patch's venous trees are the
#      venous patches' trees cut down to those venules and the veins above them. Saved as flow/venules.jls.
#   1. Pairing, per patch. `couple_veins` pairs the beds with those venules at maximal dilation, at boundary pressures
#      from pass 1 at the equal share, with the patch's own κ set so that it takes its maximal-flow share. One more solve
#      at 1.25 κ gives its resistance at any κ as a + b κ (κ scales only the microvessels). Saved as flow/pair_<k>.jls
#      and flow/vleaf_<k>.jls.
#   2. The network. Pass 1, arteries and veins at their Murray radii, with every patch as its resistance: each patch's
#      inlet and outlet pressure, and κ for the leg above the ankle and for the foot, so that each takes its maximal flow
#      (`r_total_dilated`, the foot at `FOOT_MAX_SHARE`; a patch across the ankle takes a κ in between). Saved as
#      flow/network.jls.
#   3. Final, per patch: maximal dilation at those pressures and that κ, then rest with the same κ and its tone constant
#      set by its own O2 balance, at pass 1's resting pressures. Saved as flow/final_<k>.jls (summaries) and
#      flow/flows_<k>.jls (every terminal's flow, both states, and its venous SO2 at rest).
#   4. Checks: section 7's (resting inflow, femoral venous SO2, VO2, reversed beds, tone), the flow by height, and how
#      far the final flows move pass 1's pressures. Printed and saved as flow/summary.txt.
#
# Numbers only. Everything is made from XCAT, so it stays in LIMB_CACHE_DIR (lab storage). FLOW_CONCURRENT sets how many
# patches are solved at once (default 8).

using Printf, Serialization, Statistics

for name in ("XCAT_LEG_NRB", "LIMB_CACHE_DIR")
    haskey(ENV, name) || error("set $name first: without it the notebook falls back to another machine's paths")
end

"Windows may run this process on the efficiency cores alone, about 6× slower (2026-11-03): opt it out."
function full_speed()
    Sys.iswindows() || return true
    state = UInt32[1, 1, 0] # PROCESS_POWER_THROTTLING_STATE: version 1, execution speed controlled, not throttled
    process = ccall((:GetCurrentProcess, "kernel32"), Ptr{Cvoid}, ())
    return ccall((:SetProcessInformation, "kernel32"), Cint, (Ptr{Cvoid}, Cint, Ptr{UInt32}, UInt32), process, 4, state, 12) != 0
end
full_speed() || @warn "Windows may keep this process on its efficiency cores, about 6× slower"

const STARTED = time()
include(get(ENV, "LIMB_NOTEBOOK", joinpath(@__DIR__, "limb_pipeline.jl")))

const F32 = Float32
const FLOW_DIR = mkpath(joinpath(CACHE_DIR, "flow"))
const CONCURRENT = parse(Int, get(ENV, "FLOW_CONCURRENT", "8"))
const LOG = open(joinpath(FLOW_DIR, "progress.log"), "a")
const LOCK = ReentrantLock()
function say(s)
    lock(LOCK) do
        line = @sprintf("[%6.0f s] %s", time() - STARTED, s)
        println(line)
        flush(stdout)
        println(LOG, line)
        flush(LOG)
    end
    return nothing
end
say(@sprintf("notebook loaded, %d threads, %d patches at once", Threads.nthreads(), CONCURRENT))

"Write `x` to `path` through a .part file, so an interrupted run leaves no broken file."
function save(path, x)
    serialize(path * ".part", x)
    mv(path * ".part", path; force = true)
    return x
end

# ---------------------------------------------------------------------------------------------------------------------
# trees

leaves(t) = (has = falses(length(t.parent)); foreach(p -> p > 0 && (has[p] = true), t.parent); findall(!, has))

"Tree `t` cut down to the terminals `keep_leaf` marks and the vessels above them."
function prune(t, keep_leaf)
    keep = copy(keep_leaf)
    for j in findall(keep_leaf)
        p = t.parent[j]
        while p > 0 && !keep[p]
            keep[p] = true
            p = t.parent[p]
        end
    end
    ks = findall(keep)
    idx = zeros(Int32, length(t.parent))
    idx[ks] = 1:length(ks)
    parent = Int32[t.parent[j] == 0 ? 0 : idx[t.parent[j]] for j in ks]
    return V.Flow.VTree(parent, t.a[:, ks], t.b[:, ks], t.r[ks], fill(false, length(ks)))
end

"Arterial patch `k`'s tree as a flow tree (no segment imaged)."
function arterial(k)
    A = deserialize(joinpath(CACHE_DIR, "patch_$k.jls")).tree
    return V.Flow.VTree(A.parent, A.a, A.b, A.r, fill(false, length(A.r)))
end

"Tissue (mL) of arterial patch `k`."
patch_ml(k) = PATCHES.voxels[k] * LEG.grid.h^3 / 1000

"The share of arterial patch `k`'s tissue below the ankle."
below(k) = (1 - PATCH_MAX_SHARE[k]) / (1 - FOOT_MAX_SHARE)

"Each venous patch's maximal-flow share per terminal, as `PATCH_MAX_SHARE` for the arterial ones."
const VEIN_MAX_SHARE = let label = VEIN_PATCHES.label, n = zeros(Int, length(VEIN_PATCHES.voxels))
    for k in 1:(ANKLE.slice - 1), j in axes(label, 2), i in axes(label, 1)
        l = label[i, j, k]
        l > 0 && (n[l] += 1)
    end
    1 .- (1 - FOOT_MAX_SHARE) .* n ./ max.(VEIN_PATCHES.voxels, 1)
end

"Each arterial patch's height in the leg (0 at its lowest pass-1 terminal, 1 at its highest)."
const HEIGHT = let z = pass1.arteries.tree.b[3, PATCHES.segment]
    (z .- minimum(z)) ./ (maximum(z) - minimum(z))
end

# ---------------------------------------------------------------------------------------------------------------------
# step 0: the venules of each arterial patch's tissue

"For each venous patch, its terminals by the arterial patch whose tissue they end in; for each arterial patch, the venous patches holding its venules."
function venule_map()
    path = joinpath(FLOW_DIR, "venules.jls")
    isfile(path) && return deserialize(path)
    t0 = time()
    K, M = length(PATCHES.voxels), length(VEIN_PATCHES.voxels)
    label, g = PATCHES.label, LEG.grid
    per_m = Vector{Vector{Pair{Int32, Vector{Int32}}}}(undef, M)
    Threads.@threads :dynamic for m in 1:M
        T = deserialize(joinpath(CACHE_DIR, "vein_patch_$m.jls")).tree
        groups = Dict{Int32, Vector{Int32}}()
        for j in leaves(T)
            c = ntuple(d -> round(Int, (T.b[d, j] - g.origin[d]) / g.h) + 1, 3)
            all(1 .<= c .<= size(label)) || continue
            k = label[c...]
            k > 0 && push!(get!(groups, Int32(k), Int32[]), Int32(j))
        end
        per_m[m] = sort!(collect(groups); by = first)
    end
    by_k = [Int32[] for _ in 1:K]
    n = zeros(Int, K)
    for m in 1:M, (k, js) in per_m[m]
        push!(by_k[k], m)
        n[k] += length(js)
    end
    W = save(path, (per_m = per_m, by_k = by_k))
    used = findall(>(0), PATCHES.voxels)
    r = n[used] ./ max.(PATCH_TERMINALS[used], 1)
    say(@sprintf(
        "venules: %d venous terminals placed in %.0f s; per arterial patch, venules ÷ its terminals: median %.2f (5%%: %.2f, min %.2f, max %.2f); patches without any: %d",
        sum(n), time() - t0, median(r), quantile(r, 0.05), minimum(r), maximum(r), count(==(0), n[used]),
    ))
    return W
end

"Arterial patch `k`'s venous trees: the venous patches' trees cut down to the venules of its tissue, joined (one root per venous patch, in the order of `m`)."
function venous_forest(k, W)
    pieces, ms, ns = V.Flow.VTree[], Int[], Int[]
    for m in W.by_k[k]
        js = last(W.per_m[m][findfirst(p -> first(p) == k, W.per_m[m])])
        T = deserialize(joinpath(CACHE_DIR, "vein_patch_$m.jls")).tree
        keep = falses(length(T.parent))
        keep[js] .= true
        push!(pieces, prune(T, keep))
        push!(ms, m)
        push!(ns, length(js))
    end
    tree = V.Flow.merge_vtrees(pieces...)
    count(iszero, tree.parent) == length(ms) || error("patch $k: a venous patch tree has more than one root")
    return (tree = tree, m = ms, n = ns)
end

# ---------------------------------------------------------------------------------------------------------------------
# pass 1

"Each segment's resistance (mmHg·s/mL) at radii `r`: Poiseuille with the in vivo viscosity, each radius times circuit's fixed lognormal factor (`radius_cv`, salted by `salt`) and `widen(j, r)`; a terminal segment (inside its patch's tree) is 0."
function resistances(t, r, salt; widen = (j, r) -> 1.0f0, phys = LIMB_FLOW)
    sig = sqrt(log(1 + phys.radius_cv^2))
    has = falses(length(t.parent))
    foreach(p -> p > 0 && (has[p] = true), t.parent)
    R = zeros(length(t.parent))
    for j in eachindex(t.parent)
        has[j] || continue
        rj = F32(r[j]) * F32(exp(sig * V.Flow._hashnormal(j + salt) - sig^2 / 2))
        rj *= widen(j, rj)
        L = F32(sqrt(sum(abs2, t.b[:, j] .- t.a[:, j])))
        R[j] = V.Flow.poiseuille(L, rj, phys.plasma_viscosity_mPas * V.Flow.mu_rel(2000 * rj, phys.hd))
    end
    return R
end

"Pressures (mmHg) at the end of every segment of a pass-1 tree whose terminals `ends` carry `q` (mL/s) out of it (arteries, `sign` -1) or into it (veins, +1), from `p0` at the root."
function pass_pressures(t, order, R, ends, q, p0, sign)
    Q = zeros(length(t.parent))
    for (j, x) in zip(ends, q)
        Q[j] += x
    end
    for j in Iterators.reverse(order)
        p = t.parent[j]
        p > 0 && (Q[p] += Q[j])
    end
    P = zeros(length(t.parent))
    for j in order
        p = t.parent[j]
        P[j] = (p == 0 ? p0 : P[p]) + sign * R[j] * Q[j]
    end
    return P
end

"Pass 1 as the network joins it: orders and resistances at maximal dilation (the imaged arteries wider by `epicardial_dilation`, as circuit widens them) and at rest."
const PASS1 = let f = LIMB_FLOW, A = pass1.arteries.tree, W = pass1.veins.tree
    widen(j, r) = A.isbase[j] && 2000 * r >= f.d_active ? 1 + f.epicardial_dilation : 1.0f0
    (
        a_order = first(downstream(A.parent, j -> 1.0)), v_order = first(downstream(W.parent, j -> 1.0)),
        ra_stress = resistances(A, PASS1_MURRAY.radius, 1 << 41; widen), ra_rest = resistances(A, PASS1_MURRAY.radius, 1 << 41),
        rv = resistances(W, PASS1_MURRAY_VEINS.radius, 1 << 42),
    )
end

"Pass 1's pressures at the equal share (each patch its terminal count, at maximal flow times its maximal-flow share): every arterial patch's inlet and every venous patch's outlet (mmHg)."
function equal_share(state)
    f = LIMB_FLOW
    stress = state === :stress
    content = 1.34 * f.hb / 100 * f.sao2
    q_leg = stress ? Q_MAX / 60 : f.mvo2 / (content * f.oef) * f.rho * SOFT_ML / 60 # mL/s
    wa = [max(PATCH_TERMINALS[k], 1) * (stress ? PATCH_MAX_SHARE[k] : 1.0) for k in eachindex(PATCH_TERMINALS)]
    wv = [max(VEIN_PATCH_TERMINALS[m], 1) * (stress ? VEIN_MAX_SHARE[m] : 1.0) for m in eachindex(VEIN_PATCH_TERMINALS)]
    A, W = pass1.arteries.tree, pass1.veins.tree
    Pa = pass_pressures(A, PASS1.a_order, stress ? PASS1.ra_stress : PASS1.ra_rest, PATCHES.segment, q_leg .* wa ./ sum(wa), Float64(f.p_in), -1)
    Pv = pass_pressures(W, PASS1.v_order, PASS1.rv, VEIN_PATCHES.segment, q_leg .* wv ./ sum(wv), Float64(f.p_out), +1)
    return (pin = Pa[PATCHES.segment], pv = Pv[VEIN_PATCHES.segment])
end

# ---------------------------------------------------------------------------------------------------------------------
# step 1: pairing

reversed(C) = count(i -> C.pa[C.terminals[i]] <= C.pv[C.vleaf[i]], eachindex(C.terminals))

"Pair patch `k`'s beds with its venules and measure its resistance at two κ (not saved). `phys` replaces `LIMB_FLOW` (to test other values)."
function pair_core(k, W, B; phys = LIMB_FLOW)
    f = phys
    t0 = time()
    art = arterial(k)
    terms = leaves(art)
    ml = patch_ml(k)
    terr = fill(F32(ml / length(terms)), length(terms))
    vf = venous_forest(k, W)
    pin = B.pin[k]
    pout = sum(vf.n .* B.pv[vf.m]) / sum(vf.n)
    target = f.p_in / f.r_total_dilated * PATCH_MAX_SHARE[k] # mL/min/g
    spec = V.FlowSpec(f; p_in = F32(pin), p_out = F32(pout), r_total_dilated = F32(pin / target))
    vleaf, _, kstar = V.couple_veins(art, terr, vf.tree, terms, leaves(vf.tree); phys = spec, assemble = false)
    C = V.Flow.circuit(art, terr; state = :stress, veins = vf.tree, vleaf, phys = spec, kappa = 1.25f0 * kstar)
    q1 = target * f.rho * ml / 60 # mL/s at κ*, the κ search's target
    q2 = sum(C.resist.root_flow) # mL/s at 1.25 κ*
    r1, r2 = (pin - pout) / q1, (pin - pout) / q2
    b = (r2 - r1) / (0.25 * kstar)
    qroot = Float64.(C.qv[findall(iszero, vf.tree.parent)])
    qt = C.q[C.terminals]
    return (
        summary = (;
            k, nterm = length(terms), ml, kstar = Float64(kstar), a = r1 - b * kstar, b, pin, pout, m = vf.m, n = vf.n,
            w = qroot ./ sum(qroot), outflow = sum(qroot) / q2, reversed = reversed(C), flow_cv = std(qt) / mean(qt),
            seconds = time() - t0,
        ),
        vleaf = vleaf,
    )
end

function pair_patch(k, W, B)
    r = pair_core(k, W, B)
    save(joinpath(FLOW_DIR, "vleaf_$k.jls"), r.vleaf)
    save(joinpath(FLOW_DIR, "pair_$k.jls"), r.summary)
    return r.summary
end

"Run `work(k)` for every `k` in `ks`, `CONCURRENT` at a time, with progress every 25 patches; failures are listed, not fatal."
function run_patches(work, ks, label)
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
                say("$label: patch $k failed: " * first(sprint(showerror, e), 300))
            end
            n = Threads.atomic_add!(done, 1) + 1
            if n % 25 == 0 || n == length(ks)
                dt = time() - t0
                say(@sprintf("%s: %d of %d patches in %.0f s (%.1f s each), about %.0f min left", label, n, length(ks), dt, dt / n, dt / n * (length(ks) - n) / 60))
            end
        end
    end
    return failed
end

# ---------------------------------------------------------------------------------------------------------------------
# step 2: the network

"""
Pass 1 with every patch as its resistance a + b κ: the inlet and outlet pressure of each patch and κ above and below the
ankle, so that each takes its maximal flow. `ra` replaces pass 1's arterial resistances at maximal dilation; with
`path = nothing` nothing is read or saved and nothing printed (to test other values).
"""
function network(P; ra = PASS1.ra_stress, path = joinpath(FLOW_DIR, "network.jls"))
    quiet = path === nothing
    !quiet && isfile(path) && return deserialize(path)
    f = LIMB_FLOW
    K = length(P)
    ok = [p !== nothing for p in P]
    fb = [below(k) for k in 1:K]
    mass = [patch_ml(k) * f.rho for k in 1:K]
    per_g = f.p_in / f.r_total_dilated # mL/min/g at maximal dilation
    T_up = per_g * sum(mass[k] * (1 - fb[k]) for k in 1:K if ok[k]) / 60
    T_dn = per_g * FOOT_MAX_SHARE * sum(mass[k] * fb[k] for k in 1:K if ok[k]) / 60
    A, W = pass1.arteries.tree, pass1.veins.tree
    pin = [ok[k] ? P[k].pin : 0.0 for k in 1:K]
    pout = [ok[k] ? P[k].pout : 0.0 for k in 1:K]
    leg_k = [k for k in 1:K if ok[k] && fb[k] < 0.5]
    foot_k = [k for k in 1:K if ok[k] && fb[k] >= 0.5]
    lk = (exp(median(log.([P[k].kstar for k in leg_k]))), isempty(foot_k) ? 1.0 : exp(median(log.([P[k].kstar for k in foot_k]))))
    kappa(k, l) = exp((1 - fb[k]) * log(l[1]) + fb[k] * log(l[2]))
    resist(k, l) = (r = P[k].a + P[k].b * kappa(k, l); max(r, 0.1 * (P[k].a + P[k].b * P[k].kstar)))
    Q = zeros(K)
    sweeps = 0
    for outer in 1:100
        R = [ok[k] ? resist(k, lk) : Inf for k in 1:K]
        for inner in 1:2000
            sweeps += 1
            for k in 1:K
                Q[k] = ok[k] ? max(pin[k] - pout[k], 0.0) / R[k] : 0.0
            end
            Pa = pass_pressures(A, PASS1.a_order, ra, PATCHES.segment, Q, Float64(f.p_in), -1)
            qv = zeros(length(VEIN_PATCHES.segment))
            for k in 1:K
                ok[k] || continue
                for (i, m) in enumerate(P[k].m)
                    qv[m] += P[k].w[i] * Q[k]
                end
            end
            Pv = pass_pressures(W, PASS1.v_order, PASS1.rv, VEIN_PATCHES.segment, qv, Float64(f.p_out), +1)
            pvm = Pv[VEIN_PATCHES.segment]
            d = 0.0
            for k in 1:K
                ok[k] || continue
                new_in = Pa[PATCHES.segment[k]]
                new_out = sum(P[k].w .* pvm[P[k].m])
                d = max(d, abs(new_in - pin[k]), abs(new_out - pout[k]))
                pin[k] = 0.5 * pin[k] + 0.5 * new_in
                pout[k] = 0.5 * pout[k] + 0.5 * new_out
            end
            d < 1.0e-5 && break
        end
        up = sum(Q[k] * (1 - fb[k]) for k in 1:K)
        dn = sum(Q[k] * fb[k] for k in 1:K)
        # flow falls with κ as the microvessels' share of a patch's resistance, b κ / (a + b κ)
        share(ks) = sum(P[k].b * kappa(k, lk) / resist(k, lk) * Q[k] for k in ks) / sum(Q[k] for k in ks)
        ratio = (up / T_up, isempty(foot_k) ? 1.0 : dn / T_dn)
        quiet || say(@sprintf("network %d: κ above the ankle %.4g, below %.4g; flow ÷ target above %.5f, below %.5f", outer, lk[1], lk[2], ratio...))
        all(r -> abs(r - 1) < 1.0e-4, ratio) && break
        lk = (lk[1] * ratio[1]^(1 / share(leg_k)), isempty(foot_k) ? lk[2] : lk[2] * ratio[2]^(1 / share(foot_k)))
    end
    N = (; kappa_leg = lk[1], kappa_foot = lk[2], kappa = [ok[k] ? kappa(k, lk) : NaN for k in 1:K], pin, pout, Q, sweeps)
    quiet && return N
    save(path, N)
    say(@sprintf("network: maximal flow %.0f mL/min (%.0f above the ankle, %.0f below), %d sweeps", 60 * sum(Q), 60 * sum(Q .* (1 .- fb)), 60 * sum(Q .* fb), sweeps))
    return N
end

# ---------------------------------------------------------------------------------------------------------------------
# step 3: final

"Patch `k` at maximal dilation (the network's pressures and κ) and at rest (pass 1's resting pressures, the same κ, its own tone constant)."
function final_patch(k, W, N, Br)
    f = LIMB_FLOW
    t0 = time()
    P = deserialize(joinpath(FLOW_DIR, "pair_$k.jls"))
    vleaf = deserialize(joinpath(FLOW_DIR, "vleaf_$k.jls"))
    art = arterial(k)
    terms = leaves(art)
    terr = fill(F32(P.ml / length(terms)), length(terms))
    vf = venous_forest(k, W)
    vf.m == P.m || error("patch $k: its venous patches differ from the pairing's")
    dem = V.Territory.demand_field(art.b[:, terms], terr)
    kappa = F32(N.kappa[k])
    spec_s = V.FlowSpec(f; p_in = F32(N.pin[k]), p_out = F32(N.pout[k]))
    Cs = V.Flow.circuit(art, terr; state = :stress, veins = vf.tree, vleaf, phys = spec_s, kappa, demand = dem)
    pout_rest = sum(vf.n .* Br.pv[vf.m]) / sum(vf.n)
    spec_r = V.FlowSpec(f; p_in = F32(Br.pin[k]), p_out = F32(pout_rest))
    Cr = V.Flow.circuit(art, terr; state = :rest, veins = vf.tree, vleaf, phys = spec_r, kappa, demand = dem)
    roots = findall(iszero, vf.tree.parent)
    summ(C) = (;
        q = Float64(sum(C.resist.root_flow)), mbf = C.resist.mbf, so2 = C.resist.so2_sinus, reversed = reversed(C),
        act = C.resist.mean_act, max_tone = C.resist.n_max_tone, dilated = C.resist.n_dilated, starved = C.resist.starved,
        iters = C.resist.iters, shares = (C.resist.art_share, C.resist.cap_share, C.resist.ven_share), qroot = Float64.(C.qv[roots]),
    )
    save(joinpath(FLOW_DIR, "flows_$k.jls"), (stress = Cs.q[Cs.terminals], rest = Cr.q[Cr.terminals], so2_rest = Cr.so2))
    return save(joinpath(FLOW_DIR, "final_$k.jls"), (;
        k, ml = P.ml, nterm = P.nterm, kappa, pin = N.pin[k], pout = N.pout[k], pin_rest = Br.pin[k], pout_rest, m = vf.m,
        stress = summ(Cs), rest = summ(Cr), ctone = Cr.ctone, seconds = time() - t0,
    ))
end

# ---------------------------------------------------------------------------------------------------------------------
# step 4: checks

function checks(N, Bs, Br, ks)
    f = LIMB_FLOW
    F = Dict(k => deserialize(joinpath(FLOW_DIR, "final_$k.jls")) for k in ks if isfile(joinpath(FLOW_DIR, "final_$k.jls")))
    have = sort!(collect(keys(F)))
    content = 1.34 * f.hb / 100 * f.sao2
    lines = String[]
    out(s) = (push!(lines, s); say(s))
    mass = sum(F[k].ml for k in have) * f.rho
    qs = sum(F[k].stress.q for k in have) * 60
    qr = sum(F[k].rest.q for k in have) * 60
    so2 = sum(F[k].rest.q * F[k].rest.so2 for k in have) / sum(F[k].rest.q for k in have)
    vo2 = sum(F[k].rest.q * 60 * 1.34 * f.hb / 100 * (f.sao2 - F[k].rest.so2) for k in have)
    rest_g = f.mvo2 / (content * f.oef)
    out(@sprintf("patches solved: %d of %d (tissue %.0f of %.0f mL)", length(have), length(ks), mass / f.rho, SOFT_ML))
    out(@sprintf("maximal flow: %.0f mL/min (%.3f mL/min/g; the network %.0f); κ above the ankle %.4g, below %.4g", qs, qs / mass, 60 * sum(N.Q), N.kappa_leg, N.kappa_foot))
    out(@sprintf("resting flow: %.0f mL/min (%.4f mL/min/g; the O2 balance %.4f); measured 300 ± 30 (Chiesa 2016), 353 (Klein 2003)", qr, qr / mass, rest_g))
    out(@sprintf("flow reserve %.1f; femoral venous SO2 at rest %.3f (target %.3f); leg VO2 %.1f mL/min (18-28)", qs / qr, so2, f.sao2 * (1 - f.oef), vo2))
    for (label, s) in ((:maximal, :stress), (:rest, :rest))
        rv = sum(getfield(F[k], s).reversed for k in have)
        nt = sum(F[k].nterm for k in have)
        np = count(k -> getfield(F[k], s).reversed > 0, have)
        out(@sprintf("%s: reversed beds %d of %d (%.3f%%) in %d patches", label, rv, nt, 100 * rv / nt, np))
    end
    act = [F[k].rest.act for k in have]
    out(@sprintf(
        "tone at rest: mean activation per patch median %.2f (5-95%%: %.2f-%.2f); patches with a vessel at full tone %d; with a bed starved of O2 %d",
        median(act), quantile(act, 0.05), quantile(act, 0.95), count(k -> F[k].rest.max_tone > 0, have), count(k -> F[k].rest.starved > 0, have),
    ))
    for (lo, hi) in ((0.0, 0.03), (0.03, 0.4), (0.4, 0.7), (0.7, 1.01))
        ks_b = filter(k -> lo <= HEIGHT[k] < hi, have)
        isempty(ks_b) && continue
        g = k -> F[k].ml * f.rho
        ms = [60 * F[k].stress.q / g(k) for k in ks_b]
        mr = [60 * F[k].rest.q / g(k) for k in ks_b]
        out(@sprintf(
            "height %3.0f-%3.0f%%: %4d patches; maximal flow %.3f mL/min/g (5-95%%: %.3f-%.3f), rest %.4f; inlet pressure at maximal flow %.1f mmHg",
            100lo, 100min(hi, 1), length(ks_b), median(ms), quantile(ms, 0.05), quantile(ms, 0.95), median(mr), median(N.pin[ks_b]),
        ))
    end
    # how far the final flows move pass 1's pressures
    A, W = pass1.arteries.tree, pass1.veins.tree
    for (label, s, ra, pin0, pout0) in (("maximal", :stress, PASS1.ra_stress, N.pin, N.pout), ("rest", :rest, PASS1.ra_rest, Br.pin, nothing))
        Q = zeros(length(PATCHES.segment))
        qv = zeros(length(VEIN_PATCHES.segment))
        for k in have
            st = getfield(F[k], s)
            Q[k] = st.q
            for (i, m) in enumerate(F[k].m)
                qv[m] += st.qroot[i]
            end
        end
        Pa = pass_pressures(A, PASS1.a_order, ra, PATCHES.segment, Q, Float64(f.p_in), -1)[PATCHES.segment]
        Pv = pass_pressures(W, PASS1.v_order, PASS1.rv, VEIN_PATCHES.segment, qv, Float64(f.p_out), +1)[VEIN_PATCHES.segment]
        dpin = [Pa[k] - (s === :stress ? pin0[k] : F[k].pin_rest) for k in have]
        dpout = [
            let st = getfield(F[k], s)
                sum(st.qroot ./ sum(st.qroot) .* Pv[F[k].m]) - (s === :stress ? F[k].pout : F[k].pout_rest)
            end for k in have
        ]
        out(@sprintf(
            "%s: pass 1 with the final flows moves the inlet pressures by median %+.2f (max |%.2f|) mmHg and the outlet pressures by %+.2f (max |%.2f|)",
            label, median(dpin), maximum(abs, dpin), median(dpout), maximum(abs, dpout),
        ))
    end
    secs = [F[k].seconds for k in have]
    out(@sprintf("final solves: %.1f s per patch (median)", median(secs)))
    open(joinpath(FLOW_DIR, "summary.txt"), "w") do io
        foreach(l -> println(io, l), lines)
    end
    return nothing
end

# ---------------------------------------------------------------------------------------------------------------------

function main()
    W = venule_map()
    Bs, Br = equal_share(:stress), equal_share(:rest)
    ks = [k for k in eachindex(PATCHES.voxels) if PATCHES.voxels[k] > 0 && !isempty(W.by_k[k])]
    skipped = setdiff(eachindex(PATCHES.voxels), ks)
    isempty(skipped) || say("patches without tissue or venules, left out: $(skipped)")
    todo = filter(k -> !isfile(joinpath(FLOW_DIR, "pair_$k.jls")), ks)
    failed = run_patches(k -> pair_patch(k, W, Bs), todo, "pairing")
    isempty(failed) || say("pairing failed for $(length(failed)) patches: $(failed)")
    P = Any[isfile(joinpath(FLOW_DIR, "pair_$k.jls")) ? deserialize(joinpath(FLOW_DIR, "pair_$k.jls")) : nothing for k in eachindex(PATCHES.voxels)]
    rv = sum(p.reversed for p in P if p !== nothing)
    say(@sprintf("pairing done: %d patches; reversed beds at 1.25 κ: %d", count(!isnothing, P), rv))
    N = network(P)
    todo = filter(k -> P[k] !== nothing && !isfile(joinpath(FLOW_DIR, "final_$k.jls")), ks)
    failed = run_patches(k -> final_patch(k, W, N, Br), todo, "final")
    isempty(failed) || say("final failed for $(length(failed)) patches: $(failed)")
    checks(N, Bs, Br, ks)
    say(@sprintf("done in %.1f h", (time() - STARTED) / 3600))
end

function test()
    W = venule_map()
    Bs = equal_share(:stress)
    near(h) = argmin(k -> PATCHES.voxels[k] > 0 ? abs(HEIGHT[k] - h) : Inf, eachindex(PATCHES.voxels))
    used = findall(>(0), PATCHES.voxels)
    across = findfirst(k -> 0.2 < below(k) < 0.8, used)
    ks = unique([1147, 243, 242, used[across], near(0.45), near(0.8), used[argmin(PATCHES.voxels[used])], used[argmax(PATCHES.voxels[used])]])
    say("test patches: $(ks)")
    t1 = @elapsed one = [pair_core(k, W, Bs) for k in ks]
    for r in one
        s = r.summary
        say(@sprintf(
            "patch %4d: %6d terminals, %5.1f mL, height %.2f, below the ankle %.2f; venules ÷ terminals %.2f in %d venous patches; κ* %.4g, a %.4g, b %.4g (micro share at κ* %.2f); in %.1f / out %.1f mmHg; venous outflow ÷ inflow %.4f, smallest outlet share %.3f; reversed %d; %.1f s",
            s.k, s.nterm, s.ml, HEIGHT[s.k], below(s.k), sum(s.n) / s.nterm, length(s.m), s.kstar, s.a, s.b,
            s.b * s.kstar / (s.a + s.b * s.kstar), s.pin, s.pout, s.outflow, minimum(s.w), s.reversed, s.seconds,
        ))
    end
    say(@sprintf("one at a time: %.0f s for %d patches", t1, length(ks)))
    many = Vector{Any}(undef, length(ks))
    t2 = @elapsed @sync for (i, k) in enumerate(ks)
        Threads.@spawn many[i] = pair_core(k, W, Bs)
    end
    same = all(i -> one[i].summary.kstar == many[i].summary.kstar && one[i].vleaf == many[i].vleaf, eachindex(ks))
    say(@sprintf("all at once: %.0f s (%.1f× faster); results identical: %s", t2, t1 / t2, same))
end

# run as a script; `include`d, it only defines the steps (as the sensitivity tests use it)
if abspath(PROGRAM_FILE) == @__FILE__
    get(ARGS, 1, "") == "test" ? test() : main()
end
