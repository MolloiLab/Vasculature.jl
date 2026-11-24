# The whole body at rest: where the injected iodine goes, organ by organ (after transport_leg.jl).
#
#     julia -t auto whole_body.jl
#
# transport_leg.jl ran the injection through Vasculature.jl's whole-body model with the grown leg in the right leg's
# place, feeding the leg's femoral venous outflow back until it settled. This runs that final body once more and records
# what the package's `Body.simulate` does not return: the iodine (mg) in every body segment over time.
#   - The body: the same circulation (`leg_body`: the right leg's share of the cardiac output at the grown leg's resting
#     flow), the same injection, and the leg's saved femoral venous outflow in place of the body's five right-leg beds.
#   - `simulate_masses` is `Body.simulate`'s loop, copied step for step (the package is not changed), that also records
#     every 0.1 s the iodine in each segment's blood (plug part, mixing stages, a ventricle's pure delay) and in the space
#     around its cells, the iodine the kidneys have cleared and the iodine injected so far.
#   - The grown right leg: its iodine from its own step responses (transport_rest/result.jls), on the same times, by part.
#   - By organ: the blood flow at rest, and the iodine over time (arrival, peak, and at chosen times).
# Checks: the copy against `Body.simulate` (every outflow), the aorta against transport_leg.jl's final pass, the leg
# against its saved compartments, and the iodine balance (organs + cleared = injected).
#
# Numbers only. The body's heart and the grown leg are made from XCAT, so everything stays in LIMB_CACHE_DIR\whole_body.

using Printf, Serialization, Statistics

include(joinpath(@__DIR__, "transport_leg.jl")) # the leg's transport and its body (includes the notebook; runs nothing)

const WB_DIR = mkpath(joinpath(CACHE_DIR, "whole_body"))
const WB_LOG = open(joinpath(WB_DIR, "progress.log"), "a")
function wsay(s)
    line = @sprintf("[%6.0f s] %s", time() - STARTED, s)
    println(line)
    flush(stdout)
    println(WB_LOG, line)
    flush(WB_LOG)
    return nothing
end

const SHOWN = [5, 10, 15, 20, 30, 60, 120] # s: the times the tables give
const TISSUE_NAMES = (:muscle, :skin, :fat, :bone, :other)
tissues(region) = [Symbol(t, :_, region) for t in TISSUE_NAMES]
const CLEARED = "Cleared into urine"
const GROWN_LEG = "Right leg (grown)" # (the notebook's `LEG` is the leg's geometry)

"The organs and the body segments each one holds. The right leg is the grown one: it has no body segments."
const ORGANS = [
    ("Heart chambers", [:ra, :rv, :la, :lv]),
    ("Heart muscle", [:coronary, :coronary_sinus]),
    ("Lungs", [:pulm_trunk, :lung]),
    ("Aorta", [:aortic_root, :arch, :desc_aorta]),
    ("Brain", [:brain_r, :brain_l]),
    ("Kidneys", [:kidney]),
    ("Liver and gut", [:gut, :portal_vein, :hepatic_artery, :liver]),
    (GROWN_LEG, Symbol[]),
    ("Left leg", tissues(:leg_l)),
    ("Arms", [tissues(:arm_r); tissues(:arm_l)]),
    ("Rest of the body", [tissues(:upper); tissues(:lower)]),
    (
        "Large veins", [
            :arm_vein_r, :axillary_r, :subclavian_r, :ij_r, :brachiocephalic_r, :arm_vein_l, :axillary_l, :subclavian_l,
            :ij_l, :brachiocephalic_l, :svc, :femoral_vein_r, :common_iliac_r, :femoral_vein_l, :common_iliac_l,
            :ivc_infrarenal, :ivc_renal, :ivc_hepatic,
        ],
    ),
]

"The organs folded into eight groups for the share of the iodine (the categorical palette's eight slots, in this order)."
const WHERE = [
    ("Heart", ["Heart chambers", "Heart muscle"]),
    ("Lungs", ["Lungs"]),
    ("Legs", [GROWN_LEG, "Left leg"]),
    ("Liver and gut", ["Liver and gut"]),
    ("Kidneys and urine", ["Kidneys", CLEARED]),
    ("Brain", ["Brain"]),
    ("Arms and the rest of the body", ["Arms", "Rest of the body"]),
    ("Large vessels", ["Aorta", "Large veins"]),
]

# ---------------------------------------------------------------------------------------------------------------------
# the body, with the iodine in every segment

"""
`Body.simulate` (Vasculature.jl src/perfusion/body.jl), step for step, that also records every `every` s the iodine (mg)
in each segment's blood (plug part, mixing stages, a ventricle's pure delay) and in its extracellular space, the iodine the
kidneys have cleared and the iodine injected so far, all as they stand at the start of that step. The package's own
pieces (Pipe, Delay, Run, pass!, pipe_mass, _inj, _stages, _vmix) are used as they are; a `detailed` segment holds
nothing here (its model keeps its own iodine).
"""
function simulate_masses(C, inj; tmax, dt, every, detailed)
    Bd = V.Body
    F = Float32
    dt = F(dt); segs = C.segs; ns = length(segs); co = C.co; pf = 1.0f0 - C.hct
    nodes = unique(vcat([s.from for s in segs], [s.to for s in segs]))
    inj.site in nodes || error("injection site $(inj.site) is not a node; nodes: $(nodes)")
    nid = Dict(n => k for (k, n) in enumerate(nodes))
    from = [nid[s.from] for s in segs]; to = [nid[s.to] for s in segs]
    path = Int[]; node = inj.site
    while node !== C.inlet
        k = findfirst(i -> segs[i].from === node && segs[i].kind === :vein, 1:ns)
        k === nothing && break
        push!(path, k); node = segs[k].to
    end
    onpath = falses(ns); onpath[path] .= true
    dvol = zeros(F, ns)
    line = Bd.Pipe(inj.dead_space); line.c[1] = inj.prime
    R = Bd.Run[]
    for (i, s) in enumerate(segs)
        n = Bd._stages(s, C.q[i], co); vm = Bd._vmix(s, C.q[i], co)
        push!(R, Bd.Run(Bd.Pipe(s.plug * s.v), s.delay > 0 ? Bd.Delay(s.delay, dt) : nothing, zeros(F, n), fill(F(n == 0 ? 0 : vm / n), n), 0.0f0, 0.0f0))
    end
    nt = round(Int, tmax / dt); stride = max(1, round(Int, every / dt)); nk = nt ÷ stride + 1
    rec_out = zeros(F, nk, ns); tt = zeros(F, nk)
    blood = zeros(nk, ns); ecv = zeros(nk, ns); cleared = zeros(nk); given = zeros(nk)
    det = [get(detailed, s.name, nothing) for s in segs]
    nodemass = zeros(F, length(nodes)); nodeflow = zeros(F, length(nodes)); cnode = zeros(F, length(nodes))
    qin = zeros(F, ns); qout = zeros(F, ns)
    beds = [i for i in 1:ns if segs[i].kind === :bed]; others = [i for i in 1:ns if segs[i].kind !== :bed]
    inc = [[k for k in 1:ns if to[k] == from[i]] for i in 1:ns]
    share = [C.q[i] / sum(C.q[j] for j in 1:ns if from[j] == from[i]) for i in 1:ns]
    injected = 0.0f0; excreted = 0.0f0
    for step in 0:nt
        t = F(step * dt)
        before = injected # what is in the body (or cleared) at the start of this step
        ri, ci = Bd._inj(inj, t); injected += dt * ri * ci
        ri > 0 && (ci = Bd.pass!(line, ri * dt, ci))
        extra = ri
        for i in path
            s = segs[i]; qin[i] = C.q[i] * co + extra
            if s.c > 0
                out = dvol[i] / (s.r * s.c); dvol[i] += dt * (extra - out); extra = out
            end
            qout[i] = C.q[i] * co + extra
        end
        coe = co + extra
        for i in others
            onpath[i] && continue
            qin[i] = segs[i].kind === :vein ? C.q[i] * co : C.q[i] * coe
            qout[i] = qin[i]
        end
        for i in beds
            qin[i] = share[i] * sum(qout[k] for k in inc[i]; init = 0.0f0); qout[i] = C.q[i] * co
        end
        if step % stride == 0
            kk = step ÷ stride + 1; tt[kk] = t
            for i in 1:ns
                r = R[i]
                rec_out[kk, i] = r.out
                blood[kk, i] = Bd.pipe_mass(r.pipe) + sum(r.c .* r.vs; init = 0.0f0) + (r.delay === nothing ? 0.0f0 : sum(r.delay.buf) * qin[i] * dt)
                ecv[kk, i] = r.ce * segs[i].ve
            end
            cleared[kk] = excreted; given[kk] = before
        end
        fill!(nodemass, 0.0f0); fill!(nodeflow, 0.0f0)
        for i in 1:ns
            nodemass[to[i]] += qout[i] * R[i].out; nodeflow[to[i]] += qout[i]
        end
        s0 = nid[inj.site]; nodemass[s0] += ri * ci; nodeflow[s0] += ri
        for k in eachindex(cnode)
            cnode[k] = nodeflow[k] > 0 ? nodemass[k] / nodeflow[k] : 0.0f0
        end
        for i in 1:ns
            s = segs[i]; r = R[i]; cin = cnode[from[i]]; f = qin[i]
            if det[i] !== nothing
                r.out = F(det[i](t)); continue
            end
            c = Bd.pass!(r.pipe, f * dt, cin)
            r.delay === nothing || (c = Bd.pass!(r.delay, c))
            n = length(r.c)
            if n == 0
                r.out = c; continue
            end
            psn = s.ps / n; cln = s.cl / n; xs = 0.0f0
            for k in 1:n
                old = r.c[k]; x = psn * (old / pf - r.ce)
                fo = k == n ? qout[i] : f
                r.c[k] += dt * (f * (c - old) - x - cln * old / pf) / r.vs[k]
                k == n && (r.vs[k] += dt * (f - fo))
                excreted += dt * cln * old / pf; xs += x
                c = old
            end
            s.ve > 0 && (r.ce += dt * xs / s.ve)
            r.out = c
        end
    end
    held = sum(
        Bd.pipe_mass(R[i].pipe) + sum(R[i].c .* R[i].vs; init = 0.0f0) + R[i].ce * segs[i].ve +
            (R[i].delay === nothing ? 0.0f0 : sum(R[i].delay.buf) * qin[i] * dt) for i in 1:ns
    )
    return (; t = tt, names = [s.name for s in segs], out = rec_out, blood, ecv, cleared, given, injected, excreted, held)
end

# ---------------------------------------------------------------------------------------------------------------------
# the grown leg

"The grown right leg's iodine (mg) on the input's grid, from its saved step responses and input, as transport_leg.jl's `report` takes it."
function leg_iodine(R)
    A = cumin(R.c)
    nr = length(R.regions)
    G1 = R.fp .* A
    G2 = [cumflux(A, R.steps.pin[r]) for r in 1:nr]
    G3 = [cumflux(A, R.steps.bin[r]) for r in 1:nr]
    G4 = [cumflux(A, R.steps.bout[r]) for r in 1:nr]
    G5 = [cumflux(A, R.steps.vroot[r]) for r in 1:nr]
    G6 = cumflux(A, R.steps.fem)
    art1 = G1 .- sum(G2)
    ven1 = sum(G5) .- G6
    return (;
        art1, ven1, total = G1 .- G6, part = [G2[r] .- G5[r] for r in 1:nr], # a region: its patches' arteries, beds and veins
        art2 = [G2[r] .- G3[r] for r in 1:nr], beds = [G3[r] .- G4[r] for r in 1:nr], ven2 = [G4[r] .- G5[r] for r in 1:nr],
    )
end

# ---------------------------------------------------------------------------------------------------------------------
# by organ

"Blood flow at rest (mL/min) of each organ that takes a share of the cardiac output."
function organ_flows(BC)
    idx = Dict(s.name => i for (i, s) in enumerate(BC.segs))
    q(names) = 60 * BC.co * sum(BC.q[idx[n]] for n in names)
    rows = [
        ("Heart muscle", q([:coronary])),
        ("Brain", q([:brain_r, :brain_l])),
        ("Kidneys", q([:kidney])),
        ("Liver and gut", q([:liver])), # the gut's blood reaches the liver through the portal vein: the liver's flow is both
        (GROWN_LEG, q(tissues(:leg_r))),
        ("Left leg", q(tissues(:leg_l))),
        ("Arms", q([tissues(:arm_r); tissues(:arm_l)])),
        ("Rest of the body", q([tissues(:upper); tissues(:lower)])),
    ]
    total = sum(last, rows)
    abs(total / (60 * BC.co) - 1) < 1.0e-4 || error(@sprintf("the organs take %.1f mL/min of %.1f", total, 60 * BC.co))
    return rows
end

"Arrival (first time at 10% of the peak), the peak and its time, and the values at SHOWN, for the curve `y` on times `t`."
function curve_stats(t, y)
    i = argmax(y)
    j = findfirst(>=(0.1 * y[i]), y)
    at = [y[findfirst(x -> abs(x - s) < 1.0e-3, t)] for s in SHOWN]
    return (arrive = Float64(t[j]), peak = y[i], tpeak = Float64(t[i]), at)
end

fmt(v) = v >= 1000 ? replace(@sprintf("%.0f", v), r"(\d)(?=(\d{3})+$)" => s"\1,") : v >= 10 ? @sprintf("%.0f", v) : @sprintf("%.1f", v)

function wb_main()
    R = deserialize(joinpath(T_DIR, "result.jls"))
    ks = [k for k in eachindex(PATCHES.voxels) if isfile(joinpath(FLOW_DIR, "final_$k.jls"))]
    q = zeros(length(PATCHES.voxels))
    for k in ks
        q[k] = deserialize(joinpath(FLOW_DIR, "final_$k.jls")).rest.q
    end
    B = leg_body(sum(q)) # transport_leg.jl's `main`, so the body is the same to the bit
    B.share == R.body.share || error("the body's right-leg share $(B.share) differs from transport_rest's $(R.body.share)")
    BC = B.BC
    wsay(@sprintf("%d patches, the grown leg %.1f mL/min; the body's segments: %d", length(ks), 60 * sum(q), length(BC.segs)))

    # every body segment in one organ, or one of the right-leg beds the grown leg replaces
    names = [s.name for s in BC.segs]
    owner = Dict{Symbol, String}()
    for (o, segs) in ORGANS, s in segs
        s in names || error("$s is not a body segment")
        haskey(owner, s) && error("$s is in two organs")
        owner[s] = o
    end
    rest = [n for n in names if !haskey(owner, n) && !(n in LEG_BEDS)]
    isempty(rest) || error("body segments in no organ: $(rest)")

    # the body with the leg's saved outflow fed back, as transport_leg.jl's final pass, recorded
    back = t -> PF * interp(Float32(t), TC, R.cs)
    det = Dict{Symbol, Any}(b => back for b in LEG_BEDS)
    t0 = time()
    S = simulate_masses(BC, INJECTION; tmax = TMAX, dt = 0.002f0, every = 0.1f0, detailed = det)
    wsay(@sprintf("body with the iodine in every segment: %.0f s", time() - t0))
    r = V.Body.simulate(BC, INJECTION; tmax = TMAX, dt = 0.002f0, every = 0.1f0, detailed = det)
    same = maximum(maximum(abs.(S.out[:, i] .- r.out[n])) for (i, n) in enumerate(names))
    aorta = maximum(abs.(S.out[:, findfirst(==(:desc_aorta), names)] .- R.body.aorta)) / maximum(R.body.aorta)
    wsay(@sprintf(
        "checks: the copy vs Body.simulate, every outflow max |difference| %.1e mg/mL, held %.3f vs %.3f mg, cleared %.3f vs %.3f mg; the aorta vs transport_rest's final pass %.1e of its peak",
        same, S.held, r.held, S.excreted, r.excreted, aorta,
    ))

    # the grown leg on the body's times, and against its saved compartments every second
    L = leg_iodine(R)
    fr = 1:round(Int, 1 / FINE):NG
    legcheck = max(
        maximum(abs.(L.art1[fr] .- R.art1)), maximum(abs.(L.ven1[fr] .- R.ven1)), maximum(abs.(L.total[fr] .- R.total)),
        maximum(maximum(abs.(L.art2[k][fr] .- R.art2[k])) for k in eachindex(R.regions)),
        maximum(maximum(abs.(L.beds[k][fr] .- R.beds[k])) for k in eachindex(R.regions)),
        maximum(maximum(abs.(L.ven2[k][fr] .- R.ven2[k])) for k in eachindex(R.regions)),
    )
    on = 1:STRIDE:NG # the input's grid (0.02 s) to the body's records (0.1 s)
    length(on) == length(S.t) && maximum(abs.(GRID[on] .- S.t)) < 1.0e-3 || error("the leg's and the body's times differ")
    wsay(@sprintf("the grown leg vs its saved compartments, every second: max |difference| %.1e mg", legcheck))

    # by organ
    t = Float64.(S.t)
    no = length(ORGANS)
    Mb = zeros(length(t), no); Me = zeros(length(t), no)
    for (g, (o, segs)) in enumerate(ORGANS), s in segs
        i = findfirst(==(s), names)
        Mb[:, g] .+= S.blood[:, i]
        Me[:, g] .+= S.ecv[:, i]
    end
    gl = findfirst(o -> first(o) == GROWN_LEG, ORGANS)
    Mb[:, gl] .= L.total[on] # the grown leg: its vessels, beds and the space around its cells together (not split)
    M = Mb .+ Me
    organs = first.(ORGANS)
    held = vec(sum(M; dims = 2)) .+ S.cleared
    gap = S.given .- held
    dose = Float64(S.injected)
    ig = argmax(abs.(gap))
    wsay(@sprintf(
        "iodine balance: injected %.0f mg; organs + cleared differ from the injected so far by at most %.1f mg (%.2f%% of the dose) at %.1f s, %.1f mg at %.0f s",
        dose, gap[ig], 100 * abs(gap[ig]) / dose, t[ig], gap[end], t[end],
    ))

    flows = organ_flows(BC)
    co = 60 * Float64(BC.co)
    stats = [curve_stats(t, M[:, g]) for g in 1:no]
    cstats = curve_stats(t, S.cleared)
    # the eight groups, as % of the iodine injected so far at SHOWN
    col(name) = name == CLEARED ? S.cleared : M[:, findfirst(==(name), organs)]
    iS = [findfirst(x -> abs(x - s) < 1.0e-3, t) for s in SHOWN]
    where_mg = [sum(col(m)[i] for m in members) for (_, members) in WHERE, i in iS]
    where_pct = 100 .* where_mg ./ reshape(vec(sum(where_mg; dims = 1)), 1, :)
    # the grown leg by part
    parts = [[uppercasefirst(r) for r in R.regions]; "Main vessels"]
    leg_mg = [hcat([L.part[k][on] for k in eachindex(R.regions)]...) (L.art1 .+ L.ven1)[on]]

    lines = String[]
    out(s) = (push!(lines, s); wsay(s))
    out(@sprintf(
        "WHOLE BODY AT REST: %s mg of iodine injected into the right arm (35 mL at 350 mg/mL, then 35 mL of a 30:70 chaser, 5 mL/s); the body model with the grown right leg; %d feedback passes in transport_rest",
        fmt(dose), R.body.passes,
    ))
    out(@sprintf(
        "checks: copied loop vs Body.simulate %.1e mg/mL; aorta vs transport_rest %.1e of its peak; grown leg vs its saved compartments %.1e mg; balance (organs + cleared vs injected) within %.1f mg (%.2f%% of the dose)",
        same, aorta, legcheck, abs(gap[ig]), 100 * abs(gap[ig]) / dose,
    ))
    out("")
    out(@sprintf("BLOOD FLOW AT REST (mL/min; %% of the cardiac output, %.0f mL/min)", co))
    out(@sprintf("  %-28s %6.0f  100.0%%  (all the blood passes through them)", "Heart chambers and lungs", co))
    for (o, f) in sort(flows; by = last, rev = true)
        out(@sprintf("  %-28s %6.0f  %5.1f%%", o, f, 100 * f / co))
    end
    out("")
    out(@sprintf("IODINE IN EACH ORGAN (mg, then %% of the %s mg injected): arrival at 10%% of its peak, the peak, and at %s s", fmt(dose), join(SHOWN, ", ")))
    out(@sprintf("  %-20s %7s %16s %7s  %s", "organ", "arrives", "peak", "at", join([@sprintf("%11s", "$(s) s") for s in SHOWN])))
    for (g, o) in enumerate(organs)
        s = stats[g]
        out(@sprintf(
            "  %-20s %6.1fs %8s %5.1f%% %6.1fs  %s", o, s.arrive, fmt(s.peak), 100 * s.peak / dose, s.tpeak,
            join([@sprintf("%5s %4.1f%%", fmt(v), 100 * v / dose) for v in s.at]),
        ))
    end
    out(@sprintf("  %-20s %7s %8s %6s %7s  %s", CLEARED, "", "", "", "", join([@sprintf("%5s %4.1f%%", fmt(v), 100 * v / dose) for v in cstats.at])))
    out("")
    out("WHERE THE IODINE IS (% of the iodine injected so far)")
    out(@sprintf("  %-30s %s", "", join([@sprintf("%7s", "$(s) s") for s in SHOWN])))
    for (j, (w, _)) in enumerate(WHERE)
        out(@sprintf("  %-30s %s", w, join([@sprintf("%6.1f%%", where_pct[j, i]) for i in eachindex(SHOWN)])))
    end
    out(@sprintf("  %-30s %s", "all of it (mg)", join([@sprintf("%7s", fmt(sum(where_mg[:, i]))) for i in eachindex(SHOWN)])))
    out("")
    out("THE GROWN RIGHT LEG BY PART (mg)")
    out(@sprintf("  %-14s %s", "", join([@sprintf("%7s", "$(s) s") for s in SHOWN])))
    for (k, p) in enumerate(parts)
        out(@sprintf("  %-14s %s", p, join([@sprintf("%7s", fmt(leg_mg[i, k])) for i in iS])))
    end
    out(@sprintf("  %-14s %s", "Total", join([@sprintf("%7s", fmt(L.total[on][i])) for i in iS])))

    save(joinpath(WB_DIR, "result.jls"), (;
        t, organs, M, blood = Mb, outside = Me, cleared = S.cleared, given = S.given, dose, gap,
        flows, co, stats, cleared_stats = cstats, shown = SHOWN,
        where = (names = first.(WHERE), members = last.(WHERE), mg = where_mg, pct = where_pct),
        leg = (parts, mg = leg_mg, total = L.total[on], ml_min = 60 * sum(q)),
        checks = (copy = same, aorta, leg = legcheck, balance = gap[ig], balance_t = t[ig]),
        injection = "35 mL of contrast (350 mg iodine/mL) at 5 mL/s, then 35 mL of a 30:70 contrast:saline chaser, right arm",
    ))
    open(joinpath(WB_DIR, "summary.txt"), "w") do io
        foreach(l -> println(io, l), lines)
    end
    return wsay(@sprintf("done in %.0f s", time() - STARTED))
end

if abspath(PROGRAM_FILE) == @__FILE__
    wb_main()
end
