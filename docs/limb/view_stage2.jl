# A look at stage 2 of the leg once grow_stage2.jl has grown it: every patch of the arteries or the veins read from
# LIMB_CACHE_DIR, checked against its target and drawn over the whole leg, with the notebook's checks of pass 2 for every
# terminal of the leg, not one region: the terminal radius, and for the arteries the pressure left at each terminal at
# maximal flow, above the femoral vein's.
#
#     julia -t auto view_stage2.jl [veins]           read every patch (a few minutes), then draw; arteries unless veins
#     julia -t auto view_stage2.jl [veins] redraw    draw again from what the last read saved
#     julia view_stage2.jl both                      the arteries (red) and veins (blue) on one page, from both reads
#
# Needs XCAT_LEG_NRB and LIMB_CACHE_DIR, as grow_stage2.jl does (LIMB_NOTEBOOK overrides the notebook's path). Writes
# LIMB_CACHE_DIR/stage2_view (arteries) or stage2_view_veins: overview.html (a page that opens at once: the findings, the
# tree and the figure at screen size), overview.png (the full-size figure), overview_small.png, tree.png (the tree alone,
# for slides), summary.txt, patches.csv (a row per patch) and overview_data.jls (what the figures are drawn from); `both`
# writes stage2_view_both: overview.html, trees.png, trees_small.png and summary.txt. They are made from XCAT, so they
# stay on lab storage with the caches.

for name in ("XCAT_LEG_NRB", "LIMB_CACHE_DIR")
    haskey(ENV, name) || error("set $name first: without it the notebook falls back to another machine's paths")
end

using Base64: base64encode
import Dates

const NOTEBOOK = get(ENV, "LIMB_NOTEBOOK", joinpath(@__DIR__, "limb_pipeline.jl"))
"The tree looked at: the arteries unless `veins` is given; `both` puts the two on one page."
const SYSTEM = "both" in ARGS ? :both : "veins" in ARGS ? :veins : :arteries
"Each look's folder in LIMB_CACHE_DIR."
const VIEWS = (arteries = "stage2_view", veins = "stage2_view_veins", both = "stage2_view_both")
const OUT = joinpath(ENV["LIMB_CACHE_DIR"], VIEWS[SYSTEM])
"What the last read of `system` saved, to draw from."
saved(system) = joinpath(ENV["LIMB_CACHE_DIR"], VIEWS[system], "overview_data.jls")
const DATA = saved(SYSTEM)
const REDRAW = "redraw" in ARGS || SYSTEM === :both # both draws from the two reads

if REDRAW
    import Pkg
    Pkg.activate(joinpath(dirname(NOTEBOOK), ".."))
    import CairoMakie as Mke
    using Printf: @sprintf
    using Serialization: serialize, deserialize
    using Statistics: median
else
    include(NOTEBOOK) # its functions and cached stages, and CairoMakie as Mke
end
mkpath(OUT)

"Pass-2 vessels drawn over the whole leg: this radius (mm) and up."
const DRAW_RADIUS_MM = 0.1f0

"The edge (mm) of the cube in the thigh drawn with every vessel."
const CUBE_MM = 8.0f0

"The terminal density's bins: this many 1 mm columns a side."
const BIN = 2

"A density bin is shown when its columns hold this many tissue voxels; fewer hold too few terminals to count."
const MIN_VOXELS = 10

"Segments sampled from each patch for the Murray plot."
const MURRAY_SAMPLE = 150

"Bin width (mmHg) of the pressure left at the terminals."
const PRESSURE_BIN = 0.01

"What the text says of each tree. A look saved before the veins had one is of the arteries."
const WORDS = (
    arteries = (
        name = "arteries", Name = "Arteries", tree = "arterial", root = "femoral artery",
        one = "one tree from the femoral artery", reach = "to",
    ),
    veins = (
        name = "veins", Name = "Veins", tree = "venous", root = "femoral vein", one = "one tree into the femoral vein",
        reach = "from",
    ),
)
words(d) = WORDS[get(d, :system, :arteries)]

"A radius (µm) as the text gives it: 7, 8.63."
um(x) = @sprintf("%.3g", x)

"Each segment's children, stored once: those of `j` are `kids[start[j]:(start[j + 1] - 1)]`."
function children(parent)
    n = zeros(Int, length(parent) + 1)
    n[1] = 1
    for p in parent
        p > 0 && (n[p + 1] += 1)
    end
    start = cumsum(n)
    next = start[1:(end - 1)]
    kids = Vector{Int}(undef, start[end] - 1)
    for (j, p) in enumerate(parent)
        p > 0 || continue
        kids[next[p]] = j
        next[p] += 1
    end
    return start, kids
end

"""
    below(parent, weight) -> (order, start, kids, n)

The segments in breadth-first order from the roots, their children (`children`) and the final terminals below each
segment: a terminal `j` counts `weight(j)`, every other segment the sum over its children.
"""
function below(parent, weight)
    start, kids = children(parent)
    order = findall(iszero, parent)
    sizehint!(order, length(parent))
    i = 1
    while i <= length(order)
        j = order[i]
        append!(order, view(kids, start[j]:(start[j + 1] - 1)))
        i += 1
    end
    n = zeros(length(parent))
    for j in Iterators.reverse(order)
        n[j] = start[j] == start[j + 1] ? weight(j) : sum(c -> n[kids[c]], start[j]:(start[j + 1] - 1))
    end
    return order, start, kids, n
end

"""
Pressure lost (mmHg) from the root to the end of every segment per mL/s of flow to each final terminal, every vessel at
radius `radius(j)` (mm) and dilated: Poiseuille with the Pries–Secomb in vivo viscosity, as the notebook's `pressure_drop`.
"""
function loss(tree, order, n, radius, flow)
    P = zeros(length(n))
    for j in order
        L = sqrt(sum(d -> (tree.b[d, j] - tree.a[d, j])^2, 1:3))
        r = Float32(radius(j))
        R = V.Flow.poiseuille(Float32(L), r, flow.plasma_viscosity_mPas * V.Flow.mu_rel(2000.0f0 * r, flow.hd))
        p = tree.parent[j]
        P[j] = (p == 0 ? 0.0 : P[p]) + R * n[j]
    end
    return P
end

"Segments to draw: their ends in the coronal plane (x, z) and their radius (mm)."
seglist() = (ax = Float32[], az = Float32[], bx = Float32[], bz = Float32[], r = Float32[])

function push_segment!(l, tree, j, r)
    push!(l.ax, tree.a[1, j])
    push!(l.az, tree.a[3, j])
    push!(l.bx, tree.b[1, j])
    push!(l.bz, tree.b[3, j])
    push!(l.r, r)
    return l
end

"`a` with the entries of `b`, a named tuple of vectors with the same fields, appended."
function append_fields!(a, b)
    foreach(f -> append!(a[f], b[f]), keys(a))
    return a
end

"What the figure needs from the patches one task reads."
tally(bins) = (
    radius = zeros(Int, 5001), # terminal radii in bins of 0.01 µm, the last for 50 µm and over
    terminals = zeros(Int, bins), # terminals per density bin
    vessels = seglist(), # pass-2 vessels of DRAW_RADIUS_MM and up
    cube = seglist(), # every vessel with an end in the cube
    murray = (n = Float32[], r = Float32[]), # a sample of segments: terminals below, radius (µm)
)

"The 1 mm column (x, z) growth puts the point `p` (mm) in: its truncated coordinates, as the growth indexes points."
engine_column(p, grid) = (
    trunc(Int, p[1] / grid.h) - round(Int, grid.origin[1] / grid.h) + 1,
    trunc(Int, p[3] / grid.h) - round(Int, grid.origin[3] / grid.h) + 1,
)

in_cube(m, j, cube) = all(d -> cube.lo[d] <= m[d, j] <= cube.hi[d], 1:3)

"""
Read the patch cached at `path`, add it to `t` (a `tally`) and return its row: what it grew, and with `flow` (the
arteries) the pressure lost inside it to each of its terminals per unit flow (`loss`).
"""
function read_patch!(t, path, flow, grid, cube)
    p = deserialize(path)
    tree = p.tree
    order, start, _, n = below(tree.parent, _ -> 1.0)
    ends = [j for j in eachindex(tree.r) if start[j] == start[j + 1]]
    drop = flow === nothing ? Float32[] : Float32.(loss(tree, order, n, j -> tree.r[j], flow)[ends])
    cells = zeros(Int32, length(ends)) # each terminal's density bin (linear index), 0 outside the grid
    for (i, j) in enumerate(ends)
        t.radius[min(floor(Int, 1.0e5 * tree.r[j]) + 1, length(t.radius))] += 1
        x, z = engine_column(view(tree.b, :, j), grid)
        checkbounds(Bool, t.terminals, cld(x, BIN), cld(z, BIN)) || continue
        t.terminals[cld(x, BIN), cld(z, BIN)] += 1
        cells[i] = LinearIndices(t.terminals)[cld(x, BIN), cld(z, BIN)]
    end
    for j in eachindex(tree.r)
        tree.r[j] >= DRAW_RADIUS_MM && push_segment!(t.vessels, tree, j, tree.r[j])
        (in_cube(tree.a, j, cube) || in_cube(tree.b, j, cube)) && push_segment!(t.cube, tree, j, tree.r[j])
    end
    for j in 1:max(1, length(tree.r) ÷ MURRAY_SAMPLE):length(tree.r)
        push!(t.murray.n, n[j])
        push!(t.murray.r, 1000 * tree.r[j])
    end
    return (
        nterm = p.nterm, grown = length(ends), segments = length(tree.r), ml = p.ml, seconds = p.seconds,
        root_um = 1000 * p.root_radius, drop = drop, cells = cells,
    )
end

"""
Every patch of a tree (`system`) in `cache`, `<prefix>_<k>.jls`, read once, on all threads, into what the figure, the
summary and the table need. With `flow` (the arteries), also the pressure left at each terminal at maximal flow
(`pressure_left`).
"""
function collect_stage2(;
        system, cache, prefix, patches, tree1, murray1, grid, soft, flow, soft_ml, n_total, terminal_mm, gamma,
    )
    K = length(patches.segment)
    paths = [joinpath(cache, "$(prefix)_$k.jls") for k in 1:K]
    absent = findall(!isfile, paths)
    isempty(absent) ||
        error("$(length(absent)) patches are not grown ($prefix $(first(absent)) first): run grow_stage2.jl $system")
    origin = collect(Float64, grid.origin)
    dims = collect(Int, grid.dims)
    println("leg grid: ", join(dims, " × "), " voxels of ", grid.h, " mm")

    # the cube: 70% of the way up the leg (the thigh), in the middle of the patch there
    z = [tree1.b[3, j] for j in patches.segment]
    k_cube = argmin(abs.(z .- (minimum(z) + 0.7 * (maximum(z) - minimum(z)))))
    cells = [collect(Tuple(c)) for c in patches.box[k_cube] if patches.label[c] == k_cube]
    mid = sum(cells) ./ length(cells)
    centre = origin .+ (cells[argmin([sum(abs2, c .- mid) for c in cells])] .- 1) .* grid.h
    cube = (lo = Float32.(centre .- CUBE_MM / 2), hi = Float32.(centre .+ CUBE_MM / 2), patch = k_cube)

    bins = (cld(dims[1], BIN), cld(dims[3], BIN))
    rows = Vector{Any}(undef, K)
    done = Threads.Atomic{Int}(0)
    t0 = time()
    tasks = map(Iterators.partition(1:K, cld(K, 8 * Threads.nthreads()))) do ks
        Threads.@spawn begin
            t = tally(bins)
            for k in ks
                rows[k] = read_patch!(t, paths[k], flow, grid, cube)
                n = Threads.atomic_add!(done, 1) + 1
                n % 100 == 0 && println(@sprintf("read %4d of %d patches in %3.0f s", n, K, time() - t0))
            end
            t
        end
    end
    total = tally(bins)
    for t in fetch.(tasks)
        total.radius .+= t.radius
        total.terminals .+= t.terminals
        append_fields!(total.vessels, t.vessels)
        append_fields!(total.cube, t.cube)
        append_fields!(total.murray, t.murray)
    end
    patch_rows = identity.(rows)
    zmin, zmax = extrema(z)
    table = [
        (
            patch = k, height = Float64((z[k] - zmin) / (zmax - zmin)), target = r.nterm, grown = r.grown,
            segments = r.segments, ml = r.ml, root_um = r.root_um, seconds = r.seconds,
        ) for (k, r) in enumerate(patch_rows)
    ]

    pass1_lines = seglist()
    for j in setdiff(eachindex(tree1.r), patches.segment) # the pass-1 segments no patch replaced
        push_segment!(pass1_lines, tree1, j, Float32(murray1.radius[j]))
    end
    labelled = dropdims(sum(patches.label .> 0; dims = 2); dims = 2) # tissue voxels in each 1 mm column (x, z)
    tissue = [
        sum(view(labelled, (BIN * (i - 1) + 1):min(BIN * i, dims[1]), (BIN * (j - 1) + 1):min(BIN * j, dims[3])))
        for i in 1:bins[1], j in 1:bins[2]
    ]
    found = (
        system = system,
        grid = (origin = origin, h = Float64(grid.h), dims = dims),
        silhouette = Matrix(dropdims(any(soft; dims = 2); dims = 2)),
        tissue = tissue,
        terminals = total.terminals,
        pass1 = pass1_lines,
        vessels = total.vessels,
        cube = merge(cube, (lines = total.cube,)),
        murray = total.murray,
        murray1 = (n = Float32.(murray1.below), r = Float32.(1000 .* murray1.radius)),
        radius = total.radius,
        table = table,
        n_total = n_total,
        grown = sum(r.grown for r in patch_rows),
        terminal_um = 1000 * terminal_mm,
        gamma = Float64(gamma),
        tissue_ml = sum(labelled) / 1000,
        segments = length(tree1.r) - K + sum(r.segments for r in patch_rows),
    )
    flow === nothing && return found
    return merge(found, pressure_left(table, patch_rows, patches, tree1, murray1, flow, soft_ml, bins))
end

"""
The pressure left at each terminal at maximal flow, the notebook's check of one region, here for the whole leg: pass 1
at its Murray radii carries each patch's grown terminals, and every terminal an equal share of the leg's dilated flow.
Returns `collect_stage2`'s fields for it, with each patch's pressures added to its row of `table`.
"""
function pressure_left(table, patch_rows, patches, tree1, murray1, flow, soft_ml, bins)
    # pass 1 carries each patch's grown terminals; a patch's inlet is the end of its pass-1 segment's parent
    grown = [r.grown for r in patch_rows]
    weight = zeros(length(tree1.r))
    weight[patches.segment] .= grown
    order1, _, _, n1 = below(tree1.parent, j -> weight[j])
    P1 = loss(tree1, order1, n1, j -> murray1.radius[j], flow)
    q_max = flow.p_in / flow.r_total_dilated * flow.rho * soft_ml # mL/min, the whole leg dilated (as the notebook)
    unit = q_max / 60 / sum(grown) # mL/s to each terminal
    # the pressure left at each terminal, from -p_in (below 0 the equal share cannot reach it; lower ones go in bin 1),
    # and per density bin the terminals left below the femoral vein's pressure
    left_from = -Float64(flow.p_in)
    left = zeros(Int, ceil(Int, (flow.p_in - left_from) / PRESSURE_BIN) + 1)
    below_vein = zeros(Int, bins)
    lowest = (mmHg = Inf, patch = 0)
    under_vein = 0
    rows = Vector{Any}(undef, length(table))
    for (k, r) in enumerate(patch_rows)
        p = tree1.parent[patches.segment[k]]
        before = unit * (p == 0 ? 0.0 : P1[p]) # mmHg lost on the way to the patch
        kept = flow.p_in .- before .- unit .* r.drop
        for (x, c) in zip(kept, r.cells)
            left[clamp(floor(Int, (x - left_from) / PRESSURE_BIN) + 1, 1, length(left))] += 1
            x < flow.p_out && c > 0 && (below_vein[c] += 1)
        end
        under_vein += count(<(flow.p_out), kept)
        minimum(kept) < lowest.mmHg && (lowest = (mmHg = minimum(kept), patch = k))
        rows[k] = merge(
            table[k], (before = before, inside = unit * median(r.drop), lowest = minimum(kept), median = median(kept)),
        )
    end
    return (
        table = identity.(rows), left = left, left_from = left_from, below_vein = below_vein, lowest = lowest,
        under_vein = under_vein, flow = (p_in = Float64(flow.p_in), p_out = Float64(flow.p_out), q_max = Float64(q_max)),
    )
end

"The value below which a fraction `q` of a histogram's counts lie (bins of `width` from `from`)."
function hist_quantile(h, width, q; from = 0.0)
    target = q * sum(h)
    seen = 0
    for (i, n) in enumerate(h)
        seen += n
        seen >= target && return from + (i - 0.5) * width
    end
    return from + length(h) * width
end

"The summary, a line per finding: the growth (`growth_lines`, the figure's caption), then the arteries' pressure check."
summary_lines(d) = haskey(d, :left) ? [growth_lines(d); pressure_lines(d)] : growth_lines(d)

"What grew: the patches, the terminals against their target, the segments and the terminal radius."
function growth_lines(d)
    t = d.table
    w = words(d)
    share = [r.grown / r.target for r in t]
    least = argmin(share)
    rq(q) = hist_quantile(d.radius, 0.01, q)
    return [
        "Stage 2 of the leg, $(w.name): " * @sprintf(
            "all %d patches grown, %.1f h of growth (%.1f s per patch, median).",
            length(t), sum(r.seconds for r in t) / 3600, median([r.seconds for r in t]),
        ),
        @sprintf(
            "Terminals: %.2f million, target %.2f million; the patch with the fewest of its own (%d) has %.1f%%.",
            d.grown / 1.0e6, d.n_total / 1.0e6, t[least].patch, 100 * share[least],
        ),
        uppercasefirst(w.one) * @sprintf(": %.1f million segments.", d.segments / 1.0e6),
        @sprintf(
            "Terminal radius: median %.2f µm (5–95%%: %.2f–%.2f), target ", rq(0.5), rq(0.05), rq(0.95),
        ) * um(d.terminal_um) * " µm.",
    ]
end

"The pressure left at the terminals at maximal flow, the notebook's check of pass 2."
function pressure_lines(d)
    t = d.table
    pq(q) = hist_quantile(d.left, PRESSURE_BIN, q; from = d.left_from)
    starved = [r for r in t if r.median < d.flow.p_out] # patches whose median terminal is below the vein's pressure
    return [
        @sprintf("Maximal flow, every vessel dilated: %.1f L/min for the leg, an equal share to each terminal.", d.flow.q_max / 1000),
        @sprintf(
            "Pressure left at the terminals, of the femoral artery's %.0f mmHg: median %.1f, 5%% below %.1f, lowest %.1f (patch %d).",
            d.flow.p_in, pq(0.5), pq(0.05), d.lowest.mmHg, d.lowest.patch,
        ),
        @sprintf(
            "Lost on the way to a patch (pass 1): median %.1f, most %.1f mmHg; inside a patch, to its median terminal: median %.1f mmHg.",
            median([r.before for r in t]), maximum(r.before for r in t), median([r.inside for r in t]),
        ),
        @sprintf(
            "Terminals below the femoral vein's %.0f mmHg: %d (%.2f%%), in %d patches. %s",
            d.flow.p_out, d.under_vein, 100 * d.under_vein / d.grown, count(r -> r.lowest < d.flow.p_out, t),
            isempty(starved) ? "No patch has its median terminal below it." :
            @sprintf(
                "%d patches have their median terminal below it, all in the lowest %.0f%% of the leg's height.",
                length(starved), 100 * maximum(r.height for r in starved),
            ),
        ),
    ]
end

"A row per patch: what it grew, and its pressures where the table has them (the arteries)."
function write_table(path, table)
    pressure = haskey(first(table), :before)
    open(path, "w") do io
        names = [
            "patch", "height_up_the_leg", "target_terminals", "grown_terminals", "grown_percent", "segments", "tissue_ml",
            "root_radius_um", "growth_s",
        ]
        pressure && append!(names, [
            "pressure_lost_before_patch_mmHg", "pressure_lost_in_patch_median_mmHg", "pressure_left_lowest_mmHg",
            "pressure_left_median_mmHg",
        ])
        println(io, join(names, ","))
        for r in table
            print(io, @sprintf(
                "%d,%.3f,%d,%d,%.2f,%d,%.2f,%.1f,%.1f", r.patch, r.height, r.target, r.grown, 100 * r.grown / r.target,
                r.segments, r.ml, r.root_um, r.seconds,
            ))
            pressure && print(io, @sprintf(",%.2f,%.2f,%.2f,%.2f", r.before, r.inside, r.lowest, r.median))
            println(io)
        end
    end
    return path
end

# colours: surface and ink; vessel radius on a red ramp for the arteries and a blue one for the veins (thin light, thick
# dark; each red step as light as its blue one, in the palette red's hue, so the two trees weigh the same together);
# terminal density against design on blue (fewer), grey, red (more); categorical slots 1 and 2 for the two passes
const SURFACE = "#fcfcfb"
const INK = "#0b0b0b"
const INK_2 = "#52514e"
const MUTED = "#898781"
const BASELINE = "#c3c2b7"
const GRIDLINE = "#e1e0d9"
const LEG_FILL = "#ecebe7"
const SERIES = ("#2a78d6", "#eb6834")
const BLUES = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
const REDS = ["#f4948d", "#ef7d76", "#e96660", "#e44747", "#d43538", "#bd2e31", "#a82227", "#931c20", "#7f1218", "#690e13"]
"A tree's vessel ramp: red for the arteries, blue for the veins."
ramp(d) = get(d, :system, :arteries) === :veins ? BLUES : REDS
const DIVERGING = ["#104281", "#3987e5", "#9ec5f4", "#f0efec", "#f2b5b1", "#e34948", "#a32a2a"]

"""
Draw the segments of `l` (a `seglist`), thin first so thick vessels lie on top, coloured on `colormap` by log radius over
`limits`.
"""
function vessels!(ax, l, limits; scale, thinnest, colormap = BLUES)
    isempty(l.r) && return nothing
    order = sortperm(l.r)
    points = Vector{Mke.Point2f}(undef, 2 * length(order))
    for (i, j) in enumerate(order)
        points[2i - 1] = Mke.Point2f(l.ax[j], l.az[j])
        points[2i] = Mke.Point2f(l.bx[j], l.bz[j])
    end
    r = l.r[order]
    return Mke.linesegments!(
        ax, points; color = repeat(log10.(1000 .* r); inner = 2), colormap, colorrange = limits,
        linewidth = repeat(max.(thinnest, scale .* r); inner = 2),
    )
end

"""
Two trees' segments (`lists`, each on its own of `colormaps`) in one axis, thin first whichever tree they belong to, so
the wide vessels of both lie on top: the radius range cut in `bands`, each drawn for every tree in turn.
"""
function together!(ax, lists, colormaps, limits; scale, thinnest, bands = 24)
    band(r) = clamp(floor(Int, bands * (log10(1000r) - limits[1]) / (limits[2] - limits[1])) + 1, 1, bands)
    groups = map(l -> (b = band.(l.r); [findall(==(i), b) for i in 1:bands]), lists)
    for i in 1:bands, (l, g, colormap) in zip(lists, groups, colormaps)
        isempty(g[i]) || vessels!(ax, map(v -> v[g[i]], l), limits; scale, thinnest, colormap)
    end
    return nothing
end

"""
The figure: the leg's arteries or veins, their terminal density against design, a cube in the thigh with every vessel
and Murray's law, over the summary. Saved at twice screen size as `path`, and at screen size as `small` if given.
"""
function draw_overview(d, path; small = nothing)
    g = d.grid
    w = words(d)
    xs = g.origin[1] .+ (0:(g.dims[1] - 1)) .* g.h
    zs = g.origin[3] .+ (0:(g.dims[3] - 1)) .* g.h
    tall = 1150 # the leg panels' height (px)
    wide = round(Int, tall * g.dims[1] / g.dims[3])
    radius_limits = (log10(d.terminal_um), log10(1000 * maximum(d.pass1.r)))
    radius_ticks = filter(t -> radius_limits[1] <= log10(t) <= radius_limits[2], [10, 30, 100, 300, 1000, 3000])
    plain = (
        backgroundcolor = SURFACE, xgridvisible = false, ygridvisible = false, topspinevisible = false,
        rightspinevisible = false, leftspinecolor = BASELINE, bottomspinecolor = BASELINE, xtickcolor = BASELINE,
        ytickcolor = BASELINE, xticklabelcolor = INK_2, yticklabelcolor = INK_2, xlabelcolor = INK_2, ylabelcolor = INK_2,
        titlecolor = INK, subtitlecolor = INK_2, titlealign = :left,
    )
    chart = merge(plain, (xgridvisible = true, ygridvisible = true, xgridcolor = GRIDLINE, ygridcolor = GRIDLINE))
    bar = (labelcolor = INK_2, ticklabelcolor = INK_2, tickcolor = BASELINE, spinewidth = 0, height = Mke.Relative(0.6))

    fig = Mke.Figure(size = (2wide + 1100, tall + 340), backgroundcolor = SURFACE, fontsize = 14)
    Mke.Label(
        fig[1, 1:5], "The leg's $(w.name) after stage 2: $(w.one) $(w.reach) $(um(d.terminal_um)) µm";
        fontsize = 20, font = :bold, color = INK, halign = :left, tellwidth = false,
    )

    # the leg: pass 1 and every patch's vessels from DRAW_RADIUS_MM, and where the cube is
    leg = Mke.Axis(
        fig[2, 1]; title = "$(w.Name) ≥ $(round(Int, 1000DRAW_RADIUS_MM)) µm", subtitle = "pass 1 + $(length(d.table)) patches",
        aspect = Mke.DataAspect(), xlabel = "x (mm)", ylabel = "z (mm)", plain...,
    )
    Mke.heatmap!(leg, xs, zs, Float32.(d.silhouette); colormap = [SURFACE, LEG_FILL], colorrange = (0, 1))
    vessels!(leg, map(vcat, d.pass1, d.vessels), radius_limits; scale = 1.5, thinnest = 0.3, colormap = ramp(d))
    lo, hi = d.cube.lo, d.cube.hi
    middle = [Mke.Point2f((lo[1] + hi[1]) / 2, (lo[3] + hi[3]) / 2)]
    Mke.scatter!(leg, middle; marker = :rect, markersize = 16, color = :transparent, strokecolor = SURFACE, strokewidth = 4)
    Mke.scatter!(leg, middle; marker = :rect, markersize = 16, color = :transparent, strokecolor = INK, strokewidth = 1.5)
    Mke.Colorbar(
        fig[2, 2]; colormap = ramp(d), limits = radius_limits, ticks = (log10.(radius_ticks), string.(radius_ticks)),
        label = "vessel radius (µm)", bar...,
    )

    # terminals per mL of tissue in each column, against the design's even share
    design = d.grown / d.tissue_ml
    ratio = map((n, v) -> v >= MIN_VOXELS ? n / (v / 1000) / design : NaN, d.terminals, d.tissue)
    bx = g.origin[1] .+ (BIN .* (0:(size(ratio, 1) - 1)) .+ BIN / 2) .* g.h
    bz = g.origin[3] .+ (BIN .* (0:(size(ratio, 2) - 1)) .+ BIN / 2) .* g.h
    density = Mke.Axis(
        fig[2, 3]; title = "Terminal density", subtitle = @sprintf("÷ design, %.1f thousand per mL", design / 1000),
        aspect = Mke.DataAspect(), xlabel = "x (mm)", ylabel = "z (mm)", plain...,
    )
    Mke.heatmap!(
        density, bx, bz, ratio; colormap = DIVERGING, colorrange = (0.5, 1.5), lowclip = DIVERGING[1],
        highclip = DIVERGING[end], nan_color = :transparent,
    )
    Mke.contour!(density, xs, zs, Float32.(d.silhouette); levels = [0.5], color = MUTED, linewidth = 0.8)
    Mke.Colorbar(
        fig[2, 4]; colormap = DIVERGING, limits = (0.5, 1.5), lowclip = DIVERGING[1], highclip = DIVERGING[end],
        ticks = 0.5:0.25:1.5, label = "terminals per mL of tissue ÷ design", bar...,
    )

    foreach(ax -> Mke.limits!(ax, xs[1], xs[end], zs[1], zs[end]), (leg, density))

    side = Mke.GridLayout(fig[2, 5])
    cube = Mke.Axis(
        side[1, 1]; title = @sprintf("%.0f mm cube in the thigh", CUBE_MM),
        subtitle = "every vessel, down to $(um(d.terminal_um)) µm; the square in the first map",
        aspect = Mke.DataAspect(), xlabel = "x (mm)",
        ylabel = "z (mm)", plain...,
    )
    vessels!(cube, d.cube.lines, radius_limits; scale = 40, thinnest = 0.5, colormap = ramp(d))
    Mke.limits!(cube, lo[1], hi[1], lo[3], hi[3])

    murray = Mke.Axis(
        side[2, 1]; xscale = log10, yscale = log10, title = "Murray's law", subtitle = "a sample of every patch's segments",
        xlabel = "terminals downstream", ylabel = "radius (µm)", chart...,
    )
    Mke.scatter!(murray, d.murray.n, d.murray.r; markersize = 2, color = (SERIES[1], 0.25), strokewidth = 0)
    n = 10.0 .^ range(0, log10(d.n_total); length = 50)
    Mke.lines!(murray, n, d.terminal_um .* n .^ (1 / d.gamma); color = INK, linewidth = 1.5)
    Mke.scatter!(murray, d.murray1.n, d.murray1.r; markersize = 5, color = SERIES[2], strokewidth = 0)
    Mke.axislegend(
        murray,
        [
            Mke.MarkerElement(color = SERIES[1], marker = :circle, markersize = 8),
            Mke.MarkerElement(color = SERIES[2], marker = :circle, markersize = 8),
            Mke.LineElement(color = INK, linewidth = 1.5),
        ],
        ["pass 2", "pass 1", "r = $(um(d.terminal_um)) µm · n^(1/$(round(Int, d.gamma)))"];
        position = :lt, framevisible = false, labelcolor = INK_2, backgroundcolor = (SURFACE, 0.85),
    )

    Mke.Label(
        fig[3, 1:5], join(growth_lines(d), "\n"); halign = :left, justification = :left, color = INK_2, fontsize = 14,
        tellwidth = false,
    )
    foreach(c -> Mke.colsize!(fig.layout, c, Mke.Fixed(wide)), (1, 3))
    Mke.colsize!(fig.layout, 5, Mke.Fixed(560))
    Mke.rowsize!(fig.layout, 2, Mke.Fixed(tall))
    foreach(((c, gap),) -> Mke.colgap!(fig.layout, c, gap), zip(1:4, (12, 40, 12, 50)))
    Mke.rowsize!(side, 1, Mke.Auto(1.3))
    Mke.rowgap!(side, 28)
    Mke.resize_to_layout!(fig)
    Mke.save(path, fig; px_per_unit = 2)
    small === nothing || Mke.save(small, fig; px_per_unit = 1)
    return path
end

"A 50 mm scale bar in the bottom left corner of a leg drawn over `xs`, `zs` (mm), clear of the foot."
function scale_bar!(ax, xs, zs)
    x0, z0 = xs[1] + 5, zs[1] + 10
    Mke.lines!(ax, [x0, x0 + 50], [z0, z0]; color = INK, linewidth = 3)
    Mke.text!(ax, x0 + 25, z0 + 5; text = "50 mm", align = (:center, :bottom), color = INK, fontsize = 14)
    return nothing
end

"""
The tree alone, for slides: pass 1 and every patch's vessels from `DRAW_RADIUS_MM` in front view, coloured and sized by
radius, with a scale bar instead of axes.
"""
function draw_tree(d, path)
    g = d.grid
    w = words(d)
    xs = g.origin[1] .+ (0:(g.dims[1] - 1)) .* g.h
    zs = g.origin[3] .+ (0:(g.dims[3] - 1)) .* g.h
    tall = 1500 # the leg's height (px)
    wide = round(Int, tall * g.dims[1] / g.dims[3])
    limits = (log10(1000 * DRAW_RADIUS_MM), log10(1000 * maximum(d.pass1.r)))
    ticks = filter(t -> limits[1] <= log10(t) <= limits[2], [100, 200, 500, 1000, 2000, 3000])
    fig = Mke.Figure(size = (wide + 300, tall + 200), backgroundcolor = SURFACE, fontsize = 16)
    Mke.Label(fig[1, 1:2], "Leg $(w.tree) tree"; fontsize = 24, font = :bold, color = INK, halign = :left, tellwidth = false)
    ax = Mke.Axis(fig[2, 1]; aspect = Mke.DataAspect(), backgroundcolor = SURFACE)
    Mke.hidedecorations!(ax)
    Mke.hidespines!(ax)
    Mke.heatmap!(ax, xs, zs, Float32.(d.silhouette); colormap = [SURFACE, LEG_FILL], colorrange = (0, 1))
    vessels!(ax, map(vcat, d.pass1, d.vessels), limits; scale = 1.8, thinnest = 0.35, colormap = ramp(d))
    Mke.limits!(ax, xs[1], xs[end], zs[1], zs[end])
    scale_bar!(ax, xs, zs)
    Mke.Colorbar(
        fig[2, 2]; colormap = ramp(d), limits, ticks = (log10.(ticks), string.(ticks)), label = "Vessel radius (µm)",
        labelcolor = INK_2, ticklabelcolor = INK_2, tickcolor = BASELINE, spinewidth = 0, height = Mke.Relative(0.5),
    )
    Mke.Label(
        fig[3, 1:2],
        @sprintf("Front view, vessels of %d µm and up. ", round(Int, 1000DRAW_RADIUS_MM)) * uppercasefirst(w.one) *
        " $(w.reach) $(um(d.terminal_um)) µm terminals: " *
        @sprintf("%.0f million terminals, %.0f million segments.", d.grown / 1.0e6, d.segments / 1.0e6);
        fontsize = 14, color = INK_2, halign = :left, justification = :left, tellwidth = false, word_wrap = true,
    )
    Mke.colsize!(fig.layout, 1, Mke.Fixed(wide))
    Mke.rowsize!(fig.layout, 2, Mke.Fixed(tall))
    Mke.resize_to_layout!(fig)
    Mke.save(path, fig; px_per_unit = 1.5)
    return path
end

"""
The two trees in front view, vessels from `DRAW_RADIUS_MM`: the arteries (`a`) in red and the veins (`v`) in blue, each
darker and thicker where wider, alone and together, on one radius scale. Saved at 1.5 times screen size as `path`, and at
screen size as `small` if given.
"""
function draw_both(a, v, path; small = nothing)
    g = a.grid
    xs = g.origin[1] .+ (0:(g.dims[1] - 1)) .* g.h
    zs = g.origin[3] .+ (0:(g.dims[3] - 1)) .* g.h
    tall = 1500 # the legs' height (px)
    wide = round(Int, tall * g.dims[1] / g.dims[3])
    limits = (log10(1000 * DRAW_RADIUS_MM), log10(1000 * max(maximum(a.pass1.r), maximum(v.pass1.r))))
    ticks = filter(t -> limits[1] <= log10(t) <= limits[2], [100, 200, 500, 1000, 2000, 5000])
    trees = (arteries = map(vcat, a.pass1, a.vessels), veins = map(vcat, v.pass1, v.vessels))
    fig = Mke.Figure(size = (3wide + 500, tall + 240), backgroundcolor = SURFACE, fontsize = 16)
    Mke.Label(
        fig[1, 1:5], "Leg arteries (red) and veins (blue) after stage 2";
        fontsize = 24, font = :bold, color = INK, halign = :left, tellwidth = false,
    )
    panels = map(enumerate(("Arteries", "Veins", "Together"))) do (i, title)
        ax = Mke.Axis(
            fig[2, i]; title, titlealign = :left, titlecolor = INK, titlesize = 18, aspect = Mke.DataAspect(),
            backgroundcolor = SURFACE,
        )
        Mke.hidedecorations!(ax)
        Mke.hidespines!(ax)
        Mke.heatmap!(ax, xs, zs, Float32.(a.silhouette); colormap = [SURFACE, LEG_FILL], colorrange = (0, 1))
        ax
    end
    style = (scale = 1.8, thinnest = 0.35)
    vessels!(panels[1], trees.arteries, limits; style..., colormap = REDS)
    vessels!(panels[2], trees.veins, limits; style..., colormap = BLUES)
    together!(panels[3], (trees.veins, trees.arteries), (BLUES, REDS), limits; style...)
    foreach(ax -> Mke.limits!(ax, xs[1], xs[end], zs[1], zs[end]), panels)
    scale_bar!(panels[1], xs, zs)
    for (i, (colormap, label)) in enumerate(((REDS, "Artery radius (µm)"), (BLUES, "Vein radius (µm)")))
        Mke.Colorbar(
            fig[2, 3 + i]; colormap, limits, ticks = (log10.(ticks), string.(ticks)), label, labelcolor = INK_2,
            ticklabelcolor = INK_2, tickcolor = BASELINE, spinewidth = 0, height = Mke.Relative(0.5),
        )
    end
    tree_line(d) = (
        w = words(d);
        "$(w.Name): $(w.one) $(w.reach) " * @sprintf("%.0f million terminals of ", d.grown / 1.0e6) * um(d.terminal_um) *
            @sprintf(" µm, %.0f million segments.", d.segments / 1.0e6)
    )
    Mke.Label(
        fig[3, 1:5],
        @sprintf("Front view, vessels of %d µm and up; darker and thicker means wider. ", round(Int, 1000DRAW_RADIUS_MM)) *
            tree_line(a) * " " * tree_line(v);
        fontsize = 14, color = INK_2, halign = :left, justification = :left, tellwidth = false, word_wrap = true,
    )
    foreach(c -> Mke.colsize!(fig.layout, c, Mke.Fixed(wide)), 1:3)
    Mke.rowsize!(fig.layout, 2, Mke.Fixed(tall))
    Mke.resize_to_layout!(fig)
    Mke.save(path, fig; px_per_unit = 1.5)
    small === nothing || Mke.save(small, fig; px_per_unit = 1)
    return path
end

"A count as the pages give it: 171.3 million, 2.35 million, 1605."
million(x) = x >= 1.0e8 ? @sprintf("%.1f million", x / 1.0e6) : x >= 1.0e6 ? @sprintf("%.2f million", x / 1.0e6) : string(round(Int, x))

"A PNG as a data URL, so a page carries its figures in one file."
embed(file) = "data:image/png;base64," * base64encode(read(file))

"""
What the pages say of a tree: its terminal radius (median and 5–95%, µm), the share of tissue within ±10% of the design's
terminal density, its least-grown patch's share of its terminals, the terminals short of target (in all and at most in
one patch), its hours of growth and its root radius (mm, the widest pass-1 vessel).
"""
function key_figures(d)
    design = d.grown / d.tissue_ml
    ok = d.tissue .>= MIN_VOXELS
    ratio = (d.terminals ./ (d.tissue ./ 1000) ./ design)[ok]
    rq(q) = hist_quantile(d.radius, 0.01, q)
    short = [r.target - r.grown for r in d.table]
    return (
        median = rq(0.5), low = rq(0.05), high = rq(0.95),
        even = sum(d.tissue[ok][0.9 .<= ratio .<= 1.1]) / sum(d.tissue[ok]),
        least = minimum(r.grown / r.target for r in d.table), short = sum(max.(short, 0)), most_short = maximum(short),
        hours = sum(r.seconds for r in d.table) / 3600, per_patch = median([r.seconds for r in d.table]),
        root_mm = Float64(maximum(d.pass1.r)),
    )
end

"The pages' look: the palette's surfaces and ink, light and dark."
const PAGE_STYLE = """
:root { --bg: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink-2: #52514e; --muted: #898781; --line: #e1e0d9; --accent: #2a78d6; }
@media (prefers-color-scheme: dark) { :root { --bg: #0d0d0d; --surface: #1a1a19; --ink: #f0efec; --ink-2: #c3c2b7; --line: #2c2c2a; --accent: #3987e5; } }
body { margin: 0; background: var(--bg); color: var(--ink); font: 16px/1.7 system-ui, "Segoe UI", "Microsoft YaHei", sans-serif; }
main { max-width: 1200px; margin: 0 auto; padding: 24px 16px 48px; }
h1 { font-size: 24px; font-weight: 600; margin: 0 0 4px; }
h2 { font-size: 18px; font-weight: 600; margin: 32px 0 8px; }
.meta { color: var(--muted); font-size: 13px; margin: 0 0 20px; }
.lead { max-width: 75ch; margin: 0 0 20px; }
.cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 12px; }
.card { background: var(--surface); border: 1px solid var(--line); border-radius: 10px; padding: 12px 16px; }
.label { color: var(--ink-2); font-size: 13px; }
.value { font-size: 26px; font-weight: 600; font-variant-numeric: tabular-nums; }
.note { color: var(--muted); font-size: 13px; }
.table { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; background: var(--surface); border: 1px solid var(--line); font-variant-numeric: tabular-nums; }
th, td { text-align: left; padding: 8px 16px; border-bottom: 1px solid var(--line); vertical-align: top; }
th { font-size: 14px; font-weight: 600; }
td:first-child { color: var(--ink-2); }
tr:last-child td { border-bottom: none; }
.key { display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 8px; }
figure { margin: 0; }
figure img { display: block; width: 100%; height: auto; border: 1px solid var(--line); border-radius: 8px; background: #fcfcfb; }
figure.tree { max-width: 380px; margin: 0 auto; }
figure.trees { max-width: 960px; margin: 0 auto; }
figcaption { color: var(--muted); font-size: 13px; margin-top: 4px; }
a { color: var(--accent); }"""

"""
A page that opens at once, unlike the full-size figure: what grew, the tree alone (`tree`) and the figure at screen size
(`png`), both embedded so the page is one file, and links to the files beside it.
"""
function write_page(path, d, png, tree)
    t = d.table
    w = words(d)
    f = key_figures(d)
    cards = [
        (label = "Terminals", value = million(d.grown), note = @sprintf("target %s; the least-grown patch has %.1f%%", million(d.n_total), 100 * f.least)),
        (label = "Segments", value = million(d.segments), note = w.one),
        (label = "Terminal radius (median)", value = @sprintf("%.2f µm", f.median), note = @sprintf("5–95%%: %.2f–%.2f µm; target ", f.low, f.high) * um(d.terminal_um) * " µm"),
        (label = "Even coverage", value = @sprintf("%.1f%%", 100 * f.even), note = "of the tissue has a terminal density within ±10% of design"),
    ]
    card_html = join(
        ["<div class=\"card\"><div class=\"label\">$(c.label)</div><div class=\"value\">$(c.value)</div><div class=\"note\">$(c.note)</div></div>" for c in cards],
        "\n",
    )
    html = """
    <!doctype html>
    <html lang="en">
    <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Leg $(w.name), stage 2</title>
    <style>
    $(PAGE_STYLE)
    </style>
    </head>
    <body>
    <main>
    <h1>Leg $(w.name): stage 2 results</h1>
    <p class="meta">All $(length(t)) patches grown · made $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM")) · made from XCAT, so it stays on this machine or lab storage · <a href="../$(VIEWS.both)/overview.html">arteries and veins together</a></p>
    <div class="cards">
    $(card_html)
    </div>
    <h2>The $(w.tree) tree</h2>
    <figure class="tree">
    <a href="tree.png"><img src="$(embed(tree))" alt="The leg's $(w.tree) tree in front view, vessels of $(round(Int, 1000DRAW_RADIUS_MM)) µm and up, coloured by radius"></a>
    <figcaption>The whole tree in front view, vessels of $(round(Int, 1000DRAW_RADIUS_MM)) µm and up, coloured and sized by radius. Click for the full-size image (tree.png).</figcaption>
    </figure>
    <h2>Overview</h2>
    <figure>
    <a href="overview.png"><img src="$(embed(png))" alt="Stage 2 overview: two maps of the leg and two smaller plots"></a>
    <figcaption>Click for the full-size figure (overview.png; large, slow to open).</figcaption>
    </figure>
    <ul>
    <li>Left: $(w.name) of $(round(Int, 1000DRAW_RADIUS_MM)) µm and up (pass 1 and all $(length(t)) patches); darker and thicker means wider. The small square marks the cube on the right.</li>
    <li>Second: terminal density divided by design. Grey is even; blue, fewer; red, more.</li>
    <li>Top right: the $(round(Int, CUBE_MM)) mm cube in the thigh with every vessel, down to $(um(d.terminal_um)) µm.</li>
    <li>Bottom right: Murray's law; the points should lie on the black line.</li>
    </ul>
    <h2>Files in this folder</h2>
    <ul>
    <li><a href="overview.png">overview.png</a>: the full-size figure</li>
    <li><a href="tree.png">tree.png</a>: the tree alone</li>
    <li><a href="patches.csv">patches.csv</a>: one row per patch, opens in Excel</li>
    <li><a href="summary.txt">summary.txt</a>: the summary as text</li>
    </ul>
    </main>
    </body>
    </html>
    """
    write(path, html)
    return path
end

"""
The arteries (`a`) and veins (`v`) on one page: what each grew, side by side, the two trees drawn together (`png`,
embedded so the page is one file), and links to each tree's own page.
"""
function write_both_page(path, a, v, png)
    fa, fv = key_figures(a), key_figures(v)
    row(label, value) = "<tr><td>$label</td><td>$(value(a, fa))</td><td>$(value(v, fv))</td></tr>"
    rows = [
        row("Root", (d, f) -> @sprintf("%s, radius %.2f mm", words(d).root, f.root_mm)),
        row("Patches", (d, f) -> string(length(d.table))),
        row("Terminals", (d, f) -> "$(million(d.grown)) (target $(million(d.n_total)))"),
        row("Terminals short of target", (d, f) -> f.short == 0 ? "none" : "$(f.short) in all, at most $(f.most_short) in one patch"),
        row("Segments", (d, f) -> million(d.segments)),
        row("Terminal radius, median (5–95%)", (d, f) -> @sprintf("%.2f µm (%.2f–%.2f)", f.median, f.low, f.high)),
        row("Terminal radius, target", (d, f) -> um(d.terminal_um) * " µm"),
        row("Tissue within ±10% of design density", (d, f) -> @sprintf("%.1f%%", 100 * f.even)),
        row("Growth", (d, f) -> @sprintf("%.1f h, %.1f s per patch (median)", f.hours, f.per_patch)),
    ]
    most_short = max(fa.most_short, fv.most_short)
    lead = string(
        "Both trees are complete: ",
        most_short <= 0 ? "every patch of both grew all its terminals" :
            "every patch of both grew all but at most $most_short of its terminals",
        @sprintf(", and %.1f%% (arteries) and %.1f%% (veins) of the tissue", 100 * fa.even, 100 * fv.even),
        " has a terminal density within ±10% of design. The veins end in wider terminals (",
        um(v.terminal_um), " µm against ", um(a.terminal_um), " µm) because both trees end in the same number of ",
        @sprintf("terminals and the femoral vein is wider than the artery (%.2f mm against %.2f mm)", fv.root_mm, fa.root_mm),
        ": by Murray's law their terminals are wider by the same factor.",
    )
    key(colour) = "<span class=\"key\" style=\"background: $colour\"></span>"
    html = """
    <!doctype html>
    <html lang="en">
    <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Leg arteries and veins, stage 2</title>
    <style>
    $(PAGE_STYLE)
    </style>
    </head>
    <body>
    <main>
    <h1>Leg arteries and veins: stage 2 results</h1>
    <p class="meta">Both trees grown · made $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM")) · made from XCAT, so it stays on this machine or lab storage</p>
    <p class="lead">$(lead)</p>
    <div class="table">
    <table>
    <thead><tr><th></th><th>$(key(REDS[5]))Arteries</th><th>$(key(BLUES[5]))Veins</th></tr></thead>
    <tbody>
    $(join(rows, "\n"))
    </tbody>
    </table>
    </div>
    <h2>The two trees</h2>
    <figure class="trees">
    <a href="trees.png"><img src="$(embed(png))" alt="The leg's arteries in red and veins in blue, in front view: each alone, then both together"></a>
    <figcaption>Front view, vessels of $(round(Int, 1000DRAW_RADIUS_MM)) µm and up; darker and thicker means wider. Left, the arteries; middle, the veins; right, both together. Click for the full-size image (trees.png).</figcaption>
    </figure>
    <h2>Each tree in detail</h2>
    <ul>
    <li><a href="../$(VIEWS.arteries)/overview.html">Arteries</a>: terminal density, the $(round(Int, CUBE_MM)) mm cube with every vessel, Murray's law, a row per patch</li>
    <li><a href="../$(VIEWS.veins)/overview.html">Veins</a>: the same for the veins</li>
    </ul>
    <h2>Files in this folder</h2>
    <ul>
    <li><a href="trees.png">trees.png</a>: the two trees, full size</li>
    <li><a href="summary.txt">summary.txt</a>: the summary as text</li>
    </ul>
    </main>
    </body>
    </html>
    """
    write(path, html)
    return path
end

if SYSTEM === :both
    for s in (:arteries, :veins)
        isfile(saved(s)) || error("nothing saved for the $s yet: run view_stage2.jl$(s === :veins ? " veins" : "") first")
    end
    arteries, veins = deserialize(saved(:arteries)), deserialize(saved(:veins))
    findings = [growth_lines(arteries); ""; growth_lines(veins)]
    foreach(println, findings)
    write(joinpath(OUT, "summary.txt"), join(findings, "\n") * "\n")
    println("drawing...")
    println("trees: ", draw_both(arteries, veins, joinpath(OUT, "trees.png"); small = joinpath(OUT, "trees_small.png")))
    println("page: ", write_both_page(joinpath(OUT, "overview.html"), arteries, veins, joinpath(OUT, "trees_small.png")))
else
    if REDRAW
        isfile(DATA) || error("nothing to redraw yet: run without `redraw` first")
        data = deserialize(DATA)
    else
        started = time()
        s = SYSTEMS[SYSTEM] # the notebook's pass-2 inputs for this tree
        data = collect_stage2(;
            system = SYSTEM, cache = CACHE_DIR, prefix = s.prefix, patches = s.patches, tree1 = s.tree, murray1 = s.murray,
            grid = LEG.grid, soft = LEG.soft, flow = SYSTEM === :arteries ? LIMB_FLOW : nothing, soft_ml = SOFT_ML,
            n_total = N_TOTAL, terminal_mm = s.terminal_mm, gamma = GAMMA,
        )
        println(@sprintf("read every patch in %.0f s", time() - started))
        serialize(DATA, data)
    end
    findings = summary_lines(data)
    foreach(println, findings)
    write(joinpath(OUT, "summary.txt"), join(findings, "\n") * "\n")
    write_table(joinpath(OUT, "patches.csv"), data.table)
    println("drawing...")
    println("figure: ", draw_overview(data, joinpath(OUT, "overview.png"); small = joinpath(OUT, "overview_small.png")))
    println("tree: ", draw_tree(data, joinpath(OUT, "tree.png")))
    println("page: ", write_page(joinpath(OUT, "overview.html"), data, joinpath(OUT, "overview_small.png"), joinpath(OUT, "tree.png")))
end
