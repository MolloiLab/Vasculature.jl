### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 329820de-86c9-4483-8daf-4cdb31dad119
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))

    # a writable font cache for CairoMakie (a notebook worker may have no HOME)
    ENV["XDG_CACHE_HOME"] = get(ENV, "XDG_CACHE_HOME", joinpath(homedir(), ".cache"))
end;

# ╔═╡ b3d8dbc1-fce8-44ba-b4f4-062be2f1769e
using PlutoUI: TableOfContents, LocalResource

# ╔═╡ 7cd96038-c409-4a12-8367-f8f1c28c7ae8
using Printf: @sprintf

# ╔═╡ ed211fb3-e82c-4a06-ab74-8b8cf919e773
using Serialization: serialize, deserialize

# ╔═╡ 63e77379-a0ca-4b7f-9bda-a244fb1896be
using Statistics: mean, median, quantile

# ╔═╡ 2d88882b-057d-4f1d-a810-03f14c411fae
md"""
# Limb Pipeline

**From the raw XCAT leg to a vascular tree that reaches 7 µm terminal arterioles, grown in two passes.**

This notebook carries `docs/notebooks/02_limb_pipeline.jl` further: 7 µm terminals, both trees grown in full, and then the flow, the contrast transport, a narrowed artery and the iodine in the whole body.

The leg is built the way the heart (`00_cardiac_pipeline`) and the lungs (`01_pulmonary_pipeline`) are, with one difference of scale: a leg needs about 170 million terminals, far more than one growth can hold. So it grows in two passes. The first grows the whole leg to 10 mL patches; the second grows each patch to 7 µm on its own, and the patches run in parallel.

This page runs the stages that take seconds or minutes: the geometry, the first pass, and the second pass on the patches of one calf branch, which time the full second pass. The longer steps run patch by patch, over every patch, in the scripts next to this notebook (`docs/limb/`; its README gives the order, the commands and the times): the full second pass, the flow, the contrast transport at rest, a narrowed superficial femoral artery and the iodine in the whole body. Section 7 reads what they save. The phantom follows the same steps as the heart and lungs, and is laid out at the end.

| Stage | What happens | Run here |
|---|---|---|
| **1. Load** | XCAT's leg: its skin, bones, muscles and imaged femoral artery and vein, as NURBS surfaces | yes |
| **2. Prepare** | the soft tissue to grow into and how much of it is muscle, the bones as blockers, the imaged trees with penetrating stubs | yes |
| **3. Terminal count** | Murray's law fixes the count from the femoral artery and a 7 µm terminal arteriole | yes |
| **4. Pass 1** | the whole leg grown to 10 mL patches, arteries and veins | yes |
| **5. Patches** | each pass-1 terminal's own tissue | yes |
| **6. Pass 2** | one calf branch's patches grown to 7 µm terminals; the full pass timed | one branch |
| **7. Flow, transport, phantom** | the leg's own physiology; the flow and the contrast transport at rest, solved patch by patch; a narrowed superficial femoral artery; the iodine in the whole body; the femoral AIF, the phantom | the flow, the transport, the narrowing and the whole body read; the rest laid out |
| **8. Design choices** | what carries over from the heart and lungs | — |
| **9. Completing the build** | the remaining steps, in order | — |

Every heavy step is cached in `CACHE_DIR`, so a rerun takes seconds. The XCAT phantom is licensed: its files and anything made from them stay on lab storage.
"""

# ╔═╡ 4a662e6f-8a3c-4219-a3a6-03a7c848b45e
md"""
## Setup
"""

# ╔═╡ 939015ba-d0c0-4411-9223-131ae2d336e1
# ╠═╡ show_logs = false
import Vasculature as V

# ╔═╡ c91bff6a-765c-4bbd-93bd-b81f892a174e
# ╠═╡ show_logs = false
import CairoMakie as Mke

# ╔═╡ b24bcc05-e7c8-4d75-884e-876a5bf0b4b5
import GPUSelect

# ╔═╡ dcdca1e1-7737-4f0c-b019-17179dcd9530
TableOfContents()

# ╔═╡ 4e477957-3699-4072-99ee-464609f3d3f1
begin
    # licensed: this file and anything made from it stay on lab storage
    const LEG_NRB = get(
        ENV,
        "XCAT_LEG_NRB",
        joinpath(homedir(), "Documents", "dev", "MolloiLab", "xcat-data", "legs", "both_legs_1.nrb"),
    )
    # every cached stage (pass 1, the patches, each pass-2 patch); on a cluster, a shared scratch folder
    const CACHE_DIR = get(ENV, "LIMB_CACHE_DIR", joinpath(homedir(), "xcat-perfused", "limb_pipeline_cache"))
    mkpath(CACHE_DIR)
end;

# ╔═╡ c45cbbd6-6d31-4d1d-9d54-fc5dfb7870fb
# ╠═╡ show_logs = false
"Arrays to this machine's GPU (CUDA, Metal, ROCm, oneAPI), or the CPU if there is none. A backend in the environment counts only when its hardware is here: Metal, there for Apple silicon, also loads on machines that cannot run it."
const DEVICE = GPUSelect.Storage(
    let gpu = GPUSelect.find_gpu_pkgid()
        gpu !== nothing && gpu.name in GPUSelect.installed_backends() ? :GPU : :CPU
    end,
)

# ╔═╡ 77a6aaf7-fa17-470e-9e01-34bc52fbfd4f
"Run `compute` once and keep its result in `CACHE_DIR/name.jls`; afterwards load it."
function cached(compute, name)
    path = joinpath(CACHE_DIR, name * ".jls")
    isfile(path) && return deserialize(path)

    seconds = @elapsed result = compute()
    serialize(path, result)
    @info "$name: $(round(seconds; digits = 1)) s"
    return result
end

# ╔═╡ 8393831b-051f-4982-a98f-d28277891340
"The six face neighbours of a voxel."
const SIX_STEPS = (
    CartesianIndex(1, 0, 0),
    CartesianIndex(-1, 0, 0),
    CartesianIndex(0, 1, 0),
    CartesianIndex(0, -1, 0),
    CartesianIndex(0, 0, 1),
    CartesianIndex(0, 0, -1),
)

# ╔═╡ fc85e6c1-ab4a-43a8-aa91-91d6f74dba97
md"""
## 1. Load the Raw XCAT Leg

XCAT's legs come as one grouped NURBS file: every surface of both legs, each a set of patches under one name. The right leg is at negative x. This notebook reads its skin (`leg_right` and the foot), its bones, its imaged arteries (`arteries_rleg`) and veins (`veins_rleg`), and, for the tissue's makeup (section 2), its muscles (`musc…` and the foot's `rfoot_musc`). XCAT's femoral artery enters at the top of the leg with a radius of about 3.9 mm.

`Vasculature.XCAT.read_nrb` reads XCAT's single-object format, so this grouped one has its own reader below.
"""

# ╔═╡ f7984dcd-8660-4bb2-b99e-433d60b08960
begin
    const ENVELOPE = ["leg_right", "foot_right", "foot_right1", "foot_right2", "foot_right3", "foot_right4", "foot_right5"]
    const BONES = ["r_femur", "r_tibia", "r_fibula", "r_patella", "r_talus", "r_calcaneus", "r_navicular", "r_cuboid", "rpelvis1", "rpelvis2", "rpelvis3"]
    const LONG_BONES = ["r_femur", "r_tibia", "r_fibula"] # the fine blockers
    const ARTERIES = "arteries_rleg"
    const VEINS = "veins_rleg"

    "A surface's object name, without the `__occNNN` a repeated name gets."
    group(name) = replace(name, r"__occ\d+$" => "")
    "Foot bones XCAT does not name by side: the right ones are told apart by x. (`cuneiform_cart` is the larynx's.)"
    side_bone(name) = startswith(name, "metatarsal") || startswith(name, "phalanx") || occursin(r"^cuneiform\d", name)
    "The surfaces this notebook reads."
    wanted(name) = name in ENVELOPE || name in BONES || name in (ARTERIES, VEINS) || (side_bone(name) && !endswith(name, "_inner"))
end;

# ╔═╡ 0928ebc1-1d0b-4b61-a361-637c87839d35
"""
    read_grouped_nrb(path, keep) -> Vector{NurbsSurface}

XCAT's grouped NRB: per patch a name, an object code line, a patch-scale line, `M :M`, `N :N`, then the knots and control
points of `Vasculature.XCAT.read_nrb`'s format. A name repeated over several patches gets `__occNNN`. Only names `keep`
accepts are parsed.
"""
function read_grouped_nrb(path::AbstractString, keep)
    L = readlines(path)
    n = length(L)
    out = V.XCAT.NurbsSurface[]
    seen = Dict{String, Int}()
    i = 1
    while i <= n - 4
        if isempty(strip(L[i])) || !endswith(rstrip(L[i + 3]), ":M") || !endswith(rstrip(L[i + 4]), ":N")
            i += 1
            continue
        end
        name = String(strip(L[i]))
        M = parse(Int, split(L[i + 3])[1])
        N = parse(Int, split(L[i + 4])[1])
        i += 5
        startswith(strip(L[i]), "U Knot") || error("$name: expected U knots at line $i")
        i += 1
        u = Float64[]
        while !startswith(strip(L[i]), "V Knot")
            push!(u, parse(Float64, strip(L[i])))
            i += 1
        end
        i += 1
        v = Float64[]
        while !startswith(strip(L[i]), "Control")
            push!(v, parse(Float64, strip(L[i])))
            i += 1
        end
        i += 1
        k = get(seen, name, 0) + 1
        seen[name] = k
        if keep(name)
            P = Array{Float64, 3}(undef, 3, N, M)
            for q in 1:(M * N)
                f = split(L[i + q - 1], (',', ' ', '\t'); keepempty = false)
                jj = mod1(q, N)
                ii = cld(q, N)
                for d in 1:3
                    P[d, jj, ii] = parse(Float64, f[d])
                end
            end
            push!(out, V.XCAT.NurbsSurface(k == 1 ? name : @sprintf("%s__occ%03d", name, k), u, v, P))
        end
        i += M * N
    end
    return out
end

# ╔═╡ 5776c7ba-8194-47d8-a28e-046ef1205abb
leg_surfaces = cached("leg_surfaces") do
    S = read_grouped_nrb(LEG_NRB, wanted)
    # the unsided foot bones of the right leg (negative x)
    filter!(s -> !side_bone(group(s.name)) || mean(@view s.P[1, :, :]) < 0, S)
end;

# ╔═╡ 2c4e27bb-61d1-4d72-9dfa-8e496e27c231
let count_of(f) = count(s -> f(group(s.name)), leg_surfaces)
    Markdown.parse(
        """
        | object | patches |
        |---|---|
        | skin (leg and foot) | $(count_of(in(ENVELOPE))) |
        | bones | $(count_of(n -> n in BONES || side_bone(n))) |
        | imaged arteries (`$ARTERIES`) | $(count_of(==(ARTERIES))) |
        | imaged veins (`$VEINS`) | $(count_of(==(VEINS))) |
        """,
    )
end

# ╔═╡ ff0ae1c4-2c23-4fee-b563-57ddb6d3e5ec
md"""
## 2. Prepare

**The tissue.** The leg's soft tissue: XCAT's skin surfaces minus its bones and the imaged vessels' lumens (`TISSUE`), on a 1 mm grid. It is muscle, fat and skin together, one tissue for growth. XCAT's muscles (`MUSCLE`) give its makeup. Above the ankle it is the composite of `LIMB_FLOW` in section 7; the foot is mostly fat, tendon and skin, so at maximal dilation it takes less flow (`FOOT_MAX_SHARE`).

**The blockers.** No new vessel may cross bone: the femur, tibia and fibula at 0.4 mm, every bone at 1 mm out of the domain.

**The inlet trees.** XCAT's imaged artery and vein patches, as centrelines joined into one tree each, rooted where the femoral vessel enters the leg. Their nodes are spaced 2.8 mm × the organ's size scale, as for the heart and lungs. Every imaged vessel then gets a **penetrating stub** about every spacing where none is imaged (`V.XCAT.penetrating_stubs`), so growth starts along the whole vessel and not only from its few imaged ends. The imaged lumens of both trees, inside the skin and out, join the domain as voxels vessels may cross but where no terminal ends (`V.NOTERM`). They are blood, so the tissue's volume leaves them out.
"""

# ╔═╡ d57038e3-65cf-49fc-b365-1d3888b8b148
"A NURBS surface's triangle mesh, finer for surfaces with more control points."
mesh(s) = V.XCAT.tessellate(s; nu = clamp(4 * size(s.P, 2), 64, 480), nv = clamp(4 * size(s.P, 3), 64, 480))

# ╔═╡ 82d39e1f-f872-4449-838e-81678f93d077
"""
The leg on a 1 mm grid: its grid, the skin envelope, the bones, the soft tissue (envelope minus bone), and the fine
blockers (femur, tibia, fibula at 0.4 mm) as `(mask, origin, h)`.
"""
LEG = cached("leg") do
    skin = [mesh(s) for s in leg_surfaces if group(s.name) in ENVELOPE]
    bones = [s for s in leg_surfaces if group(s.name) in BONES || side_bone(group(s.name))]
    around = V.grid_around(reduce(hcat, first.(skin)), 1.0; margin = 3)
    grid = V.VoxelGrid(Float64.(round.(Int, around.origin)), 1.0, around.dims)

    envelope = reduce(.|, (V.voxelize(X, T, grid) for (X, T) in skin))
    bone = reduce(.|, (V.voxelize(mesh(s)..., grid) for s in bones))

    long = [mesh(s) for s in bones if group(s.name) in LONG_BONES]
    fine_grid = V.grid_around(reduce(hcat, first.(long)), 0.4; margin = 2)
    fine = reduce(.|, (V.voxelize(X, T, fine_grid) for (X, T) in long))

    (grid = grid, envelope = envelope, bone = bone, soft = envelope .& .!bone, blockers = (fine, fine_grid.origin, 0.4))
end;

# ╔═╡ bbfd5e8b-3f19-4db3-969f-eb73cb98ae4d
"""
XCAT's muscles in the soft tissue, on the leg's grid: every right-leg muscle surface (`musc…`, and the foot's
`rfoot_musc`) that reaches the grid, voxelized and joined. The rest of the soft tissue is fat, skin, tendon and
connective tissue. Section 7 takes the foot's maximal flow from it.
"""
MUSCLE = cached("muscle") do
    g = LEG.grid
    top = g.origin .+ g.h .* (g.dims .- 1)
    reaches(s) = all(d -> maximum(@view s.P[d, :, :]) >= g.origin[d] && minimum(@view s.P[d, :, :]) <= top[d], 1:3)
    muscle(name) = occursin(r"^musc\d+$", name) || name == "rfoot_musc"
    # the right leg's: negative x
    surfaces = filter(s -> mean(@view s.P[1, :, :]) < 0 && reaches(s), read_grouped_nrb(LEG_NRB, muscle))
    mask = falses(size(LEG.soft))
    for s in surfaces
        mask .|= V.voxelize(mesh(s)..., g)
    end
    mask .& LEG.soft
end;

# ╔═╡ 41d670f3-4243-48ce-b340-b0a6b8115b57
"""
The ankle: the narrowest slice of soft tissue between the foot's widest slice (in the lowest tenth of the leg) and a
quarter of the way up. `slice` is its index along z, `z` the height (mm) the foot lies below, `mm` that height above
the sole, and `ml` the soft tissue below it.
"""
ANKLE = let soft = LEG.soft, g = LEG.grid
    area = [count(view(soft, :, :, k)) for k in axes(soft, 3)]
    ks = findall(>(0), area)
    sole, top = first(ks), last(ks)
    widest = sole - 1 + argmax(view(area, sole:(sole + (top - sole) ÷ 10)))
    slice = widest - 1 + argmin(view(area, widest:(sole + (top - sole) ÷ 4)))
    (slice = slice, z = g.origin[3] + (slice - 1.5) * g.h, mm = (slice - sole) * g.h, ml = sum(area[1:(slice - 1)]) * g.h^3 / 1000)
end

# ╔═╡ 1b2a4544-4b53-4564-8039-0c4feac36a0e
let soft = LEG.soft, ml = LEG.grid.h^3 / 1000, leg = axes(LEG.soft, 3), up = ANKLE.slice:size(LEG.soft, 3), down = 1:(ANKLE.slice - 1)
    volume(r) = count(view(soft, :, :, r)) * ml
    muscle(r) = 100 * count(view(MUSCLE, :, :, r)) / count(view(soft, :, :, r))
    Markdown.parse(
        @sprintf(
            """
            | soft tissue | volume | XCAT's muscle |
            |---|---|---|
            | the whole leg | %.0f mL | %.0f%% |
            | above the ankle | %.0f mL | %.0f%% |
            | below the ankle (the foot, up to %.0f mm above the sole) | %.0f mL | %.0f%% |

            Above the ankle the soft tissue is the composite of section 7 (muscle 0.62). The foot is not: most of it is fat, tendon and skin, so at maximal dilation it takes less flow (`FOOT_MAX_SHARE`, section 7).
            """,
            volume(leg), muscle(leg), volume(up), muscle(up), ANKLE.mm, volume(down), muscle(down),
        ),
    )
end

# ╔═╡ ca00e0f2-712a-4f9e-8b82-2ba6543e3f18
"Whether the point `p` (mm) lies in the voxels of `mask` on the leg grid."
function inside(mask, p; grid = LEG.grid)
    k = round.(Int, (p .- grid.origin) ./ grid.h) .+ 1
    return all(1 .<= k .<= size(mask)) && mask[k...]
end

# ╔═╡ 2c77ebbd-5cc6-4081-aacf-1f6e10e9a78a
"""
Node spacing of the imaged centrelines (mm): the heart's 2.8 mm × the leg's size scale (its soft tissue over the
148.1 mL reference, cube root). A branch leaves its host on the host's centreline, and the collision test exempts only
the host, so nodes closer than the host's final radius would make every new branch hit the next imaged piece.
"""
const STEP = 2.8 * cbrt(count(LEG.soft) / 1000 / 148.1)

# ╔═╡ fa99ca18-b2a3-4a33-8cc1-54c84219aced
"A centreline resampled every `step` mm by arc length."
function resample(c::V.XCAT.Centerline, step::Real)
    arc = [0.0; cumsum([sqrt(sum(abs2, c.X[:, k + 1] .- c.X[:, k])) for k in 1:(size(c.X, 2) - 1)])]
    n = max(2, round(Int, arc[end] / step) + 1)
    X = zeros(3, n)
    r = zeros(n)
    for (q, t) in enumerate(range(0, arc[end], length = n))
        k = clamp(searchsortedlast(arc, t), 1, length(arc) - 1)
        w = (t - arc[k]) / max(arc[k + 1] - arc[k], eps())
        X[:, q] .= (1 - w) .* c.X[:, k] .+ w .* c.X[:, k + 1]
        r[q] = (1 - w) * c.r[k] + w * c.r[k + 1]
    end
    return V.XCAT.Centerline(c.name, X, r)
end

# ╔═╡ 6536e61a-f90e-4a09-973b-0f3f15bf2542
"""
    imaged_tree(surfaces) -> (tree, r_root, attach_tol)

The imaged patches of one vessel system joined into one tree: each patch's centreline (u and v swapped where the patch
runs the other way), the highest patch as the root, cut where it enters the leg's skin, resampled every `STEP` mm and
joined with the smallest attachment tolerance that joins them all. `r_root` is the root's mean radius over its first
10 mm.
"""
function imaged_tree(surfaces)
    lines = map(surfaces) do s
        c = V.XCAT.centerline(s; step = 1.0)
        # parametrized the other way round (u along the vessel): swap u and v
        median(c.r) > 15 ? V.XCAT.centerline(V.XCAT.NurbsSurface(s.name, s.v, s.u, permutedims(s.P, (1, 3, 2))); step = 1.0) : c
    end
    top = argmax([maximum(@view l.X[3, :]) for l in lines])
    c = lines[top]
    c.X[3, 1] < c.X[3, end] && (c = V.XCAT.Centerline(c.name, c.X[:, end:-1:1], reverse(c.r)))
    k0 = findfirst(k -> inside(LEG.envelope, c.X[:, k]), axes(c.X, 2))
    lines[top] = V.XCAT.Centerline(c.name, c.X[:, k0:end], c.r[k0:end])
    r_root = mean(lines[top].r[1:min(10, end)])
    pieces = [resample(l, STEP) for l in lines]
    for tol in (3.0, 5.0, 8.0, 12.0, 20.0)
        try
            return (tree = only(V.XCAT.base_forest(pieces; roots = [c.name], attach_tol = tol)), r_root, attach_tol = tol)
        catch e
            (e isa ErrorException && occursin("base_forest", e.msg)) || rethrow()
        end
    end
    return error("the imaged $(group(c.name)) patches do not join into one tree within 20 mm")
end

# ╔═╡ 4fb8cda2-0411-4a9e-a2f5-f4156b25a698
"""
The growth inputs of one vessel system: its imaged tree with penetrating stubs, the stub report, the root radius, and
its lumen voxels, which vessels may cross but where no terminal ends: `lumen` in the soft tissue, at the vessels' own
radii, and `pass` outside it, 1.5 mm wider.
"""
inlets = cached("inlets") do
    map((arteries = ARTERIES, veins = VEINS)) do name
        base = imaged_tree([s for s in leg_surfaces if group(s.name) == name])
        tree, report = V.XCAT.penetrating_stubs(
            base.tree,
            LEG.grid,
            LEG.soft;
            scale = STEP / 2.8,
            forbid = LEG.bone,
            blockers = (LEG.bone, LEG.grid.origin, LEG.grid.h),
        )
        lumen = V.XCAT.vessel_tube(tree, LEG.grid; margin = 0.0) .& LEG.soft
        pass = V.XCAT.vessel_tube(tree, LEG.grid; margin = 1.5) .& .!LEG.soft .& .!LEG.bone
        (tree = tree, report = report, r_root = base.r_root, attach_tol = base.attach_tol, lumen = lumen, pass = pass)
    end
end;

# ╔═╡ 5c1e8f3a-7b42-4d96-a0e5-3f9b2d6c8e17
begin
    "The imaged lumens in the soft tissue, both trees': blood, not tissue."
    const LUMENS = inlets.arteries.lumen .| inlets.veins.lumen
    "The tissue: the soft tissue less the imaged lumens. Terminals end only in it, and the patches and the flow count only it."
    const TISSUE = LEG.soft .& .!LUMENS
end;

# ╔═╡ d3bfa261-2daf-42e5-9bfe-4bc1d4a72f3c
"The tissue's volume (mL): what every growth pass and the flow are sized by."
const SOFT_ML = count(TISSUE) / 1000

# ╔═╡ 6e3760ce-8bde-4094-b3ba-d19d5c6cb0cb
let a = inlets.arteries, v = inlets.veins
    Markdown.parse(
        @sprintf(
            """
            | | arteries | veins |
            |---|---|---|
            | root radius (mm) | %.2f | %.2f |
            | imaged nodes | %d | %d |
            | penetrating stubs | %d of %d sites | %d of %d sites |
            | lumen in the soft tissue (mL) | %.0f | %.0f |
            | lumen outside it, 1.5 mm wider (mL) | %.0f | %.0f |

            Tissue %.0f mL: the skin envelope's %.0f mL less %.0f mL of bone and %.0f mL of imaged lumen. Centrelines resampled every %.1f mm. XCAT's femoral artery is %.1f mm across; its femoral vein (%.1f mm) is narrower than the measured supine 13.6 ± 3.0 mm (Keiler et al. 2019).
            """,
            a.r_root, v.r_root,
            length(a.tree.parent), length(v.tree.parent),
            a.report.stubs, a.report.sites, v.report.stubs, v.report.sites,
            count(a.lumen) / 1000, count(v.lumen) / 1000,
            count(a.pass) / 1000, count(v.pass) / 1000,
            SOFT_ML, count(LEG.envelope) / 1000, count(LEG.bone .& LEG.envelope) / 1000, count(LUMENS) / 1000, STEP,
            2a.r_root, 2v.r_root,
        ),
    )
end

# ╔═╡ 6193958a-729a-4632-b491-4e2c0ca5edf1
"The leg's front (x–z) and side (y–z) silhouettes, with segments drawn over them: one `linesegments!` per view."
function leg_figure(trees; title = "", size = (900, 900))
    fig = Mke.Figure(; size)
    g = LEG.grid
    xs = [g.origin[d] .+ (0:(g.dims[d] - 1)) .* g.h for d in 1:3]
    for (column, (view, label)) in enumerate(((1, "front (x–z)"), (2, "side (y–z)")))
        ax = Mke.Axis(fig[1, column]; title = label, aspect = Mke.DataAspect(), xlabel = view == 1 ? "x (mm)" : "y (mm)", ylabel = "z (mm)")
        silhouette = dropdims(any(LEG.soft; dims = 3 - view); dims = 3 - view)
        Mke.heatmap!(ax, xs[view], xs[3], silhouette; colormap = [:white, :gray88])
        for (segments, colour, scale) in trees
            a, b, r = segments
            points = Mke.Point2f[]
            for j in eachindex(r)
                push!(points, Mke.Point2f(a[view, j], a[3, j]), Mke.Point2f(b[view, j], b[3, j]))
            end
            widths = repeat(max.(0.4, scale .* r); inner = 2)
            Mke.linesegments!(ax, points; color = colour, linewidth = widths)
        end
        Mke.hidedecorations!(ax; label = false, ticklabels = false, ticks = false)
    end
    Mke.Label(fig[0, :], title; fontsize = 18, font = :bold)
    return fig
end

# ╔═╡ 2cc35552-467a-418b-9ebd-cd4a927e5b61
"An inlet tree's segments as `(a, b, r)`: each node with its parent."
function base_segments(tree)
    k = 2:length(tree.parent)
    return (tree.X[:, tree.parent[k]], tree.X[:, k], tree.r[k])
end

# ╔═╡ 826f47f1-d661-4ae1-83e7-a6a0fe1c1e62
leg_figure(
    [(base_segments(inlets.veins.tree), :royalblue, 0.6), (base_segments(inlets.arteries.tree), :firebrick, 0.6)];
    title = "Imaged arteries (red) and veins (blue), with their penetrating stubs",
)

# ╔═╡ 2a5c666d-73d6-49aa-b3c2-285b61e90e26
md"""
## 3. The Terminal Count

**Murray's law fixes the count.** With the bifurcation exponent γ = 3 (Murray) and every terminal carrying an equal share of the flow, a tree's terminal radius is

```math
r_\text{terminal} = r_\text{root} \, N^{-1/3}.
```

The femoral artery's radius is imaged. A terminal arteriole, the last arteriole before the capillaries, is about 7 µm in radius in skeletal muscle when fully dilated (below). So the leg needs N = (r\_root / 7 µm)³ terminals: about 170 million, about twice the 82 million that the heart's terminal density (957,045 in 148.1 mL) would give over this leg.

**Why 7 µm.** A grown radius is a vessel's passive size, without tone, so the terminal is the last arteriole before the capillaries at full dilation. No study has measured it in human muscle; three measurements in animal skeletal muscle agree on about 7 µm:
- **Rat spinotrapezius**: terminal arterioles are 7–13 µm in internal diameter at rest and widen by up to 55% after a tetanic contraction, as much as with papaverine (Marshall & Tandon 1984). Widened by 55%, that is 5.4–10.1 µm in radius.
- **Hamster skeletal muscle**: the arterioles that give off capillaries (Strahler order 1) are 9.0 µm in diameter at their widest in vasomotion, and the next order widens 1.55× with L-arginine (Lapi et al. 2019). If order 1 widens as much, its radius is 7.0 µm.
- **Hamster retractor**: a microvascular unit, one terminal arteriole and its capillaries, is 1 × 0.5 × 0.1 mm, 0.05 mm³ (Emerson & Segal 1997). One terminal per unit over this leg's soft tissue, or over its muscle alone, gives a terminal radius by Murray's law of 6.1–7.2 µm.

Fat and skin take the same terminal radius.

**Why γ = 3 and not the coronary 2.33.** A leg grown with the coronary exponent has conduits at 0.54–0.76 of XCAT's imaged radii. They take almost all of the arteriovenous pressure drop, and 39% of the capillary beds reverse. With γ = 3 the imaged vessels get 0.87–0.95 of their XCAT radii and no bed reverses. Lowering γ does give small terminals (9 µm at γ = 2.33), but only by starving them of pressure: a 10 mL block of thigh muscle grown that way loses about 800 mmHg at maximal flow, against 8 mmHg at γ = 3 with the count for 10 µm terminals.

**Why two passes.** 170 million terminals means 340 million segments, too many for one growth on a workstation. Growth is local, though: a terminal's patch of tissue does not depend on how the rest of the leg is drawn. So the leg grows in two passes:
- **Pass 1** grows the whole leg to about 10 mL patches, one terminal each.
- **Pass 2** grows each patch from its pass-1 terminal down to 7 µm, with the patch's share of the 170 million.

Every radius in both passes then follows Murray's law from the terminals it feeds, so the result is one Murray tree from the femoral artery to the terminals. The patches share no vessels, so they run in batches on a GPU or in parallel on a cluster.
"""

# ╔═╡ 413b0149-ca96-4678-956d-5e67b98f412a
begin
    "Radius (mm) of a terminal arteriole, the size growth must reach: skeletal muscle's last arteriole before the capillaries, fully dilated (section 3)."
    const TERMINAL_RADIUS_MM = 0.007
    "Murray's exponent: the leg's conduits need it (γ = 2.33 starves the beds)."
    const GAMMA = 3.0
    "The terminal count Murray's law gives for the imaged femoral artery and a terminal of radius `TERMINAL_RADIUS_MM`."
    const N_TOTAL = round(Int, (inlets.arteries.r_root / TERMINAL_RADIUS_MM)^GAMMA)
    "Tissue (mL) one pass-1 terminal feeds: one pass-2 patch."
    const PATCH_ML = 10.0
    "Pass-1 terminals: one per patch."
    const N_PASS1 = round(Int, SOFT_ML / PATCH_ML)
end;

# ╔═╡ 47c099f1-608a-42d3-8ea7-a4cc0b76c4bf
let heart_density = 957_045 / 148.1, per_patch = N_TOTAL / N_PASS1
    Markdown.parse(
        @sprintf(
            """
            | | value |
            |---|---|
            | femoral artery radius | %.2f mm |
            | terminal radius | %.0f µm |
            | terminals, Murray (γ = %.0f) | %.1f million |
            | terminals at the heart's density, for comparison | %.1f million |
            | pass 1: patches (%.0f mL each) | %d |
            | pass 2: terminals per patch | about %.0f thousand |
            | pass-1 terminal radius, Murray | %.0f µm |
            """,
            inlets.arteries.r_root,
            1000TERMINAL_RADIUS_MM,
            GAMMA,
            N_TOTAL / 1.0e6,
            heart_density * SOFT_ML / 1.0e6,
            PATCH_ML,
            N_PASS1,
            per_patch / 1000,
            1000 * inlets.arteries.r_root * N_PASS1^(-1 / GAMMA),
        ),
    )
end

# ╔═╡ 13d24d4f-32e5-4c8f-8ff2-9542863323b8
md"""
## 4. Pass 1: The Whole Leg to 10 mL Patches

HCT-CCO's staged growth (the only schedule the package grows with), from the imaged femoral artery and vein with their stubs, into the soft tissue, around the bones. The root radius is fixed: the imaged femoral artery's, and the measured femoral vein's (`ROOT_RADII`); γ = 3. Veins grow the same way from the femoral vein, with the same count, but from their own random stream (`SEEDS`).
"""

# ╔═╡ 1357edb4-59f6-4f23-a522-3d67f4eb503e
"""
The leg's growth preset: `HEART_GROWTH`'s generic values with the leg's own tissue, imaged root and Murray's exponent.
A draft: once the leg is built it belongs in `src/organs.jl` (as `LIMB_GROWTH`), with its sources in
`test/reference/`.
"""
limb_growth(nterm, tissue_ml, root_radius; seed = V.HEART_GROWTH.seed) = V.GrowthSpec(
    V.HEART_GROWTH;
    nterm,
    flow_ml_min = nothing,
    tissue_ml,
    root_caps = [Float64(root_radius)],
    root_fixed = true,
    symmetry_floor = Tuple{Float64, Float64}[],
    gamma = GAMMA,
    seed,
    device = DEVICE,
)

# ╔═╡ 63c202d8-e781-45bf-8efc-340e18870f5b
"""
Grow one tree into `tissue` (a mask on `grid`): `pass` voxels may be crossed but take no terminal. Returns the grown
tree (`V.vtree`) and the seconds it took.
"""
function grow_tree(spec, inlet, tissue, pass, grid, blockers, name)
    dir = mktempdir()
    allow = zeros(UInt16, size(tissue))
    allow[pass] .= V.NOTERM
    domain = V.write_vol(joinpath(dir, "$name.vol"), tissue .| pass, Tuple(round.(Int, grid.origin ./ grid.h)))
    gxl = V.XCAT.write_base_gxl(joinpath(dir, "$name.xml"), inlet)
    seconds = @elapsed engine = V.grow(spec, domain, [gxl]; blockers, territories = allow)
    return (tree = V.vtree(engine), seconds = seconds)
end

# ╔═╡ 63b42cf4-caf6-4c98-ba66-b4a77123e8a2
"""
The random stream each tree's growth draws its candidate terminals from, in pass 1 and pass 2. The arteries keep
`HEART_GROWTH`'s. The veins need their own: with the same stream in the same tissue, they would draw the same candidates
and grow onto the arteries, one vessel inside the other.
"""
const SEEDS = (arteries = V.HEART_GROWTH.seed, veins = V.HEART_GROWTH.seed + 1)

# ╔═╡ b7e2c9a4-5d13-4f6e-9a8b-2c4d6e8f0a13
"""
The radius (mm) each tree grows from. The arteries take XCAT's imaged femoral artery. The veins take the measured supine
femoral vein, 13.6 mm across (Keiler et al. 2019), rather than XCAT's 9.6 mm: every venous radius scales with it, and
from XCAT's the venous pressure at maximal flow leaves about 10% of the capillaries with no pressure across them.
"""
const ROOT_RADII = (arteries = inlets.arteries.r_root, veins = 13.6 / 2)

# ╔═╡ 1e7050c2-417a-4674-9f47-a887894d24bd
pass1 = cached("pass1") do
    map(inlets, SEEDS, ROOT_RADII) do inlet, seed, r_root
        spec = limb_growth(N_PASS1, SOFT_ML, r_root; seed)
        grow_tree(spec, inlet.tree, TISSUE, inlet.pass .| LUMENS, LEG.grid, LEG.blockers, "pass1")
    end
end;

# ╔═╡ 2d36963b-6887-4e12-a57e-99d254c5c9cb
"The terminal segments of a grown tree: the segments without children."
function terminals(tree)
    has_child = falses(length(tree.r))
    foreach(p -> p > 0 && (has_child[p] = true), tree.parent)
    return findall(!, has_child)
end

# ╔═╡ 849eebe7-c721-44ca-b04b-8b75316b95be
let row(name, grown) = (
        t = grown.tree;
        ends = terminals(t);
        r = 1000 .* t.r[ends];
        @sprintf(
            "| %s | %d | %d | %.0f µm (%.0f–%.0f) | %.0f µm | %.1f s |",
            name, length(t.r), length(ends), median(r), quantile(r, 0.05), quantile(r, 0.95),
            1000 * ROOT_RADII[name] * length(ends)^(-1 / GAMMA), grown.seconds,
        )
    ),
        rows = [row(name, grown) for (name, grown) in pairs(pass1)]
    Markdown.parse(
        """
        | tree | segments | terminals | terminal radius, median (5–95%) | Murray: r\\_root · N^(-1/3) | growth |
        |---|---|---|---|---|---|
        $(join(rows, "\n"))

        Each pass-1 terminal is the root of one pass-2 patch. Pass 1 sets the tree's geometry; section 6 replaces its radii with Murray's, from the terminals each vessel feeds.
        """,
    )
end

# ╔═╡ dd3a7c54-c5ed-4139-9211-411c81f568f4
let segments(t) = (t.a, t.b, t.r)
    leg_figure(
        [(segments(pass1.veins.tree), (:royalblue, 0.6), 0.6), (segments(pass1.arteries.tree), (:firebrick, 0.8), 0.6)];
        title = "Pass 1: arteries (red) and veins (blue), one terminal per 10 mL patch",
    )
end

# ╔═╡ ff93d210-7eaa-41ca-979d-233217bc534c
md"""
## 5. The Patches

Each pass-1 terminal's patch is the tissue nearer to it, through the tissue, than to any other terminal: a flood from every terminal at once, one voxel step at a time, which goes around bone instead of through it. It crosses the imaged lumens, since many pass-1 terminals end inside one (the imaged vessels' own ends and stubs), but the lumens belong to no patch. Staged growth gives every terminal an equal share of the tissue, so the patches come out near 10 mL.
"""

# ╔═╡ a1dd86e8-9921-431a-81c1-cc68ddb7311b
"""
    flood(mask, grid, points; keep = mask) -> (label, voxels, box)

Each voxel of `mask` labelled with the nearest of `points` (mm) through the mask (a breadth-first flood from all of them
at once, six neighbours per step), the voxel count per label, and each label's bounding box (`CartesianIndices`). Each
point starts from the nearest free voxel of the mask within one voxel of its own, so a point beside bone still gets its
patch. Only the voxels of `keep` keep their label: the flood crosses the rest.
"""
function flood(mask, grid, points; keep = mask)
    label = zeros(Int32, size(mask))
    queue = CartesianIndex{3}[]
    for (k, p) in enumerate(eachcol(points))
        c0 = CartesianIndex(Tuple(clamp.(round.(Int, (p .- grid.origin) ./ grid.h) .+ 1, 1, size(mask))))
        free = [c0 + d for d in CartesianIndices((-1:1, -1:1, -1:1)) if checkbounds(Bool, mask, c0 + d) && mask[c0 + d] && label[c0 + d] == 0]
        isempty(free) && continue
        c = argmin(c -> sum(abs2, grid.origin .+ (Tuple(c) .- 1) .* grid.h .- p), free)
        label[c] = k
        push!(queue, c)
    end
    head = 1
    while head <= length(queue)
        c = queue[head]
        head += 1
        for s in SIX_STEPS
            n = c + s
            checkbounds(Bool, mask, n) && mask[n] && label[n] == 0 || continue
            label[n] = label[c]
            push!(queue, n)
        end
    end
    label[.!keep] .= 0
    n = size(points, 2)
    voxels = zeros(Int, n)
    lo = fill(CartesianIndex(size(mask)), n)
    hi = fill(CartesianIndex(1, 1, 1), n)
    for c in CartesianIndices(label)
        k = label[c]
        k == 0 && continue
        voxels[k] += 1
        lo[k] = min(lo[k], c)
        hi[k] = max(hi[k], c)
    end
    return (label = label, voxels = voxels, box = [lo[k]:hi[k] for k in 1:n])
end

# ╔═╡ 0a4968ed-9eec-4358-9e0e-9f3ffb03c846
"The arterial patches: each pass-1 arterial terminal's segment and its tissue."
PATCHES = cached("patches") do
    t = pass1.arteries.tree
    ends = terminals(t)
    f = flood(LEG.soft, LEG.grid, t.b[:, ends]; keep = TISSUE)
    (segment = ends, f...)
end;

# ╔═╡ 76695632-e727-470d-bc81-908179e8fe6f
"The venous patches: each pass-1 venous terminal's segment and its tissue, flooded as the arterial ones."
VEIN_PATCHES = cached("vein_patches") do
    t = pass1.veins.tree
    ends = terminals(t)
    f = flood(LEG.soft, LEG.grid, t.b[:, ends]; keep = TISSUE)
    (segment = ends, f...)
end;

# ╔═╡ 15144905-f653-42a0-97f5-67f3963ee428
let ml = PATCHES.voxels ./ 1000
    seeded = count(>(0), PATCHES.voxels)
    Markdown.parse(
        @sprintf(
            """
            | | value |
            |---|---|
            | patches | %d (%d with tissue) |
            | patch volume, median (5–95%%) | %.1f mL (%.1f–%.1f) |
            | tissue in a patch | %.1f%% |
            """,
            length(ml), seeded, median(ml[ml .> 0]), quantile(ml[ml .> 0], 0.05), quantile(ml[ml .> 0], 0.95),
            100 * sum(PATCHES.voxels) / count(TISSUE),
        ),
    )
end

# ╔═╡ dd46fe50-0a2e-41cd-a507-bceab8d0e381
let g = LEG.grid, j = round(Int, (mean(pass1.arteries.tree.b[2, :]) - g.origin[2]) / g.h) + 1
    fig = Mke.Figure(size = (700, 900))
    ax = Mke.Axis(fig[1, 1]; title = "Patches, coronal slice at y = $(round(Int, g.origin[2] + (j - 1) * g.h)) mm", aspect = Mke.DataAspect(), xlabel = "x (mm)", ylabel = "z (mm)")
    slice = Float32.(PATCHES.label[:, j, :])
    colour = Float32.(mod.(slice .* 0.6180339f0, 1.0f0))
    colour[slice .== 0] .= NaN32
    Mke.heatmap!(ax, g.origin[1] .+ (0:(g.dims[1] - 1)) .* g.h, g.origin[3] .+ (0:(g.dims[3] - 1)) .* g.h, colour; colormap = :tab20, nan_color = :white)
    fig
end

# ╔═╡ bf480e2e-8ebb-4030-81fd-6a5e062bdfd2
md"""
## 6. Pass 2: Patches to 7 µm, Stitched into One Tree

**Every radius follows Murray's law from the terminals it feeds.** Pass 1 gives the tree's geometry and its patches, not its radii. Growth balances each terminal's path resistance, so pass-1 terminals that carry the same flow come out up to a third thinner or thicker than each other. A patch grown from a thin one would end in terminals too thin for their flow. So every vessel takes the radius Murray's law gives for the 7 µm terminals below it, r = 7 µm · n^(1/3): each patch's root from the patch's own count, each pass-1 vessel from the counts of all the patches below it. The femoral artery comes out at its imaged radius, because the total count was chosen from it.

**Each patch grows on its own, then joins the tree.** Its tissue is the patch. Its inlet is the pass-1 terminal segment that feeds it, at its Murray radius. Its count is its share of the 170 million by tissue volume, so every terminal carries the same flow. The segment's lumen, and the neighbouring patches' tissue and the imaged lumens within 3 mm, may be crossed but take no terminal, so terminals reach the patch's border and neighbouring patches meet without a seam. The grown tree then replaces that pass-1 segment, so the two passes become one tree.

Here every patch under one pass-1 branch in the calf is grown and stitched in: a contiguous piece of the leg, complete to 7 µm. The checks:
- **terminal radius** near 7 µm;
- **pressure lost** from the femoral artery to each terminal at maximal flow (every vessel dilated, each terminal an equal share of the leg's flow, Poiseuille with the Pries–Secomb in vivo viscosity). The capillaries must keep a pressure above the femoral vein's 11 mmHg.
"""

# ╔═╡ ac61ef6f-c424-4d8f-90fa-41906373d29d
"The segments of a tree (its `parent` vector) in breadth-first order from the root, and each segment's children."
function walk(parent)
    children = [Int[] for _ in eachindex(parent)]
    order = Int[]
    for j in eachindex(parent)
        parent[j] == 0 ? push!(order, j) : push!(children[parent[j]], j)
    end
    for j in order # `order` grows as it is read
        append!(order, children[j])
    end
    return order, children
end

# ╔═╡ ca33ac06-932e-4f4f-a34d-d39d651aefce
"""
    downstream(parent, weight) -> (order, children, n)

The terminals below each segment: a terminal `j` counts `weight(j)` (the final terminals it stands for), every other
segment the sum over its children.
"""
function downstream(parent, weight)
    order, children = walk(parent)
    n = zeros(length(parent))
    for j in reverse(order)
        n[j] = isempty(children[j]) ? weight(j) : sum(n[c] for c in children[j])
    end
    return order, children, n
end

# ╔═╡ 5b216de6-c728-457d-a182-01c990235112
"Each patch's final terminal count: its share of `N_TOTAL` by tissue volume, so every terminal carries the same flow."
const PATCH_TERMINALS = round.(Int, N_TOTAL .* PATCHES.voxels ./ sum(PATCHES.voxels))

# ╔═╡ 6e378f0f-405d-48fe-ac44-36b1621c638f
"""
The pass-1 arteries at their Murray radii: each segment's final terminal count below it (`below`) and its radius
`TERMINAL_RADIUS_MM · below^(1/γ)` (mm). A pass-1 terminal counts its patch's terminals (at least one).
"""
PASS1_MURRAY = let t = pass1.arteries.tree, patch_of = Dict(j => k for (k, j) in enumerate(PATCHES.segment))
    _, _, below = downstream(t.parent, j -> max(PATCH_TERMINALS[patch_of[j]], 1))
    (below = below, radius = TERMINAL_RADIUS_MM .* below .^ (1 / GAMMA))
end;

# ╔═╡ 432afba0-dfb6-4965-9ed4-fbbd231f1e26
begin
    "Each venous patch's final terminal count: its share of `N_TOTAL` by tissue volume, as for the arteries."
    const VEIN_PATCH_TERMINALS = round.(Int, N_TOTAL .* VEIN_PATCHES.voxels ./ sum(VEIN_PATCHES.voxels))

    """
    The pass-1 veins at their Murray radii, as `PASS1_MURRAY` for the arteries. The veins drain as many terminals as the
    arteries feed, so their terminal radius (`terminal`, mm) is the one that brings the femoral vein out at the radius it
    grows from: `ROOT_RADII.veins · N_TOTAL^(-1/γ)`.
    """
    PASS1_MURRAY_VEINS = let t = pass1.veins.tree, patch_of = Dict(j => k for (k, j) in enumerate(VEIN_PATCHES.segment))
        terminal = ROOT_RADII.veins * N_TOTAL^(-1 / GAMMA)
        _, _, below = downstream(t.parent, j -> max(VEIN_PATCH_TERMINALS[patch_of[j]], 1))
        (below = below, radius = terminal .* below .^ (1 / GAMMA), terminal = terminal)
    end
end;

# ╔═╡ 59b8c1c5-e00d-4e74-901c-1ac10ddf5f87
let t = pass1.arteries.tree, ends = PATCHES.segment, root = findfirst(==(0), t.parent)
    ratio = PASS1_MURRAY.radius ./ t.r
    µm(x) = 1000 .* x
    Markdown.parse(
        @sprintf(
            """
            | | value |
            |---|---|
            | femoral artery: imaged / Murray from the counts | %.2f / %.2f mm |
            | pass-1 terminals as grown, median (5–95%%) | %.0f µm (%.0f–%.0f) |
            | pass-1 terminals at their Murray radius, median (5–95%%) | %.0f µm (%.0f–%.0f) |
            | every pass-1 vessel: Murray ÷ grown, median (5–95%%) | %.2f (%.2f–%.2f) |
            """,
            inlets.arteries.r_root, PASS1_MURRAY.radius[root],
            median(µm(t.r[ends])), quantile(µm(t.r[ends]), 0.05), quantile(µm(t.r[ends]), 0.95),
            median(µm(PASS1_MURRAY.radius[ends])), quantile(µm(PASS1_MURRAY.radius[ends]), 0.05), quantile(µm(PASS1_MURRAY.radius[ends]), 0.95),
            median(ratio), quantile(ratio, 0.05), quantile(ratio, 0.95),
        ),
    )
end

# ╔═╡ 857f7ea1-1fc7-4939-b2a0-1b7f334ec058
let t = pass1.veins.tree, ends = VEIN_PATCHES.segment, root = findfirst(==(0), t.parent), ml = VEIN_PATCHES.voxels ./ 1000
    µm(x) = 1000 .* x
    Markdown.parse(
        @sprintf(
            """
            | veins | value |
            |---|---|
            | femoral vein: imaged / grown from / Murray from the counts | %.2f / %.2f / %.2f mm |
            | venous terminal radius, Murray | %.1f µm |
            | patches | %d (%d with tissue), %.1f mL median (5–95%%: %.1f–%.1f) |
            | pass-1 terminals at their Murray radius, median (5–95%%) | %.0f µm (%.0f–%.0f) |
            """,
            inlets.veins.r_root, ROOT_RADII.veins, PASS1_MURRAY_VEINS.radius[root], µm(PASS1_MURRAY_VEINS.terminal),
            length(ml), count(>(0), ml), median(ml[ml .> 0]), quantile(ml[ml .> 0], 0.05), quantile(ml[ml .> 0], 0.95),
            median(µm(PASS1_MURRAY_VEINS.radius[ends])), quantile(µm(PASS1_MURRAY_VEINS.radius[ends]), 0.05),
            quantile(µm(PASS1_MURRAY_VEINS.radius[ends]), 0.95),
        ),
    )
end

# ╔═╡ 641271f3-d904-4d16-893c-848557c02b94
begin
    "Neighbouring tissue (mm) a patch may grow through but not end in, so patches meet without a seam."
    const PATCH_MARGIN_MM = 3.0

    "`mask` grown `steps` voxels through its six face neighbours."
    function dilate(mask, steps)
        out = copy(mask)
        for _ in 1:steps
            grown = copy(out)
            for c in CartesianIndices(out)
                out[c] || continue
                for s in SIX_STEPS
                    n = c + s
                    checkbounds(Bool, out, n) && (grown[n] = true)
                end
            end
            out = grown
        end
        return out
    end
end

# ╔═╡ 670aa3fe-ffe9-42d3-941e-a10c8dcf72dd
"""
Pass 2's inputs for each tree: its pass-1 tree, patches, Murray radii, terminal counts and terminal radius (mm), and the
name its patches are cached under (`<prefix>_<k>`).
"""
SYSTEMS = (
    arteries = (
        tree = pass1.arteries.tree, patches = PATCHES, murray = PASS1_MURRAY, counts = PATCH_TERMINALS,
        terminal_mm = TERMINAL_RADIUS_MM, prefix = "patch",
    ),
    veins = (
        tree = pass1.veins.tree, patches = VEIN_PATCHES, murray = PASS1_MURRAY_VEINS, counts = VEIN_PATCH_TERMINALS,
        terminal_mm = PASS1_MURRAY_VEINS.terminal, prefix = "vein_patch",
    ),
);

# ╔═╡ 4f1209a3-a2aa-4d81-95e8-2f56c5862681
"""
    growth_finds_clear(domain, grid, a, b, clearance) -> Bool

Whether growth finds the segment `a`–`b` (mm) clear in `domain`, a mask on `grid`: the line between the voxels of the
ends' truncated coordinates, where growth puts points, at least `clearance` voxels inside the domain all along. Growth's
own domain and test, as it tests each new vessel and the segment it joins.
"""
function growth_finds_clear(domain, grid, a, b, clearance)
    lower = Tuple(round.(Int, grid.origin ./ grid.h))
    voxel(p) = V.itrunc3(Tuple(Float32.(p ./ grid.h)))
    return mktempdir() do dir
        vol = V.write_vol(joinpath(dir, "check.vol"), domain, lower)
        dom = V.MaskDomain{Float32}(vol; min_dist = clearance, metric = :euclidean)
        V.check_no_intersect(dom, voxel(a), voxel(b))
    end
end

# ╔═╡ eedb0097-f6eb-4984-bcd0-0a994f928cd2
"""
    grow_patch(k; system = :arteries) -> (tree, seconds, nterm, ml, root_radius, segment)

Pass 2 for patch `k` of the arteries or the veins (`system`, a key of `SYSTEMS`): the patch cropped to its box (and the
pass-1 segment that feeds or drains it), grown from that segment at its Murray radius with the patch's share of
`N_TOTAL` terminals. The neighbouring patches' tissue and the imaged lumens within `PATCH_MARGIN_MM`, and the segment's
lumen, may be crossed but take no terminal, so the patch's terminals reach its border and neighbouring patches meet
without a seam.

Bone is taken out of the lumen, unless growth would then find the feeding segment itself blocked (`growth_finds_clear`).
Growth puts a point in the voxel of its truncated coordinates, as upstream HCT-CCO does, which can be a voxel away from
the nearest one, and it tests every new vessel from both ends of the segment it joins. With an end outside the domain,
or the line between the ends shifted into bone (a vein beside the toe bones), no vessel could attach, and the patch would
end with its feeding segment alone. In that case the lumen keeps its bone and both end voxels join it; the fine bone mask
(`LEG.blockers`) still keeps new vessels out of the long bones, and elsewhere only the lumen's few bone voxels may be
crossed.
"""
function grow_patch(k; system = :arteries)
    s = SYSTEMS[system]
    j = s.patches.segment[k]
    a = Float64.(s.tree.a[:, j])
    b = Float64.(s.tree.b[:, j])
    r = Float64(s.murray.radius[j])
    g = LEG.grid
    # the patch's box, grown to hold its feeding segment and the margin
    pad = ceil(Int, r / g.h) + ceil(Int, PATCH_MARGIN_MM / g.h) + 2
    ka = round.(Int, (a .- g.origin) ./ g.h) .+ 1
    lo = max.(Tuple(min(first(s.patches.box[k]), CartesianIndex(Tuple(ka)))) .- pad, 1)
    hi = min.(Tuple(max(last(s.patches.box[k]), CartesianIndex(Tuple(ka)))) .+ pad, g.dims)
    box = CartesianIndices(Tuple(UnitRange.(lo, hi)))
    grid = V.VoxelGrid(g.origin .+ (Tuple(first(box)) .- 1) .* g.h, g.h, size(box))

    tissue = s.patches.label[box] .== k
    soft = LEG.soft[box]
    margin = dilate(tissue, round(Int, PATCH_MARGIN_MM / g.h)) .& soft .& .!tissue
    inlet = V.XCAT.BaseTree("$(s.prefix)_$k", [0, 1], hcat(a, b), [r, r], ["pass1", "pass1"])
    lumen = V.XCAT.vessel_tube(inlet, grid; margin = 1.5)
    pass = (lumen .& .!LEG.bone[box] .| margin) .& .!tissue
    ml = count(tissue) / 1000
    nterm = s.counts[k]
    spec = limb_growth(nterm, ml, r; seed = SEEDS[system])
    # the voxels growth puts the segment's ends in (see the docstring)
    ends = [CartesianIndex(Tuple(trunc.(Int, p ./ grid.h) .- round.(Int, grid.origin ./ grid.h) .+ 1)) for p in (a, b)]
    if !growth_finds_clear(tissue .| pass, grid, a, b, spec.wall_clearance)
        pass .|= lumen .& .!tissue
        foreach(c -> tissue[c] || (pass[c] = true), ends)
    end
    grown = grow_tree(spec, inlet, tissue, pass, grid, LEG.blockers, "$(s.prefix)_$k")
    return (tree = grown.tree, seconds = grown.seconds, nterm = nterm, ml = ml, root_radius = r, segment = j)
end

# ╔═╡ b57241e2-449d-4e78-b18e-44d5b6bc793d
"""
The region grown here: the pass-1 branch with 8 to 16 patches below it whose end lies nearest 30% of the way up the leg
(the calf), and the patches below it.
"""
REGION = let t = pass1.arteries.tree, z = t.b[3, PATCHES.segment]
    _, _, patches = downstream(t.parent, j -> 1)
    height = minimum(z) + 0.3 * (maximum(z) - minimum(z))
    candidates = findall(n -> 8 <= n <= 16, patches)
    branch = candidates[argmin(abs.(t.b[3, candidates] .- height))]
    function under(j)
        while j != 0 && j != branch
            j = t.parent[j]
        end
        return j == branch
    end
    (branch = branch, patches = [k for (k, j) in enumerate(PATCHES.segment) if under(j) && PATCHES.voxels[k] > 0])
end

# ╔═╡ 78ceb802-8b99-480e-9767-bd24ca40f15b
"Pass 2 for every patch of the region, each cached on its own."
pass2 = [
    cached("patch_$k") do
        grow_patch(k)
    end for k in REGION.patches
];

# ╔═╡ 2bdfc5d1-1351-4b5a-b678-90995faf4dd1
"""
    stitch(tree, radius, grown) -> (parent, a, b, r, pass, from)

One tree from the pass-1 `tree` (at its Murray `radius`) and the pass-2 trees `grown`: each grown tree replaces the
pass-1 terminal segment it grew from and hangs from that segment's parent. `pass` is each segment's pass (1 or 2) and
`from` a pass-1 segment's index in `tree` (0 for a grown one).
"""
function stitch(tree, radius, grown)
    replaced = Set(p.segment for p in grown)
    keep = [j for j in eachindex(tree.r) if !(j in replaced)]
    index = zeros(Int, length(tree.r))
    index[keep] = eachindex(keep)
    parent = Int[p == 0 ? 0 : index[p] for p in tree.parent[keep]]
    a = [tree.a[:, keep]]
    b = [tree.b[:, keep]]
    r = [Float32.(radius[keep])]
    pass = [fill(1, length(keep))]
    from = [keep]
    n = length(keep)
    for p in grown
        t = p.tree
        up = index[tree.parent[p.segment]]
        append!(parent, [q == 0 ? up : q + n for q in t.parent])
        push!(a, t.a)
        push!(b, t.b)
        push!(r, t.r)
        push!(pass, fill(2, length(t.r)))
        push!(from, zeros(Int, length(t.r)))
        n += length(t.r)
    end
    return (
        parent = parent, a = reduce(hcat, a), b = reduce(hcat, b), r = reduce(vcat, r),
        pass = reduce(vcat, pass), from = reduce(vcat, from),
    )
end

# ╔═╡ 6ae0879c-c30e-450d-be22-f831d2ac317d
"The leg's arterial tree with the region complete to 7 µm: pass 1 at its Murray radii, the region's patches grafted in."
stitched = stitch(pass1.arteries.tree, PASS1_MURRAY.radius, pass2);

# ╔═╡ e5c89aca-c566-48d4-9ccf-59bbb4ca9e8b
"The final terminals a terminal of `stitched` stands for: one if it was grown in pass 2, its patch's count if not."
terminal_weight = let patch_of = Dict(j => k for (k, j) in enumerate(PATCHES.segment))
    j -> stitched.pass[j] == 2 ? 1.0 : Float64(max(PATCH_TERMINALS[patch_of[stitched.from[j]]], 1))
end

# ╔═╡ 5b657f30-24f2-4481-88c0-1c23e8732280
let fig = Mke.Figure(size = (1500, 860)), t = stitched, g = LEG.grid
    region = findall(==(2), t.pass)
    ends = let has_child = falses(length(t.r))
        foreach(p -> p > 0 && (has_child[p] = true), t.parent)
        filter(j -> !has_child[j], region)
    end
    # the region's box (its terminals) with 10 mm of surrounding leg, and an 8 mm detail at the region's middle
    box_lo = vec(minimum(t.b[:, ends]; dims = 2)) .- 10
    box_hi = vec(maximum(t.b[:, ends]; dims = 2)) .+ 10
    centre = vec(median(t.b[:, ends]; dims = 2))
    detail_lo = centre .- 4
    detail_hi = centre .+ 4
    # the pass-1 path from the femoral artery to the region's branch
    path = Int[findfirst(==(REGION.branch), t.from)]
    while t.parent[path[end]] != 0
        push!(path, t.parent[path[end]])
    end
    within(j, lo, hi) = lo[1] <= t.b[1, j] <= hi[1] && lo[3] <= t.b[3, j] <= hi[3]
    function draw!(ax, segments; colour = nothing, scale = 1.0)
        isempty(segments) && return
        segments = segments[sortperm(t.r[segments])] # thin first, so thick vessels lie on top
        points = Mke.Point2f[]
        for j in segments
            push!(points, Mke.Point2f(t.a[1, j], t.a[3, j]), Mke.Point2f(t.b[1, j], t.b[3, j]))
        end
        return Mke.linesegments!(
            ax, points; color = colour === nothing ? repeat(log10.(1000 .* t.r[segments]); inner = 2) : colour,
            colormap = :viridis, colorrange = (log10(1000TERMINAL_RADIUS_MM), 3.6), linewidth = repeat(max.(0.4, scale .* t.r[segments]); inner = 2),
        )
    end
    xs = [g.origin[d] .+ (0:(g.dims[d] - 1)) .* g.h for d in 1:3]
    silhouette = dropdims(any(LEG.soft; dims = 2); dims = 2)
    outline!(ax, lo, hi) = Mke.lines!(ax, [lo[1], hi[1], hi[1], lo[1], lo[1]], [lo[3], lo[3], hi[3], hi[3], lo[3]]; color = :black, linewidth = 1.5)

    # 1. the leg: pass 1, the path to the region, the region's box
    leg = Mke.Axis(fig[1, 1]; title = "The leg", aspect = Mke.DataAspect(), xlabel = "x (mm)", ylabel = "z (mm)")
    Mke.heatmap!(leg, xs[1], xs[3], silhouette; colormap = [:white, :gray92])
    draw!(leg, findall(==(1), t.pass); colour = :gray55, scale = 2.0)
    draw!(leg, path; colour = :firebrick, scale = 2.0)
    outline!(leg, box_lo, box_hi)

    # 2. the region and 10 mm around it: pass 1 nearby in grey, the region's vessels of 20 µm or more in colour
    zoom = Mke.Axis(fig[1, 2]; title = "The region: $(length(pass2)) patches stitched to pass 1 (vessels ≥ 20 µm)", aspect = Mke.DataAspect(), xlabel = "x (mm)", ylabel = "z (mm)")
    Mke.heatmap!(zoom, xs[1], xs[3], silhouette; colormap = [:white, :gray92])
    draw!(zoom, filter(j -> t.pass[j] == 1 && within(j, box_lo .- 20, box_hi .+ 20), eachindex(t.r)); colour = :gray60, scale = 6.0)
    draw!(zoom, [path; filter(j -> t.r[j] >= 0.02f0, region)]; scale = 6.0)
    outline!(zoom, detail_lo, detail_hi)
    Mke.limits!(zoom, box_lo[1], box_hi[1], box_lo[3], box_hi[3])

    # 3. 8 mm at the region's middle: every vessel, down to the terminals
    detail = Mke.Axis(fig[1, 3]; title = "8 mm detail: every vessel, down to $(round(Int, 1000TERMINAL_RADIUS_MM)) µm", aspect = Mke.DataAspect(), xlabel = "x (mm)", ylabel = "z (mm)")
    draw!(detail, filter(j -> within(j, detail_lo .- 2, detail_hi .+ 2), region); scale = 40.0)
    Mke.limits!(detail, detail_lo[1], detail_hi[1], detail_lo[3], detail_hi[3])
    Mke.Colorbar(fig[1, 4]; colormap = :viridis, limits = (log10(1000TERMINAL_RADIUS_MM), 3.6), label = "log₁₀ radius (µm)")
    Mke.colsize!(fig.layout, 1, Mke.Relative(0.16))
    Mke.colgap!(fig.layout, 1, 60)
    Mke.colgap!(fig.layout, 2, 60)

    # lines from each box's right corners to the next panel's left corners, in figure pixels (they follow the layout),
    # on a scene made after the axes so they are drawn over them
    overlay = Mke.Scene(fig.scene; clear = false, camera = Mke.campixel!)
    function connect!(from, lo, hi, to)
        points = Mke.lift(from.scene.viewport, from.finallimits, to.scene.viewport) do vf, lf, vt
            pixel(x, z) = Mke.Point2f(vf.origin .+ (Mke.Point2f(x, z) .- lf.origin) ./ lf.widths .* vf.widths)
            [pixel(hi[1], hi[3]), Mke.Point2f(vt.origin[1], vt.origin[2] + vt.widths[2]), pixel(hi[1], lo[3]), Mke.Point2f(vt.origin...)]
        end
        return Mke.linesegments!(overlay, points; color = :gray30, linewidth = 1)
    end
    connect!(leg, box_lo, box_hi, zoom)
    connect!(zoom, detail_lo, detail_hi, detail)
    fig
end

# ╔═╡ 5f97677c-390b-41e6-85ab-ceff75a4c140
# Murray's law across the stitched tree: each segment's radius against the final terminals below it, pass 1 and pass 2
# together; the line is r = TERMINAL_RADIUS_MM · n^(1/3)
let fig = Mke.Figure(size = (760, 560)), t = stitched
    _, _, n = downstream(t.parent, terminal_weight)
    ax = Mke.Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "terminals downstream", ylabel = "radius (µm)", title = "One Murray tree from the femoral artery to $(round(Int, 1000TERMINAL_RADIUS_MM)) µm")
    grown = findall(==(2), t.pass)
    keep = grown[1:max(1, length(grown) ÷ 40000):end]
    Mke.scatter!(ax, n[keep], 1000 .* t.r[keep]; markersize = 2, color = (:seagreen, 0.3), label = "pass 2 (the region)")
    first_pass = findall(==(1), t.pass)
    Mke.scatter!(ax, n[first_pass], 1000 .* t.r[first_pass]; markersize = 4, color = (:firebrick, 0.6), label = "pass 1")
    line = 10.0 .^ range(0, log10(N_TOTAL), length = 50)
    Mke.lines!(ax, line, 1000TERMINAL_RADIUS_MM .* line .^ (1 / GAMMA); color = :black, label = "r = $(round(Int, 1000TERMINAL_RADIUS_MM)) µm · n^(1/3)")
    Mke.axislegend(ax; position = :lt)
    fig
end

# ╔═╡ cc1cfb89-1c01-4fe9-ba4a-0abe7041e4ec
let seconds_per_terminal = sum(p.seconds for p in pass2) / sum(p.nterm for p in pass2), segments = 2N_TOTAL
    Markdown.parse(
        @sprintf(
            """
            **The full pass 2, from this region.** %.0f µs of growth per terminal, so the %d patches (%.1f million terminals) take about **%.1f hours** one after another on this machine, or about %.0f minutes spread over 40 GPUs. The arterial tree holds about %.0f million segments, %.1f GB at 40 bytes a segment; the veins as many again.
            """,
            1.0e6seconds_per_terminal, N_PASS1, N_TOTAL / 1.0e6, seconds_per_terminal * N_TOTAL / 3600,
            seconds_per_terminal * N_TOTAL / 60 / 40, segments / 1.0e6, 40 * segments / 1.0e9,
        ),
    )
end

# ╔═╡ ba9579ba-1796-47e6-b072-c4466c0e68cb
md"""
## 7. Flow, Transport and the Phantom

### The leg's own physiology

`LIMB_FLOW` takes every value from measurements in human legs where one exists, and derives the rest from cited measurements. The leg is supine and at rest. Above the ankle its soft tissue is one composite tissue: muscle 0.62, fat 0.34 and skin 0.04 by volume; XCAT's own muscles give the same muscle share there (section 2). The foot below the ankle is mostly fat, tendon and skin, and takes its own maximal flow (`FOOT_MAX_SHARE`, below). Seventeen fields differ from the heart's. These set the oxygen and the flow:

- **O₂ demand** `mvo2` 0.00146 mL O₂/min/g, from the measured resting leg flow (2.6 mL/min/100 mL, Koch Esteves et al. 2023) and the femoral venous saturation;
- **O₂ extraction** `oef` 0.315: femoral venous SO₂ 67.1% against arterial 97.9% (Pearson et al. 2011);
- **venous pressure** `p_out` 11 mmHg (common femoral vein, supine);
- **maximal dilation** `r_total_dilated` 300 mmHg·min·g/mL: 90 mmHg over the PET flow under adenosine (muscle 0.40, fat 0.105 mL/min/g; Heinonen et al. 2010, 2014); the foot's tissue takes `FOOT_MAX_SHARE` of that flow, from its own muscle share;
- **tone** `d_active` 300 µm: skeletal-muscle feed arteries regulate;
- **blood and exchange** `fb` 0.037, `ve` 0.11, `ps` 0.014 mL/min/g (EDTA-derived), `rho` 1.02 g/mL;
- **venous radii** `venous_scale` `:grown`: the veins grow from the measured femoral vein (`ROOT_RADII`), so they keep their grown radii instead of taking the heart's venous-to-arterial volume ratio.

Flow is pressure-driven: growth knows nothing of these values, and the leg's flow comes out of the arterial and venous pressures, the trees and the O₂-driven tone. The resting flow below is what the O₂ balance implies; the solve must find it, not be tuned to it.
"""

# ╔═╡ 51296634-0886-4732-a48f-f1313041b7da
"""
The leg's flow physiology: the fields that differ from the heart's, each with its source. It moves into `src/organs.jl`
as a preset once its sources are in `test/reference/`, as every preset's are. The heart-only fields take the leg's
no-change values: no measured blood-volume rise under adenosine (`fb_dilated = fb`) and a constant PS.
"""
const LIMB_FLOW = V.FlowSpec(
    V.HEART_FLOW;
    p_out = 11.0f0, # mmHg. Common femoral vein, supine rest, catheter (Pearson 2011)
    fb = 0.037f0, # mL/mL. Muscle 3.3 mL/100 g by PET (Raitakari 1995), fat and skin (Leggett 1991), + conduits
    fb_dilated = 0.037f0, # no measured rise under adenosine
    share_art = 0.15f0, # systemic arteries : veins = 13 : 64 (Guyton & Hall), applied to 1 - share_cap
    share_cap = 0.11f0, # muscle capillary volume 0.48 of muscle blood (Gidlöf 1988; Raitakari 1995)
    share_ven = 0.74f0, # the remainder
    cap_length_mm = 0.44f0, # functional capillary path 438 µm (Sarelius 1986)
    cap_diameter_um = 5.3f0, # human quadriceps capillary lumen (Gidlöf 1988)
    ve = 0.11f0, # muscle 0.07, fat 0.14, skin 0.46 by composition (Goodall 2020; Thompson 2005; Levitt 2003)
    ps = 0.014f0, # forearm EDTA PS (Keiding 1988) scaled to iodinated contrast (Levitt 2002)
    ps_slope = 0.0f0,
    rho = 1.02f0, # ICRU-44 muscle, adipose, skin by composition (Hubbell & Seltzer 1995)
    r_total_dilated = 300.0f0, # 90 mmHg / PET flow under adenosine (Heinonen 2010, 2014)
    d_active = 300.0f0, # muscle feed arteries regulate (Williams & Segal 1993; Ives 2012)
    venous_scale = :grown, # the veins grow from the measured femoral vein (ROOT_RADII), their own anatomical scale
    oef = 0.315f0, # supine femoral venous SO2 67.1% vs arterial 97.9% (Pearson 2011)
    mvo2 = 0.00146f0, # 2.6 mL/min/100 mL (Koch Esteves 2023) × O2 content × oef / rho
)

# ╔═╡ 18747ccf-b3d2-4ded-a8dc-b4d7e0637c1a
begin
    """
    The foot's maximal flow per gram, as a share of the composite's. Under adenosine muscle takes 0.40 mL/min/g and fat
    0.105 (Heinonen et al. 2010, 2014, as `r_total_dilated`), so tissue that is a share m muscle takes
    0.40 m + 0.105 (1 − m), the rest at fat's rate. The composite (m = 0.62) is the leg above the ankle; the foot below
    it has XCAT's own muscle share (`MUSCLE`, section 2).
    """
    const FOOT_MAX_SHARE = let below = 1:(ANKLE.slice - 1), flow(m) = 0.40 * m + 0.105 * (1 - m)
        flow(count(view(MUSCLE, :, :, below)) / count(view(LEG.soft, :, :, below))) / flow(0.62)
    end

    "The whole leg's flow at maximal dilation (mL/min): the composite's per gram, with the foot's tissue at `FOOT_MAX_SHARE` of it."
    const Q_MAX = LIMB_FLOW.p_in / LIMB_FLOW.r_total_dilated * LIMB_FLOW.rho * (SOFT_ML - (1 - FOOT_MAX_SHARE) * ANKLE.ml)
end

# ╔═╡ 33b28f6d-f2ca-44a1-b2ca-78e2bbc067c1
begin
    "Each arterial patch's maximal flow per terminal, as a share of the composite's: `FOOT_MAX_SHARE` over its tissue below the ankle, 1 above it."
    const PATCH_MAX_SHARE = let label = PATCHES.label, below = zeros(Int, length(PATCHES.voxels))
        for k in 1:(ANKLE.slice - 1), j in axes(label, 2), i in axes(label, 1)
            l = label[i, j, k]
            l > 0 && (below[l] += 1)
        end
        1 .- (1 - FOOT_MAX_SHARE) .* below ./ max.(PATCHES.voxels, 1)
    end

    "A terminal of `stitched` at maximal flow: its final terminals (`terminal_weight`) times its tissue's share of the composite's flow per terminal."
    flow_weight = let patch_of = Dict(j => k for (k, j) in enumerate(PATCHES.segment))
        j -> stitched.pass[j] == 2 ? (stitched.b[3, j] < ANKLE.z ? FOOT_MAX_SHARE : 1.0) :
            terminal_weight(j) * PATCH_MAX_SHARE[patch_of[stitched.from[j]]]
    end
end

# ╔═╡ 016fdbaa-8dd3-4459-b9a3-f141d21e2b05
"""
    pressure_drop(tree, q, weight) -> (terminals, drop)

Pressure lost (mmHg) from the root to each terminal when the tree carries `q` mL/min, each terminal taking a share in
proportion to `weight(j)`, the final terminals it stands for. Every vessel at its radius (fully dilated), Poiseuille
with the Pries–Secomb in vivo viscosity at the discharge hematocrit.
"""
function pressure_drop(tree, q, weight; phys = LIMB_FLOW)
    order, children, n = downstream(tree.parent, weight)
    unit = q / 60 / n[order[1]] # mL/s per final terminal
    P = zeros(length(tree.r))
    for j in order
        L = Float32(sqrt(sum(abs2, tree.b[:, j] .- tree.a[:, j])))
        r = tree.r[j]
        drop = V.Flow.poiseuille(L, r, phys.plasma_viscosity_mPas * V.Flow.mu_rel(2000 * r, phys.hd)) * unit * n[j]
        P[j] = (tree.parent[j] == 0 ? 0.0 : P[tree.parent[j]]) + drop
    end
    ends = findall(isempty, children)
    return (terminals = ends, drop = P[ends])
end

# ╔═╡ 45930316-b6ba-4992-b30b-ead047135d61
let f = LIMB_FLOW
    rows = map(zip(REGION.patches, pass2)) do (k, p)
        r = 1000 .* p.tree.r[terminals(p.tree)]
        @sprintf(
            "| %d | %.1f mL | %.0f µm | %d | %.1f µm (%.1f–%.1f) | %.0f s |",
            k, p.ml, 1000p.root_radius, p.nterm, median(r), quantile(r, 0.05), quantile(r, 0.95), p.seconds,
        )
    end
    loss = pressure_drop(stitched, Q_MAX, flow_weight)
    grown = loss.drop[stitched.pass[loss.terminals] .== 2]
    Markdown.parse(
        @sprintf(
            """
            | patch | tissue | root radius (Murray) | terminals | terminal radius, median (5–95%%) | growth |
            |---|---|---|---|---|---|
            %s

            The region: %d patches, %.1f mL, %.2f million terminals, %.1f million segments, grown in %.0f s. At maximal flow (%.1f L/min for the whole leg, the foot's tissue at %.2f of the composite's flow per gram), the pressure lost from the femoral artery to the region's terminals is %.0f mmHg (median; 95%%: %.0f), which leaves the capillaries %.0f mmHg of the %.0f mmHg arterial pressure, above the femoral vein's %.0f mmHg.
            """,
            join(rows, "\n"),
            length(pass2), sum(p.ml for p in pass2), sum(p.nterm for p in pass2) / 1.0e6,
            count(==(2), stitched.pass) / 1.0e6, sum(p.seconds for p in pass2),
            Q_MAX / 1000, FOOT_MAX_SHARE, median(grown), quantile(grown, 0.95), f.p_in - median(grown), f.p_in, f.p_out,
        ),
    )
end

# ╔═╡ 04cb85af-10a2-46b7-b415-f5879a214023
let f = LIMB_FLOW, content = 1.34 * f.hb / 100 * f.sao2 # mL O2 per mL blood
    rest = f.mvo2 / (content * f.oef) # mL/min/g
    maximal = f.p_in / f.r_total_dilated
    mass = f.rho * SOFT_ML
    fp(q) = q * (1 - f.hd)
    Markdown.parse(
        @sprintf(
            """
            | | the O₂ balance implies | measured |
            |---|---|---|
            | resting flow | %.4f mL/min/g, %.0f mL/min for this %.0f g leg | 300 ± 30 mL/min (Chiesa 2016), 353 mL/min (Klein 2003) |
            | leg VO₂ | %.1f mL/min | 18–28 mL/min |
            | maximal flow | %.2f mL/min/g (the foot %.2f), %.1f L/min | nutritive, under adenosine |
            | flow reserve | %.1f | — |
            | contrast extraction, rest / maximal | %.2f / %.2f | — |
            | mean blood transit at rest, fb / flow | %.0f s | — |

            The resting transit is long: at rest a leg holds iodine for over a minute, so very little leaves between two volume scans a few seconds apart. Under adenosine it is about %.0f s, and venous outflow matters as it did for the heart at stress.
            """,
            rest, rest * mass, mass, f.mvo2 * mass, maximal, FOOT_MAX_SHARE * maximal, Q_MAX / 1000, maximal / rest,
            1 - exp(-f.ps / fp(rest)), 1 - exp(-f.ps / fp(maximal)),
            60 * f.fb / (rest * f.rho), 60 * f.fb / (maximal * f.rho),
        ),
    )
end

# ╔═╡ 7d3a9e21-4c6b-4f0e-9a52-1b8e6c4d2f90
md"""
### The leg's flow

Each tree has 342 million segments, too many for one `circuit` on a workstation. So the flow is solved patch by patch over the full trees, the way the second pass is grown: a driver next to this notebook (`flow_leg.jl`) solves every patch and keeps the results in `CACHE_DIR/flow`, and the tables below read them.

1. **Venules.** Every venous terminal goes to the arterial patch whose tissue it ends in. A patch's veins are the venous patches' trees cut down to those venules, so no venule drains two patches.
2. **Pairing.** In each patch, `couple_veins` pairs every bed with one of its venules at maximal dilation, at pass 1's pressures for an equal share of the flow. One more solve at 1.25 κ gives the patch's resistance as a + b κ, since κ scales only its microvessels.
3. **The network.** Pass 1, both trees at their Murray radii, with every patch as its resistance, gives each patch's inlet and outlet pressure. κ is set once above the ankle and once for the foot, so that each takes its maximal flow (`r_total_dilated`, the foot at `FOOT_MAX_SHARE`). A patch across the ankle takes a κ in between.
4. **Each patch again**, at the network's pressures: maximal dilation, then rest with the same κ and a tone constant set by the patch's own O₂ balance.

On the arterial side the patches share no vessels, so the split loses only the straight-line fit of each patch's resistance in κ: the patches' final flows add up to the network's. The venous side is approximate. Each patch drains at one outlet pressure, the mean over its venous patches, and a venous trunk that several patches share carries only the solved patch's flow. The veins take little of the pressure drop inside a patch (below), so the error is small. The demand field and the rest tone are also set per patch, not once for the leg.
"""

# ╔═╡ a4f1c7e8-2b93-4d5a-8e61-9c0b3f7a5d24
let dir = joinpath(CACHE_DIR, "flow"), f = LIMB_FLOW
    ks = filter(k -> isfile(joinpath(dir, "final_$k.jls")), collect(eachindex(PATCHES.voxels)))
    if isempty(ks) || !isfile(joinpath(dir, "network.jls"))
        md"_Not solved here yet: the flow needs every patch of the second pass, arteries and veins (section 9)._"
    else
        F = Dict(k => deserialize(joinpath(dir, "final_$k.jls")) for k in ks)
        N = deserialize(joinpath(dir, "network.jls"))
        g(k) = F[k].ml * f.rho # grams
        content = 1.34 * f.hb / 100 * f.sao2 # mL O₂ per mL blood
        mass = sum(g, ks)
        qs, qr = 60 * sum(k -> F[k].stress.q, ks), 60 * sum(k -> F[k].rest.q, ks) # mL/min
        so2 = sum(k -> F[k].rest.q * F[k].rest.so2, ks) / sum(k -> F[k].rest.q, ks)
        vo2 = sum(k -> 60 * F[k].rest.q * 1.34 * f.hb / 100 * (f.sao2 - F[k].rest.so2), ks)
        nterm = sum(k -> F[k].nterm, ks)
        starved = sum(k -> F[k].rest.starved * F[k].nterm, ks) # beds below 5% SO₂
        act = [F[k].rest.act for k in ks]
        share(i) = 100 * median([F[k].stress.shares[i] for k in ks])
        z = pass1.arteries.tree.b[3, PATCHES.segment]
        height = (z .- minimum(z)) ./ (maximum(z) - minimum(z))
        foot(k) = PATCH_MAX_SHARE[k] < (1 + FOOT_MAX_SHARE) / 2 # most of its tissue below the ankle
        regions = (
            ("foot, below the ankle", foot),
            ("calf, up to 40% of the height", k -> !foot(k) && height[k] < 0.4),
            ("knee, 40–70%", k -> !foot(k) && 0.4 <= height[k] < 0.7),
            ("thigh, 70–100%", k -> !foot(k) && height[k] >= 0.7),
        )
        rows = map(regions) do (name, inside)
            r = filter(inside, ks)
            @sprintf(
                "| %s | %d | %.0f mL | %.3f | %.4f | %.0f mmHg |",
                name, length(r), sum(k -> F[k].ml, r), median([60 * F[k].stress.q / g(k) for k in r]),
                median([60 * F[k].rest.q / g(k) for k in r]), median(N.pin[r]),
            )
        end
        Markdown.parse(
            @sprintf(
                """
                | | solved | from the physiology, or measured |
                |---|---|---|
                | maximal flow | %.0f mL/min, %.3f mL/min/g (the network: %.0f mL/min) | `r_total_dilated`, the foot at `FOOT_MAX_SHARE` |
                | resting flow | %.0f mL/min, %.4f mL/min/g | %.4f mL/min/g by the O₂ balance; 300 ± 30 (Chiesa 2016), 353 mL/min (Klein 2003) |
                | femoral venous SO₂ at rest | %.3f | %.3f |
                | leg VO₂ | %.1f mL/min | 18–28 mL/min |
                | flow reserve | %.1f | — |
                | reversed beds, maximal / rest | %d / %d of %.0f million | none |
                | tone at rest: mean activation per patch | %.2f (5–95%%: %.2f–%.2f); %d patches with a vessel at full tone | not saturated |
                | beds below 5%% SO₂ at rest | %d (%.4f%%) | — |
                | pressure drop inside a patch at maximal flow | arteries %.0f%%, capillaries %.0f%%, veins %.0f%% (medians) | — |
                | κ, above the ankle / the foot | %.3f / %.3f | — |

                | region | patches | tissue | maximal flow (mL/min/g) | rest (mL/min/g) | inlet pressure at maximal flow |
                |---|---|---|---|---|---|
                %s

                %d of %d patches solved, %.0f of %.0f mL of tissue. The regions give the median over their patches; the height is the patch's pass-1 terminal, between the lowest and the highest.
                """,
                qs, qs / mass, 60 * sum(N.Q), qr, qr / mass, f.mvo2 / (content * f.oef),
                so2, f.sao2 * (1 - f.oef), vo2, qs / qr,
                sum(k -> F[k].stress.reversed, ks), sum(k -> F[k].rest.reversed, ks), nterm / 1.0e6,
                median(act), quantile(act, 0.05), quantile(act, 0.95), count(k -> F[k].rest.max_tone > 0, ks),
                starved, 100 * starved / nterm, share(1), share(2), share(3), N.kappa_leg, N.kappa_foot,
                join(rows, "\n"), length(ks), count(>(0), PATCHES.voxels), mass / f.rho, SOFT_ML,
            ),
        )
    end
end

# ╔═╡ 3c9f5b7a-1d2e-4f60-8a4b-6e7d9c0b1a25
md"""
### The leg's contrast transport

The package's `transport` runs the injection through the whole-body model with a grown heart or lungs in that organ's place (`ORGAN_SLOTS`). The leg takes the place of the body's right leg. A driver next to this notebook (`transport_leg.jl`) runs the same loop with the leg there, at rest, calling the package's body model and circuit as they are. It keeps the results in `CACHE_DIR/transport_rest`, and the tables below read them.

1. **The body.** `Body.human_circulation` with the right leg's share of the cardiac output set to the grown leg's resting flow, the other regions keeping their proportions, so the body and the leg carry the same blood, as `coronary_share` does for the heart. The leg's input is the plasma leaving the abdominal aorta. Its femoral venous outflow replaces the body's right-leg beds (`Body.simulate`'s `detailed`), repeated until that outflow changes by at most 10⁻³, as `transport` iterates.
2. **Each patch** at rest: its circuit solved again as for the flow. Pass 1 at the resting flows gives the transit to each patch's inlet, and from each venous patch's root to the femoral vein.
3. **Shared veins.** A vein that drains several arterial patches is in each of their circuits, with only that patch's flow, so each circuit counts its whole volume, at too long a transit. Here its transit is its plasma volume over the plasma flow of every bed it drains, summed over the patches, so the leg holds each vein once.
4. **The leg's response.** Transport is linear in the input. So every bed's arrival at its patch's inlet, at its own inlet and outlet, at its venous patch's root and at the femoral vein is binned on a 0.02 s grid into step responses, and the iodine in each compartment is the input convolved with them. On single patches this matches the package's `frame_masses` and `sinus_outflow` (below). The beds differ by the package's own imbalance between its bed masses and their inflow less outflow. Its input tables (`InputCurve`) are computed in single precision, which loses digits at the leg's slow back-flux from the interstitium (k·dt below about 10⁻³). The step responses here are in double precision and hold each bed's iodine as exactly its inflow less its outflow.

Maximal dilation waits on its adenosine route (*Transport*, below).
"""

# ╔═╡ 8b2e4d6f-5a71-4c3e-9f08-2d4b6a8c0e13
let path = joinpath(CACHE_DIR, "transport_rest", "result.jls")
    if !isfile(path)
        md"_Not run here yet: the transport needs the flow of every patch (above)._"
    else
        R = deserialize(path)
        b = R.body
        pf = 1 - LIMB_FLOW.hd
        timing(t, y) = (i = argmax(y); j = findfirst(>=(0.1 * y[i]), y); (t[j], y[i], t[i]))
        at(t, x) = findfirst(s -> abs(s - x) < 1.0e-3, t)
        a = timing(R.grid, pf .* R.c) # blood
        v = timing(R.tc, pf .* R.cs)
        p = R.vp
        nr = length(R.regions)
        rows_t = map((15, 20, 30, 60, 90, 120)) do tq
            i = at(R.frames, tq)
            @sprintf(
                "| %d s | %.0f | %.1f | %.1f | %.0f | %.1f | %.1f |", tq, R.total[i], R.art1[i], sum(r -> R.art2[r][i], 1:nr),
                sum(r -> R.beds[r][i], 1:nr), sum(r -> R.ven2[r][i], 1:nr), R.ven1[i],
            )
        end
        rows_r = map(1:nr) do r
            conc = (R.art2[r] .+ R.beds[r] .+ R.ven2[r]) ./ R.ml[r] # mg per mL of tissue
            c(tq) = conc[at(R.frames, tq)]
            @sprintf(
                "| %s | %.0f mL | %.3f | %.3f | %.3f | %.3f | %.3f | %.0f%% |", R.regions[r], R.ml[r], c(20), c(30), c(60), c(90), c(120),
                100 * R.beds[r][end] / (conc[end] * R.ml[r]),
            )
        end
        worst(f) = maximum(getfield(x, f) for x in R.compare)
        Markdown.parse(
            @sprintf(
                """
                | | at rest |
                |---|---|
                | the body's right leg | %.0f mL/min, %.2f%% of the cardiac output's %.0f mL/min (the body's own right leg: %.0f) |
                | feedback passes | %d (the leg's outflow then changes by %.0e, weighted by its share) |
                | the leg's input, the abdominal aorta's outflow | arrives at %.1f s (10%% of its peak), peaks at %.2f mg/mL of blood at %.1f s |
                | femoral venous outflow | arrives at %.1f s, peaks at %.3f mg/mL of blood at %.1f s |
                | the leg's plasma | %.0f mL: pass-1 arteries %.1f, patch arteries %.1f, beds %.0f, patch veins %.0f, pass-1 veins %.1f. The patches' own circuits would hold %.0f mL in their veins |
                | mean transit, plasma volume ÷ plasma flow | %.0f s (the beds' flow-weighted arrival at the femoral vein: %.0f s) |
                | iodine at %.0f s | %.0f mg in the leg, of %.0f mg injected. The body outside it holds %.0f mg and has excreted %.0f, which leaves %.0f; with its own lumped leg in place, the body leaves %.0f mg unaccounted |
                | three patches alone, against the package | arteries %.0e, beds %.0e, veins %.0e, venous outflow %.0e (the largest difference over each one's peak). The package's own compartments against its inflow less its outflow: %.0e |

                The iodine in the leg (mg):

                | time | leg | pass-1 arteries | patch arteries | beds | patch veins | pass-1 veins |
                |---|---|---|---|---|---|---|
                %s

                The patches' iodine (their arteries, beds and veins) per mL of tissue (mg/mL), by region (as for the flow):

                | region | tissue | 20 s | 30 s | 60 s | 90 s | 120 s | in the beds at 120 s |
                |---|---|---|---|---|---|---|---|
                %s
                """,
                60 * b.share * b.co, 100 * b.share, 60 * b.co, 60 * b.lumped, b.passes, b.change, a[1], a[2], a[3], v[1], v[2], v[3],
                p.art1 + p.art2 + p.beds + p.ven2 + p.ven1, p.art1, p.art2, p.beds, p.ven2, p.ven1, p.own,
                (p.art1 + p.art2 + p.beds + p.ven2 + p.ven1) / R.fp, R.mean_arrival,
                R.frames[end], R.total[end], b.injected, b.held, b.excreted, b.injected - b.held - b.excreted, b.gap0,
                worst(:arteries), worst(:tissue), worst(:veins), worst(:outflow), worst(:package_balance),
                join(rows_t, "\n"), join(rows_r, "\n"),
            ),
        )
    end
end

# ╔═╡ 847c2060-2a86-4308-8ff6-8a4ba65aae63
md"""
### A narrowed superficial femoral artery

The leg's counterpart of the heart's LAD stenosis: the superficial femoral artery narrowed over 20 mm in the distal adductor canal, below the deep femoral artery, where superficial femoral stenoses are most common. Its diameter is narrowed from 0 to 95%, at rest and at maximal dilation, and every part of the leg below it, down to the foot, is compared with the healthy leg. A driver next to this notebook (`stenosis_leg.jl`) solves it on the flow's patches:

1. **Each patch's curve.** Its flow against its inlet pressure: at rest with its healthy tone constant, so its arterioles dilate as the pressure falls; at maximal dilation with its κ.
2. **Pass 1 on the curves.** The tree's flows and pressures, exact on the tree, with the narrowing's loss from the package circuit's own `stenosis` (Young–Tsai); the throat does not dilate.
3. **Exact patches.** At 50, 70 and 90%, every patch whose pressures moved is solved again in full, until nothing moves.

`stenosis_figs.jl` draws the figures below from its result.
"""

# ╔═╡ 7ac4cb73-9e0d-424d-bbe2-455732481b4a
let dir = joinpath(CACHE_DIR, "stenosis", "figures")
    paths = joinpath.(dir, ["flow_vs_stenosis.png", "flow_by_part.png", "pressure_along_leg.png", "flow_maps.png"])
    if !all(isfile, paths)
        md"_Not run here yet: `stenosis_leg.jl` needs the flow of every patch (above); then `stenosis_figs.jl` draws its figures._"
    else
        img(p) = LocalResource(p, :style => "max-width: 100%; height: auto;")
        md"""
        $(img(paths[1]))

        $(img(paths[2]))

        $(img(paths[3]))

        $(img(paths[4]))

        Every number behind them: `CACHE_DIR/stenosis/summary.txt`, and the tables in `CACHE_DIR/stenosis/figures/figures.html`.
        """
    end
end

# ╔═╡ 15efba44-39d9-490c-987b-4c21692c5b15
md"""
### Iodine in the whole body

Where the injected iodine is, organ by organ, from the injection to 2 minutes, with the grown leg in the body's right-leg place, at rest. A driver next to this notebook (`whole_body.jl`) runs the transport's last pass through the body model again, with a step-for-step copy of `Body.simulate` that also records every segment's iodine (its blood and the space around its cells) every 0.1 s and the iodine the kidneys clear, and checks the copy against `Body.simulate`. The grown leg's iodine comes from its step responses. The heart and the lungs are the body model's own lumped ones. `whole_body_figs.jl` draws the figures below.
"""

# ╔═╡ 4da54df2-91aa-4eb5-a84b-305167963b31
let dir = joinpath(CACHE_DIR, "whole_body", "figures")
    paths = joinpath.(dir, ["blood_flow_by_organ.png", "iodine_by_organ.png", "where_the_iodine_is.png"])
    if !all(isfile, paths)
        md"_Not run here yet: `whole_body.jl` needs the transport (above); then `whole_body_figs.jl` draws its figures._"
    else
        img(p) = LocalResource(p, :style => "max-width: 100%; height: auto;")
        md"""
        $(img(paths[1]))

        $(img(paths[2]))

        $(img(paths[3]))

        Every number behind them: `CACHE_DIR/whole_body/summary.txt`, and the tables in `CACHE_DIR/whole_body/figures/figures.html`.
        """
    end
end

# ╔═╡ ccb1198f-120e-49a1-aa97-0de84552442e
md"""
### Transport, compartments and the phantom

These follow the heart and lungs step for step, on the full pass-2 trees.

1. **Veins, pass 2.** The same, from the femoral vein's pass-1 terminals: `grow_patch(k; system = :veins)` on the venous patches (`VEIN_PATCHES`). The veins drain the same count, so the femoral vein's radius (`ROOT_RADII.veins`) sets their terminal radius.
2. **Flow.** Solved patch by patch above (*The leg's flow*), at rest and at maximal dilation.
3. **Transport.** The injection (`INJECTION`, below) through the whole-body model, with the leg in the body's right-leg place: at rest above (*The leg's contrast transport*). Maximal dilation needs its adenosine route first. The maximal flows behind `r_total_dilated` were measured with adenosine infused into the femoral artery, which dilates the leg alone; the heart's stress protocol infuses it into a vein, which dilates every bed and raises the cardiac output 1.71× (`Body.adenosine()`).
4. **The AIF.** A tracking point in the imaged **femoral artery**, at the top of the scanned slab. The package gives a tracking point in a vessel its own label at the segment's blood concentration, so its curve is the true arterial input. Put the AIF ROI inside that one segment's lumen: the other imaged segments are painted at their mass over XCAT's lumen, which differs from the grown radius.
5. **Compartments.** The muscles (or thigh, calf and foot), each also grouped with the conduit vessels that feed and drain it, as the heart's "myocardium + epicardial vessels" and the lungs' "lungs + hilar vessels": the ideal first-pass compartment counts the iodine in transit.
6. **The phantom.** XCAT has no 0.2 mm label volume of the leg, and the whole leg at 0.2 mm is 8 billion voxels. Build the phantom for one scanned slab (the 160 mm beam: the thigh, or the calf), voxelizing XCAT's leg surfaces there.
"""

# ╔═╡ f6355eaf-f55a-4d71-99f2-f68b2dd5e88d
"""
The injection, the same as the heart's and lungs' (Hadjiabdolhamid et al. 2024, Shojazadeh et al. 2026): 0.5 mL/kg of
contrast at 5 mL/s, then 0.5 mL/kg of a 30:70 contrast:saline chaser at the same rate, for the body model's 70 kg.
"""
const INJECTION = let kg = 70, conc = 350
    V.Body.Injection(; phases = [(0.5 * kg, 5, conc), (0.5 * kg, 5, 0.3 * conc)])
end

# ╔═╡ c6d2e2bb-6ace-4a49-b0de-ebfe038aab20
md"""
## 8. Design Choices Carried Over from the Heart and Lungs

**One source of truth**
- One code path for every organ: grow → flow → transport → phantom. Organs differ only in their spec values (`GrowthSpec`, `FlowSpec`, `TransportSpec`, `PhantomSpec`). The leg's presets live in this notebook until they join the others in `src/organs.jl`.
- Every organ has its own physiology. The leg takes none of the heart's values where a leg measurement exists: γ = 3, its own O₂ demand, extraction, pressures and blood volumes.
- Flow is pressure-driven. Ground truth comes out of the trees, the pressures and the tone; measured flow is only a comparator. Growth and flow are never tuned to make a perfusion measurement agree.
- A tree too large for one circuit is solved patch by patch: each patch's resistance joins the pass-1 network, κ is set per region from its maximal flow, and each patch is solved again at the network's pressures. Its transport runs on the same patches: a vein that several patches share takes the flow of every bed it drains, and the leg's response to its input is a set of step responses.

**Growth**
- HCT-CCO's staged schedule only.
- Imaged vessels get penetrating stubs, so growth starts along their whole length (the heart's and lungs' imaged trees had too few ends without them).
- The imaged lumens are crossable but take no terminal; blockers keep vessels out of bone (the lungs' fissures, the heart's chambers).
- Terminal size follows the count, r\_root · N^(-1/γ). Reach real terminal arterioles by growing more terminals (two passes), never by lowering γ.

**Anatomy and timing**
- Check XCAT against measurements. Its aorta is about three standard deviations narrow, so the aorta's timing was taken from the body model (`V.XCAT.scale_tree`). Its main pulmonary artery is large. Its femoral vein is narrower than measured, so the veins grow from the measured one (`ROOT_RADII`).
- The AIF must be in time with the compartment's inlet. The ideal compartment counts the iodine still in transit in its conduit vessels (epicardial, hilar; here the femoral and popliteal branches).
- V2 about 1 s before the arterial peak loses the least iodine out of the veins. At rest the leg barely loses any; under adenosine it will.
- One injection for every organ: the lab's 0.5 mL/kg protocol.

**Painting**
- An imaged vessel's radius is set by its flow, not by XCAT's lumen. Its iodine is spread over XCAT's lumen, so its painted concentration can differ from its blood's.
- A tracking point's segment is painted at its true blood concentration. Keep the CT AIF ROI inside that segment.

**Notebooks**
- Loops over notebook globals are slow; they run inside functions.
- Every heavy step is cached.
- Many segments are drawn with one `linesegments!`.
"""

# ╔═╡ 333ce447-64ae-4ec5-898b-8a8d86427091
md"""
## 9. Completing the Build

The scripts next to this notebook (`docs/limb/`) run the long steps; its README gives the order, the commands and the times.

1. **Pass 2 for every arterial patch** (done, `grow_stage2.jl`), as parallel jobs. Each patch is independent and cached on its own; the patches are then grafted onto pass 1 (`view_stage2.jl` draws and checks them).
2. **Veins** (done, `grow_stage2.jl veins`): the same two passes from the femoral vein (pass 2: `grow_patch(k; system = :veins)`).
3. **Flow** (done, `flow_leg.jl`): each patch, then the leg, at rest and at maximal dilation, against the checks in section 7 (*The leg's flow*).
4. **The leg's presets in the package**: `LIMB_GROWTH` and `LIMB_FLOW` in `src/organs.jl`, their sources (and those of `TERMINAL_RADIUS_MM`) in `test/reference/`, with tests.
5. **Transport** (done at rest, `transport_leg.jl`, section 7): the right leg's share of the cardiac output from the grown leg's flow, and the injection through the body and the leg. Maximal dilation once its adenosine route is chosen; then the femoral AIF.
6. **A narrowed superficial femoral artery** (done at rest and at maximal dilation, `stenosis_leg.jl`) and **the iodine in the whole body** (done at rest, `whole_body.jl`), section 7.
7. **Compartments and a phantom slab**, then idealized first-pass perfusion against ground truth, as for the heart and lungs.
"""


# ╔═╡ Cell order:
# ╟─2d88882b-057d-4f1d-a810-03f14c411fae
# ╟─4a662e6f-8a3c-4219-a3a6-03a7c848b45e
# ╠═329820de-86c9-4483-8daf-4cdb31dad119
# ╠═b3d8dbc1-fce8-44ba-b4f4-062be2f1769e
# ╠═7cd96038-c409-4a12-8367-f8f1c28c7ae8
# ╠═ed211fb3-e82c-4a06-ab74-8b8cf919e773
# ╠═63e77379-a0ca-4b7f-9bda-a244fb1896be
# ╠═939015ba-d0c0-4411-9223-131ae2d336e1
# ╠═c91bff6a-765c-4bbd-93bd-b81f892a174e
# ╠═b24bcc05-e7c8-4d75-884e-876a5bf0b4b5
# ╠═dcdca1e1-7737-4f0c-b019-17179dcd9530
# ╠═4e477957-3699-4072-99ee-464609f3d3f1
# ╠═c45cbbd6-6d31-4d1d-9d54-fc5dfb7870fb
# ╠═77a6aaf7-fa17-470e-9e01-34bc52fbfd4f
# ╠═8393831b-051f-4982-a98f-d28277891340
# ╟─fc85e6c1-ab4a-43a8-aa91-91d6f74dba97
# ╠═f7984dcd-8660-4bb2-b99e-433d60b08960
# ╠═0928ebc1-1d0b-4b61-a361-637c87839d35
# ╠═5776c7ba-8194-47d8-a28e-046ef1205abb
# ╟─2c4e27bb-61d1-4d72-9dfa-8e496e27c231
# ╟─ff0ae1c4-2c23-4fee-b563-57ddb6d3e5ec
# ╠═d57038e3-65cf-49fc-b365-1d3888b8b148
# ╠═82d39e1f-f872-4449-838e-81678f93d077
# ╠═bbfd5e8b-3f19-4db3-969f-eb73cb98ae4d
# ╠═41d670f3-4243-48ce-b340-b0a6b8115b57
# ╟─1b2a4544-4b53-4564-8039-0c4feac36a0e
# ╠═ca00e0f2-712a-4f9e-8b82-2ba6543e3f18
# ╠═2c77ebbd-5cc6-4081-aacf-1f6e10e9a78a
# ╠═fa99ca18-b2a3-4a33-8cc1-54c84219aced
# ╠═6536e61a-f90e-4a09-973b-0f3f15bf2542
# ╠═4fb8cda2-0411-4a9e-a2f5-f4156b25a698
# ╠═5c1e8f3a-7b42-4d96-a0e5-3f9b2d6c8e17
# ╠═d3bfa261-2daf-42e5-9bfe-4bc1d4a72f3c
# ╟─6e3760ce-8bde-4094-b3ba-d19d5c6cb0cb
# ╠═6193958a-729a-4632-b491-4e2c0ca5edf1
# ╠═2cc35552-467a-418b-9ebd-cd4a927e5b61
# ╟─826f47f1-d661-4ae1-83e7-a6a0fe1c1e62
# ╟─2a5c666d-73d6-49aa-b3c2-285b61e90e26
# ╠═413b0149-ca96-4678-956d-5e67b98f412a
# ╟─47c099f1-608a-42d3-8ea7-a4cc0b76c4bf
# ╟─13d24d4f-32e5-4c8f-8ff2-9542863323b8
# ╠═1357edb4-59f6-4f23-a522-3d67f4eb503e
# ╠═63c202d8-e781-45bf-8efc-340e18870f5b
# ╠═63b42cf4-caf6-4c98-ba66-b4a77123e8a2
# ╠═b7e2c9a4-5d13-4f6e-9a8b-2c4d6e8f0a13
# ╠═1e7050c2-417a-4674-9f47-a887894d24bd
# ╠═2d36963b-6887-4e12-a57e-99d254c5c9cb
# ╟─849eebe7-c721-44ca-b04b-8b75316b95be
# ╟─dd3a7c54-c5ed-4139-9211-411c81f568f4
# ╟─ff93d210-7eaa-41ca-979d-233217bc534c
# ╠═a1dd86e8-9921-431a-81c1-cc68ddb7311b
# ╠═0a4968ed-9eec-4358-9e0e-9f3ffb03c846
# ╠═76695632-e727-470d-bc81-908179e8fe6f
# ╟─15144905-f653-42a0-97f5-67f3963ee428
# ╟─dd46fe50-0a2e-41cd-a507-bceab8d0e381
# ╟─bf480e2e-8ebb-4030-81fd-6a5e062bdfd2
# ╠═ac61ef6f-c424-4d8f-90fa-41906373d29d
# ╠═ca33ac06-932e-4f4f-a34d-d39d651aefce
# ╠═5b216de6-c728-457d-a182-01c990235112
# ╠═6e378f0f-405d-48fe-ac44-36b1621c638f
# ╟─59b8c1c5-e00d-4e74-901c-1ac10ddf5f87
# ╠═432afba0-dfb6-4965-9ed4-fbbd231f1e26
# ╟─857f7ea1-1fc7-4939-b2a0-1b7f334ec058
# ╠═670aa3fe-ffe9-42d3-941e-a10c8dcf72dd
# ╠═4f1209a3-a2aa-4d81-95e8-2f56c5862681
# ╠═eedb0097-f6eb-4984-bcd0-0a994f928cd2
# ╠═641271f3-d904-4d16-893c-848557c02b94
# ╠═b57241e2-449d-4e78-b18e-44d5b6bc793d
# ╠═78ceb802-8b99-480e-9767-bd24ca40f15b
# ╠═2bdfc5d1-1351-4b5a-b678-90995faf4dd1
# ╠═6ae0879c-c30e-450d-be22-f831d2ac317d
# ╠═e5c89aca-c566-48d4-9ccf-59bbb4ca9e8b
# ╠═016fdbaa-8dd3-4459-b9a3-f141d21e2b05
# ╟─45930316-b6ba-4992-b30b-ead047135d61
# ╟─5b657f30-24f2-4481-88c0-1c23e8732280
# ╟─5f97677c-390b-41e6-85ab-ceff75a4c140
# ╟─cc1cfb89-1c01-4fe9-ba4a-0abe7041e4ec
# ╟─ba9579ba-1796-47e6-b072-c4466c0e68cb
# ╠═51296634-0886-4732-a48f-f1313041b7da
# ╠═18747ccf-b3d2-4ded-a8dc-b4d7e0637c1a
# ╠═33b28f6d-f2ca-44a1-b2ca-78e2bbc067c1
# ╟─04cb85af-10a2-46b7-b415-f5879a214023
# ╟─7d3a9e21-4c6b-4f0e-9a52-1b8e6c4d2f90
# ╟─a4f1c7e8-2b93-4d5a-8e61-9c0b3f7a5d24
# ╟─3c9f5b7a-1d2e-4f60-8a4b-6e7d9c0b1a25
# ╟─8b2e4d6f-5a71-4c3e-9f08-2d4b6a8c0e13
# ╟─847c2060-2a86-4308-8ff6-8a4ba65aae63
# ╟─7ac4cb73-9e0d-424d-bbe2-455732481b4a
# ╟─15efba44-39d9-490c-987b-4c21692c5b15
# ╟─4da54df2-91aa-4eb5-a84b-305167963b31
# ╟─ccb1198f-120e-49a1-aa97-0de84552442e
# ╠═f6355eaf-f55a-4d71-99f2-f68b2dd5e88d
# ╟─c6d2e2bb-6ace-4a49-b0de-ebfe038aab20
# ╟─333ce447-64ae-4ec5-898b-8a8d86427091
