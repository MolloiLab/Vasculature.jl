# Slide figures and tables of the leg with a narrowed superficial femoral artery, from stenosis_leg.jl's result.jls, drawn by
# stenosis_figs_plot.jl:
#   flow_by_part.png        each part's flow at 50, 70 and 90%, at rest and at maximal dilation
#   flow_vs_stenosis.png    the tissue fed through the narrowing and above it, from 0 to 95%
#   pressure_along_leg.png  the main artery's pressure from the top of the leg to the foot
#   flow_maps.png           every tissue patch's flow as a share of healthy
#   figures.html            the four, and the numbers as tables
#
#     julia stenosis_figs.jl
#
# Numbers only. The figures are made from XCAT: they stay in LIMB_CACHE_DIR\stenosis\figures (lab storage).

haskey(ENV, "LIMB_CACHE_DIR") || error("set LIMB_CACHE_DIR first")
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
import CairoMakie as Mke
using Serialization, Statistics, Printf
include(joinpath(@__DIR__, "stenosis_figs_plot.jl"))

const DIR = joinpath(ENV["LIMB_CACHE_DIR"], "stenosis")
const OUT = mkpath(joinpath(DIR, "figures"))

R = deserialize(joinpath(DIR, "result.jls"))
dets = R.detail_severities
states = ((:rest, "Rest"), (:stress, "Maximal dilation"))

println(fig_flow_by_part(R, joinpath(OUT, "flow_by_part.png")))
println(fig_flow_vs_stenosis(R, joinpath(OUT, "flow_vs_stenosis.png")))
println(fig_pressure(R, joinpath(OUT, "pressure_along_leg.png")))
println(fig_maps(R, joinpath(OUT, "flow_maps.png")))

# the tables
cell(q, q0) = @sprintf("%.1f (%.0f%%)", 60q, 100q / q0)
flow_rows = map(eachindex(R.part_names)) do p
    row = [R.part_names[p], @sprintf("%.0f", R.nets[(:rest, 0.0)].part_ml[p])]
    for (state, _) in states
        h = R.nets[(state, 0.0)].part_q[p]
        push!(row, @sprintf("%.1f", 60h))
        append!(row, [cell(R.nets[(state, ds)].part_q[p], h) for ds in dets])
    end
    row
end
leg_row = let row = ["Whole leg", @sprintf("%.0f", sum(R.nets[(:rest, 0.0)].part_ml))]
    for (state, _) in states
        h = R.nets[(state, 0.0)].leg_q
        push!(row, @sprintf("%.1f", 60h))
        append!(row, [cell(R.nets[(state, ds)].leg_q, h) for ds in dets])
    end
    row
end
heads = ["Part of the leg", "Tissue (mL)"]
for (_, s) in states
    append!(heads, ["$s: healthy"; ["$s: " * sev_label(ds) for ds in dets]])
end
t_flow = html_table("Blood flow (mL/min, and % of healthy)", heads, [flow_rows; [leg_row]])

p_rows = map([0.0; dets]) do ds
    row = [ds == 0 ? "Healthy" : sev_label(ds)]
    for (state, _) in states
        n = R.nets[(state, ds)]
        append!(row, [@sprintf("%.1f", n.above), @sprintf("%.1f", n.below), @sprintf("%.1f", n.ankle), @sprintf("%.2f", n.ankle / R.p_root), @sprintf("%.1f", 60 * n.through)])
    end
    row
end
p_heads = ["Narrowing"]
for (_, s) in states
    append!(p_heads, ["$s: above it", "$s: below it", "$s: ankle", "$s: ankle ÷ inlet", "$s: flow through it (mL/min)"])
end
t_pressure = html_table(@sprintf("Mean pressure in the main artery (mmHg; %.0f mmHg enters the leg)", R.p_root), p_heads, p_rows)

o_rows = map(eachindex(R.part_names)) do p
    row = [R.part_names[p]]
    for ds in [0.0; dets]
        o = R.oxygen[ds][p]
        push!(row, @sprintf("%.2f / %.1f%%", o.so2, 100 * o.starved))
    end
    row
end
t_oxygen = html_table("At rest: venous O₂ saturation / beds short of oxygen (venous SO₂ below 0.05)", ["Part of the leg"; "Healthy"; sev_label.(dets)], o_rows)

items = [
    ("parts", "Flow in each part of the leg", "Each part's blood flow with the narrowing, as a share of healthy, at rest and at maximal dilation", joinpath(OUT, "flow_by_part.png")),
    ("curve", "Flow against the narrowing", "The flow to the tissue fed through the narrowing, and above it, from 0 to 95% narrowing", joinpath(OUT, "flow_vs_stenosis.png")),
    ("pressure", "Pressure along the main artery", "The main artery's mean pressure from the top of the leg to the foot, healthy and narrowed", joinpath(OUT, "pressure_along_leg.png")),
    ("maps", "Flow in every tissue patch", "The right leg, every tissue patch coloured by its blood flow as a share of healthy", joinpath(OUT, "flow_maps.png")),
]
page = write_stenosis_page(
    joinpath(OUT, "figures.html"), "A narrowed thigh artery: blood flow down to the foot",
    "Four figures for slides and the numbers behind them. 右键图片 → 复制图片，可以直接贴进 PPT（原分辨率）。", items, [t_flow, t_pressure, t_oxygen];
    foot = "Made from XCAT: keep it on lab storage.",
)
@printf("%s (%.1f MB)\n", page, filesize(page) / 1.0e6)
