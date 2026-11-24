# Slide figures of the leg's contrast transport at rest, from transport_rest's saved data (transport_view.jl's view_data.jls
# and transport_leg.jl's result.jls), drawn by transport_figs_plot.jl:
#   leg_snapshots.png          the leg at 15, 20, 30, 60 and 120 s
#   iodine_by_compartment.png  the iodine in the leg's arteries, tissue and veins over time
#   input_and_outflow.png      the leg's arterial input and its femoral venous outflow
#   arrival_map.png            when the iodine reaches each part of the leg
#
#     julia transport_figs.jl
#
# Numbers only. The figures are made from XCAT: they stay in LIMB_CACHE_DIR\transport_rest\figures (lab storage).

haskey(ENV, "LIMB_CACHE_DIR") || error("set LIMB_CACHE_DIR first")
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
import CairoMakie as Mke
using Serialization, Statistics, Printf
include(joinpath(@__DIR__, "transport_figs_plot.jl"))

const DIR = joinpath(ENV["LIMB_CACHE_DIR"], "transport_rest")
const OUT = mkpath(joinpath(DIR, "figures"))

V = deserialize(joinpath(DIR, "view_data.jls"))
R = deserialize(joinpath(DIR, "result.jls"))
pf = R.fp / (R.body.share * R.body.co) # plasma fraction: the leg's plasma flow over its blood flow
ks = [s.k for s in R.patches]
region = [s.region for s in R.patches]
Dp = R.D[ks] # each patch's arrival at its inlet, in V.tis order
order = [(4, "Thigh"), (3, "Knee"), (2, "Calf"), (1, "Foot")]
medians = [(name, median(Dp[region .== r])) for (r, name) in order]
cin = pf .* R.c
t_arrive = R.grid[findfirst(>=(0.1 * maximum(cin)), cin)]

println(fig_snapshots(V, [15, 20, 30, 60, 120], joinpath(OUT, "leg_snapshots.png")))
println(fig_compartments(V, joinpath(OUT, "iodine_by_compartment.png")))
println(fig_curves(collect(R.grid), cin, collect(R.tc), pf .* R.cs, joinpath(OUT, "input_and_outflow.png")))
println(fig_arrival(V, Dp, medians, t_arrive, joinpath(OUT, "arrival_map.png")))
@printf("plasma fraction %.3f; input arrives at %.1f s; median arrival %s\n", pf, t_arrive, join([@sprintf("%s %.1f s", n, m) for (n, m) in medians], ", "))

# the four on one page, at full resolution
items = [
    ("snapshots", "The leg over time", "The right leg at 15, 20, 30, 60 and 120 s after the injection starts, its main arteries, veins and tissue patches coloured by their iodine", joinpath(OUT, "leg_snapshots.png")),
    ("compartments", "Where the iodine in the leg is", "Iodine in the leg's arteries, tissue and veins, and the total, over time", joinpath(OUT, "iodine_by_compartment.png")),
    ("curves", "Iodine entering and leaving the leg", "The iodine concentration in the blood entering the leg and leaving it through the femoral vein, over time", joinpath(OUT, "input_and_outflow.png")),
    ("arrival", "When the iodine arrives", "The right leg coloured by when the iodine reaches its main arteries and its tissue patches", joinpath(OUT, "arrival_map.png")),
]
page = write_figures_page(
    joinpath(OUT, "figures.html"), "Contrast through the leg at rest",
    "Four figures for slides. 右键图片 → 复制图片，可以直接贴进 PPT（原分辨率）。", items;
    foot = "Interactive version, to play the injection over time: <a href=\"../transport_view.html\">transport_view.html</a>. Made from XCAT: keep it on lab storage.",
)
@printf("%s (%.1f MB)\n", page, filesize(page) / 1.0e6)
