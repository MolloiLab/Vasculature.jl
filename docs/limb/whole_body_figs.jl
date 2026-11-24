# Slide figures and tables of where the blood and the injected iodine go in the whole body at rest, from whole_body.jl's
# result.jls, drawn by whole_body_figs_plot.jl:
#   blood_flow_by_organ.png   the cardiac output's split between the organs that take it
#   iodine_by_organ.png       the iodine in each organ over time, one panel each
#   where_the_iodine_is.png   the share of all the iodine in eight groups at 5, 10, 15, 20, 30, 60 and 120 s
#   figures.html              the three at full resolution, and the tables
#
#     julia whole_body_figs.jl
#
# Numbers only. Made from XCAT: everything stays in LIMB_CACHE_DIR\whole_body\figures (lab storage).

haskey(ENV, "LIMB_CACHE_DIR") || error("set LIMB_CACHE_DIR first")
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
import CairoMakie as Mke
using Serialization, Statistics, Printf
include(joinpath(@__DIR__, "whole_body_figs_plot.jl"))

const DIR = joinpath(ENV["LIMB_CACHE_DIR"], "whole_body")
const OUT = mkpath(joinpath(DIR, "figures"))

R = deserialize(joinpath(DIR, "result.jls"))
items = [
    ("flow", "Where the blood goes", "The cardiac output at rest split between the organs that take it, in mL/min and as a share of the output", fig_blood_flow(R, joinpath(OUT, "blood_flow_by_organ.png"))),
    ("organs", "The iodine in each organ", "Iodine in each organ over the first two minutes after the injection, one panel per organ, the most each holds marked", fig_iodine_by_organ(R, joinpath(OUT, "iodine_by_organ.png"))),
    ("where", "Where the iodine is", "The share of all the iodine in the heart, lungs, legs, liver and gut, kidneys and urine, brain, arms and the rest of the body, and the large vessels, at 5 to 120 s", fig_where(R, joinpath(OUT, "where_the_iodine_is.png"))),
]
foreach(it -> println(last(it)), items)
page = write_stenosis_page(
    joinpath(OUT, "figures.html"), "Where the blood and the iodine go: the whole body at rest",
    "Three figures for slides, and the numbers behind them. 右键图片 → 复制图片，可以直接贴进 PPT（原分辨率）。", items, whole_body_tables(R);
    foot = @sprintf(
        "The whole-body model with the grown right leg in its place; %s. Checks: the iodine in the organs plus what the kidneys cleared matches the iodine injected to within %.0f mg (%.1f%%). Made from XCAT: keep it on lab storage.",
        R.injection, abs(R.checks.balance), 100 * abs(R.checks.balance) / R.dose,
    ),
)
@printf("%s (%.1f MB)\n", page, filesize(page) / 1.0e6)
