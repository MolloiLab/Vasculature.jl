# Slide figures and tables of where the blood and the injected iodine go in the whole body at rest (whole_body_figs.jl
# draws them from whole_body.jl's result.jls). The functions take the result as a NamedTuple and need only CairoMakie (as
# `Mke`), so made-up data can check the layout. The helpers and the colours come from stenosis_figs_plot.jl and
# transport_figs_plot.jl.

include(joinpath(@__DIR__, "stenosis_figs_plot.jl"))

"The categorical palette's eight slots (light mode), in their fixed order: adjacent pairs validated."
const SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
"A label inside a fill: white, except primary ink on yellow."
inside_ink(hex) = hex == "#eda100" ? INK : "#ffffff"
"Milligrams: one decimal below 10, thousands separated above 1000."
fmt1(v) = v < 10 ? @sprintf("%.1f", v) : fmt_mg(v)
"The index of time `s` (s) in `t`."
time_index(t, s) = findfirst(x -> abs(x - s) < 1.0e-3, t)

"The cardiac output's split between the organs that take it: ranked bars, one hue."
function fig_blood_flow(R, path)
    rows = sort(R.flows; by = last) # the largest at the top
    n = length(rows)
    v = last.(rows)
    fig = Mke.Figure(backgroundcolor = SURFACE, fontsize = 15, figure_padding = 24)
    title!(fig[1, 1], "Where the blood goes at rest")
    ax = plain_axis(
        fig[2, 1]; xlabel = "Blood flow (mL/min)", width = 640, height = 38n + 10, yticks = (1:n, first.(rows)),
        ygridvisible = false, yticksvisible = false, leftspinevisible = false, yticklabelsize = 15,
    )
    Mke.barplot!(ax, 1:n, v; direction = :x, width = 0.62, color = SLOTS[1], strokewidth = 0)
    Mke.text!(ax, v, 1:n; text = [@sprintf("%s mL/min  (%.1f%%)", fmt_mg(x), 100x / R.co) for x in v], align = (:left, :center), offset = (6, 0), fontsize = 14, color = INK_2)
    Mke.xlims!(ax, 0, 1.42 * maximum(v))
    Mke.ylims!(ax, 0.4, n + 0.6)
    caption!(fig[3, 1], @sprintf(
        "At rest the heart pumps %s mL/min. All of it passes through the heart's chambers and the lungs; the bars split it between the organs that take it, as %% of that output. The right leg is the grown model at its own resting flow; the other organs keep the whole-body model's standard shares. Liver and gut: the liver's own artery plus the blood that reaches it from the gut.",
        fmt_mg(R.co),
    ))
    Mke.resize_to_layout!(fig)
    Mke.save(path, fig; px_per_unit = 2)
    return path
end

"The iodine in each organ over time, one panel each (own scale), the peak marked."
function fig_iodine_by_organ(R, path; cols = 4)
    n = length(R.organs)
    rows = cld(n, cols)
    fig = Mke.Figure(backgroundcolor = SURFACE, fontsize = 14, figure_padding = 24)
    title!(fig[1, 1], "The iodine in each organ, from the injection to 2 minutes")
    g = fig[2, 1] = Mke.GridLayout()
    for (k, o) in enumerate(R.organs)
        i, j = fldmod1(k, cols)
        s = R.stats[k]
        ax = plain_axis(
            g[i, j]; title = o, titlesize = 16, width = 270, height = 160, xticks = 0:30:120,
            xlabel = i == rows ? "Time after the injection starts (s)" : "", ylabel = j == 1 ? "Iodine (mg)" : "",
            xgridvisible = false, xticklabelsize = 12, yticklabelsize = 12, xlabelsize = 13, ylabelsize = 14,
        )
        Mke.lines!(ax, R.t, R.M[:, k]; color = SLOTS[1], linewidth = 2)
        Mke.scatter!(ax, [s.tpeak], [s.peak]; color = SLOTS[1], markersize = 9, strokecolor = SURFACE, strokewidth = 2)
        late = s.tpeak > 60 # label on the side away from the plot's edge
        Mke.text!(
            ax, s.tpeak, s.peak; text = @sprintf("%s mg (%.1f%%) at %.0f s", fmt1(s.peak), 100 * s.peak / R.dose, s.tpeak),
            align = late ? (:right, :bottom) : (:left, :bottom), offset = late ? (-2, 5) : (4, 5), fontsize = 12, color = INK_2,
        )
        Mke.xlims!(ax, 0, 120)
        Mke.ylims!(ax, 0, 1.3 * s.peak)
    end
    Mke.colgap!(g, 26)
    Mke.rowgap!(g, 16)
    caption!(fig[3, 1], @sprintf(
        "Iodine in each organ's blood and in the space around its cells, after the injection into the right arm starts; %% of the %s mg injected. Each panel has its own scale; the dot marks the most the organ holds. Heart muscle: the heart's own tissue and its veins. Large veins: the injected arm's veins, the veins to the heart and the venae cavae. The kidneys also clear %s mg (%.1f%%) into the urine by %.0f s.",
        fmt_mg(R.dose), fmt1(R.cleared[end]), 100 * R.cleared[end] / R.dose, R.t[end],
    ))
    Mke.resize_to_layout!(fig)
    Mke.save(path, fig; px_per_unit = 2)
    return path
end

"The share of all the iodine in each of the eight groups at the shown times: one stacked bar per time."
function fig_where(R, path)
    W = R.where
    nt = length(R.shown)
    ng = length(W.names)
    ng <= length(SLOTS) || error("more groups than the palette's slots")
    y0 = collect(nt:-1:1) # the first time at the top
    fig = Mke.Figure(backgroundcolor = SURFACE, fontsize = 15, figure_padding = 24)
    title!(fig[1, 1], "Where the iodine is, from the injection to 2 minutes")
    Mke.Legend(
        fig[2, 1], [Mke.PolyElement(color = c) for c in SLOTS[1:ng]], W.names;
        orientation = :horizontal, nbanks = 2, framevisible = false, labelcolor = INK_2, labelsize = 14, patchsize = (16, 12),
        colgap = 20, rowgap = 2, halign = :left, tellwidth = false, backgroundcolor = :transparent, padding = (0, 0, 0, 0),
    )
    ax = plain_axis(
        fig[3, 1]; xlabel = "Share of all the iodine at that time (%)", width = 860, height = 44nt + 10,
        yticks = (y0, ["$(s) s" for s in R.shown]), xticks = 0:20:100, ygridvisible = false, xgridvisible = false,
        yticksvisible = false, leftspinevisible = false, yticklabelsize = 15,
    )
    for i in 1:nt
        x = 0.0
        for k in 1:ng
            w = W.pct[k, i]
            w <= 0 && continue
            Mke.poly!(ax, Mke.Rect(x, y0[i] - 0.32, w, 0.64); color = SLOTS[k], strokecolor = SURFACE, strokewidth = 1.2)
            w >= 6 && Mke.text!(ax, x + w / 2, y0[i]; text = @sprintf("%.0f%%", w), align = (:center, :center), fontsize = 13, color = inside_ink(SLOTS[k]))
            x += w
        end
    end
    Mke.xlims!(ax, 0, 100)
    Mke.ylims!(ax, 0.4, nt + 0.6)
    caption!(fig[4, 1], "Each bar is all the iodine injected up to that time (the injection runs for the first 14 s), split by where it is: in each organ's blood and the space around its cells, and for the kidneys also what they have cleared into the urine. Legs: the grown right leg and the body model's left leg. Large vessels: the aorta and the large veins, the injected arm's among them. Shares of 6% and more are printed on the bars; the table below gives every one.")
    Mke.rowgap!(fig.layout, 1, 8)
    Mke.rowgap!(fig.layout, 2, 10)
    Mke.resize_to_layout!(fig)
    Mke.save(path, fig; px_per_unit = 2)
    return path
end

"The page's tables, as HTML strings."
function whole_body_tables(R)
    dose_pct(v) = @sprintf("%.1f%%", 100v / R.dose)
    times = ["$(s) s" for s in R.shown]
    flow = html_table(
        "Blood flow at rest", ["Organ", "mL/min", "% of the heart's output"],
        vcat(
            [["Heart chambers and lungs (all the blood passes through)", fmt_mg(R.co), "100%"]],
            [[o, fmt_mg(f), @sprintf("%.1f%%", 100f / R.co)] for (o, f) in sort(R.flows; by = last, rev = true)],
        ),
    )
    organ_rows = [
        [o, @sprintf("%.1f", s.arrive), fmt1(s.peak), dose_pct(s.peak), @sprintf("%.1f", s.tpeak), [fmt1(v) for v in s.at]...]
            for (o, s) in zip(R.organs, R.stats)
    ]
    push!(organ_rows, ["Cleared into urine", "", "", "", "", [fmt1(v) for v in R.cleared_stats.at]...])
    organs = html_table(
        @sprintf("Iodine in each organ, mg (the peak also as %% of the %s mg injected); arrival: when it first holds 10%% of its peak", fmt_mg(R.dose)),
        ["Organ", "Arrives (s)", "Peak (mg)", "Peak (%)", "Peak at (s)", times...], organ_rows,
    )
    where = html_table(
        "Where the iodine is: % of all the iodine at that time", ["Group", times...],
        [[w, [@sprintf("%.1f%%", R.where.pct[k, i]) for i in eachindex(R.shown)]...] for (k, w) in enumerate(R.where.names)],
    )
    iS = [time_index(R.t, s) for s in R.shown]
    leg = html_table(
        @sprintf("The grown right leg by part: mg of iodine (the leg takes %.0f mL/min at rest)", R.leg.ml_min), ["Part", times...],
        vcat(
            [[p, [fmt1(R.leg.mg[i, k]) for i in iS]...] for (k, p) in enumerate(R.leg.parts)],
            [["Total", [fmt1(R.leg.total[i]) for i in iS]...]],
        ),
    )
    return [flow, organs, where, leg]
end
