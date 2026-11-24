# Slide figures of the leg with a narrowed superficial femoral artery (stenosis_figs.jl draws them from stenosis_leg.jl's
# result.jls). The functions take the result as a NamedTuple and need only CairoMakie (as `Mke`), so made-up data can
# check the layout. The helpers and the colours come from transport_figs_plot.jl.

include(joinpath(@__DIR__, "transport_figs_plot.jl"))

"The stenoses shown in detail, on the sequential blue ramp's ordinal steps 250, 450 and 650."
const SEV_COLOURS = [BLUES[4], BLUES[8], BLUES[12]]
const HEALTHY = TOTAL           # secondary ink: the reference, not a series of its own
const THROUGH = ARTERY          # categorical slot 2: the tissue fed through the narrowing
const ABOVE = VEIN              # slot 1: the tissue fed above it
const DIV_LOW = ["#7a1f1f", "#e34948", "#cfcdc5"] # diverging, red arm: less flow, to the neutral midpoint
const DIV_HIGH = ["#cfcdc5", "#2a78d6"]           # blue arm: more flow (the midpoint a shade darker than the leg, so unchanged dots show)

pct(x) = @sprintf("%.0f%%", 100x)
sev_label(ds) = @sprintf("%.0f%% narrower", 100ds)
"The detail stenoses' indices among the severities."
detail_nets(R, state) = [R.nets[(state, ds)] for ds in R.detail_severities]

"An axis with the house style: recessive grid and spines, secondary-ink labels."
function plain_axis(pos; kw...)
    return Mke.Axis(
        pos; backgroundcolor = SURFACE, titlesize = 19, titlefont = :bold, titlecolor = INK, titlealign = :left,
        xlabelcolor = INK_2, ylabelcolor = INK_2, xticklabelcolor = INK_2, yticklabelcolor = INK_2, xgridcolor = GRIDLINE,
        ygridcolor = GRIDLINE, bottomspinecolor = BASELINE, leftspinecolor = BASELINE, topspinevisible = false,
        rightspinevisible = false, xtickcolor = BASELINE, ytickcolor = BASELINE, xlabelsize = 16, ylabelsize = 16,
        xticklabelsize = 14, yticklabelsize = 14, kw...,
    )
end

caption!(pos, text) = Mke.Label(pos, text; fontsize = 14, color = INK_2, halign = :left, justification = :left, word_wrap = true, tellwidth = false)
title!(pos, text) = Mke.Label(pos, text; fontsize = 21, font = :bold, color = INK, halign = :left, tellwidth = false)

"Each part's flow at the detail stenoses as a share of its healthy flow, at rest and at maximal dilation (grouped bars)."
function fig_flow_by_part(R, path)
    fig = Mke.Figure(backgroundcolor = SURFACE, fontsize = 15, figure_padding = 24)
    title!(fig[1, 1:2], "Blood flow in each part of the leg, with the thigh artery narrowed")
    nparts = length(R.part_names)
    y0 = collect(nparts:-1:1) # the thigh at the top, the foot at the bottom
    bw = 0.26
    offs = [bw, 0.0, -bw]
    for (c, (state, heading)) in enumerate(((:rest, "At rest"), (:stress, "At maximal dilation")))
        h = R.nets[(state, 0.0)].part_q
        ax = plain_axis(
            fig[3, c]; title = heading, xlabel = "Blood flow (% of the healthy leg's)", width = 470, height = 430,
            yticks = (y0, [@sprintf("%s\nhealthy %.0f mL/min", n, 60q) for (n, q) in zip(R.part_short, h)]),
            xticks = 0:25:100, ygridvisible = false, yticksvisible = false, leftspinevisible = false,
        )
        for (i, n) in enumerate(detail_nets(R, state))
            v = 100 .* n.part_q ./ h
            Mke.barplot!(ax, y0 .+ offs[i], v; direction = :x, width = bw * 0.92, color = SEV_COLOURS[i], strokewidth = 0)
            Mke.text!(ax, v, y0 .+ offs[i]; text = [@sprintf("%.0f%%", x) for x in v], align = (:left, :center), offset = (5, 0), fontsize = 13, color = INK_2)
        end
        Mke.vlines!(ax, [100]; color = BASELINE, linewidth = 1.5, linestyle = :dash)
        Mke.xlims!(ax, 0, 128)
        Mke.ylims!(ax, 0.45, nparts + 0.55)
    end
    Mke.Legend(
        fig[2, 1:2], [Mke.PolyElement(color = c) for c in SEV_COLOURS], sev_label.(R.detail_severities);
        orientation = :horizontal, framevisible = false, labelcolor = INK_2, labelsize = 15, patchsize = (18, 12), colgap = 24,
        halign = :left, tellwidth = false, backgroundcolor = :transparent, padding = (0, 0, 0, 0),
    )
    caption!(fig[4, 1:2], @sprintf(
        "The right leg's superficial femoral artery narrowed over %.0f mm, %.0f cm above the knee joint (the distal adductor canal); the narrowing's diameter as a share of the artery's. The thigh fed above it takes its blood from the deep femoral artery and the branches above. Rest: the arterioles regulate the flow to the tissue's oxygen need. Maximal dilation: every arteriole fully open, as under adenosine.",
        sum(s.mm for s in R.site), R.site_above_knee / 10,
    ))
    Mke.colgap!(fig.layout, 36)
    Mke.rowgap!(fig.layout, 1, 8)
    Mke.resize_to_layout!(fig)
    Mke.save(path, fig; px_per_unit = 2)
    return path
end

"The flow to the tissue fed through the narrowing, and above it, against the stenosis, at rest and at maximal dilation."
function fig_flow_vs_stenosis(R, path)
    fig = Mke.Figure(backgroundcolor = SURFACE, fontsize = 15, figure_padding = 24)
    title!(fig[1, 1:2], "How much the narrowing takes, from mild to severe")
    x = 100 .* R.severities
    lo = 100.0
    hi = 100.0
    axes_ = Any[]
    for (c, (state, heading)) in enumerate(((:rest, "At rest"), (:stress, "At maximal dilation")))
        h = R.nets[(state, 0.0)]
        share(n, sel) = 100 * sum(n.Q[sel]) / sum(h.Q[sel])
        thr = [share(R.nets[(state, ds)], R.fed) for ds in R.severities]
        abv = [share(R.nets[(state, ds)], .!R.fed) for ds in R.severities]
        lo, hi = min(lo, minimum(thr)), max(hi, maximum(abv))
        ax = plain_axis(fig[3, c]; title = heading, xlabel = "Narrowing (% of the artery's diameter)", ylabel = c == 1 ? "Blood flow (% of healthy)" : "", xticks = 0:10:90, yticks = 0:20:120, width = 470, height = 360)
        push!(axes_, ax)
        Mke.hlines!(ax, [100]; color = BASELINE, linewidth = 1.5, linestyle = :dash)
        Mke.lines!(ax, x, abv; color = ABOVE, linewidth = 2.5)
        Mke.lines!(ax, x, thr; color = THROUGH, linewidth = 2.5)
        for (i, ds) in enumerate(R.detail_severities)
            j = findfirst(==(ds), R.severities)
            Mke.scatter!(ax, [x[j]], [thr[j]]; color = THROUGH, markersize = 9, strokecolor = SURFACE, strokewidth = 2)
            low = thr[j] < 15 # near the axis: label above and to the right
            Mke.text!(ax, x[j], thr[j]; text = @sprintf("%.0f%%", thr[j]), align = low ? (:left, :bottom) : (:right, :top), offset = low ? (6, 4) : (-6, -4), fontsize = 13, color = INK_2)
        end
    end
    for ax in axes_
        Mke.xlims!(ax, 0, 95)
        Mke.ylims!(ax, 0, max(115, 5 * ceil(hi / 5) + 5))
    end
    Mke.Legend(
        fig[2, 1:2], [Mke.LineElement(color = THROUGH, linewidth = 3), Mke.LineElement(color = ABOVE, linewidth = 3)],
        ["Tissue fed through the narrowing (lower thigh, knee, calf, foot)", "Tissue fed above it (upper thigh)"];
        orientation = :horizontal, framevisible = false, labelcolor = INK_2, labelsize = 15, patchsize = (24, 3), colgap = 28,
        halign = :left, tellwidth = false, backgroundcolor = :transparent, padding = (0, 0, 0, 0),
    )
    caption!(fig[4, 1:2], "Flow as a share of the same tissue's healthy flow. At rest the arterioles below the narrowing open up as the pressure falls, so the flow holds until they run out of dilation. At maximal dilation they are already open, so every bit of pressure the narrowing takes costs flow. With the flow through the narrowing cut, less pressure is lost above it, so the tissue above it gains a little.")
    Mke.colgap!(fig.layout, 36)
    Mke.rowgap!(fig.layout, 1, 8)
    Mke.resize_to_layout!(fig)
    Mke.save(path, fig; px_per_unit = 2)
    return path
end

"The trunk's pressure as a polyline over distance from the top of the leg (cm): each segment's inlet and end pressure."
function trunk_line(R, n)
    xs, ys = Float64[0.0], Float64[n.trunk_pin[1]]
    for i in eachindex(R.trunk)
        x0 = i == 1 ? 0.0 : R.trunk_mm[i - 1]
        push!(xs, x0 / 10, R.trunk_mm[i] / 10)
        push!(ys, n.trunk_pin[i], n.trunk_pend[i])
    end
    return xs, ys
end

"The pressure in the main artery from the top of the leg to the foot, healthy and at the detail stenoses, at rest and at maximal dilation."
function fig_pressure(R, path)
    fig = Mke.Figure(backgroundcolor = SURFACE, fontsize = 15, figure_padding = 24)
    title!(fig[1, 1:2], "Pressure along the main artery, from the top of the leg to the foot")
    xmax = R.trunk_mm[end] / 10
    for (c, (state, heading)) in enumerate(((:rest, "At rest"), (:stress, "At maximal dilation")))
        ax = plain_axis(fig[3, c]; title = heading, xlabel = "Distance along the artery from the top of the leg (cm)", ylabel = c == 1 ? "Mean pressure (mmHg)" : "", yticks = 0:20:100, width = 470, height = 360)
        Mke.vspan!(ax, [R.site_from / 10], [R.site_to / 10]; color = (INK, 0.10))
        for (d, name) in ((R.knee_mm, "knee joint"), (R.ankle_mm, "ankle"))
            Mke.vlines!(ax, [d / 10]; color = BASELINE, linewidth = 1.2, linestyle = :dot)
            Mke.text!(ax, d / 10, 99; text = name, align = (:left, :top), offset = (4, 0), fontsize = 13, color = INK_2)
        end
        Mke.hlines!(ax, [R.p_vein]; color = BASELINE, linewidth = 1.2, linestyle = :dash)
        Mke.text!(ax, 0, R.p_vein; text = "vein", align = (:left, :top), offset = (4, -3), fontsize = 13, color = INK_2)
        xs, ys = trunk_line(R, R.nets[(state, 0.0)])
        Mke.lines!(ax, xs, ys; color = HEALTHY, linewidth = 2.5)
        for (i, n) in enumerate(detail_nets(R, state))
            xs, ys = trunk_line(R, n)
            Mke.lines!(ax, xs, ys; color = SEV_COLOURS[i], linewidth = 2.5)
        end
        Mke.text!(ax, R.site_from / 10, 99; text = "narrowing", align = (:right, :top), offset = (-4, 0), fontsize = 13, color = INK_2)
        Mke.xlims!(ax, 0, xmax)
        Mke.ylims!(ax, 0, 100)
    end
    Mke.Legend(
        fig[2, 1:2], [Mke.LineElement(color = c, linewidth = 3) for c in [HEALTHY; SEV_COLOURS]], ["Healthy"; sev_label.(R.detail_severities)];
        orientation = :horizontal, framevisible = false, labelcolor = INK_2, labelsize = 15, patchsize = (24, 3), colgap = 28,
        halign = :left, tellwidth = false, backgroundcolor = :transparent, padding = (0, 0, 0, 0),
    )
    caption!(fig[4, 1:2], @sprintf(
        "The femoral, popliteal and tibial arteries, followed to the foot. Grey band: the narrowing. %.0f mmHg enters the leg; blood flows only while the arterial pressure is above the vein's (%.0f mmHg). Rest needs little flow, so little pressure is lost; maximal flow loses more even when healthy.",
        R.p_root, R.p_vein,
    ))
    Mke.colgap!(fig.layout, 36)
    Mke.rowgap!(fig.layout, 1, 8)
    Mke.resize_to_layout!(fig)
    Mke.save(path, fig; px_per_unit = 2)
    return path
end

"A diverging colour map for flow as % of healthy over (0, top): red below 100%, the neutral grey at 100%, blue above."
function flow_colours(top; n = 256)
    lo = [Mke.RGBf(hexrgb(h)...) for h in DIV_LOW]
    hi = [Mke.RGBf(hexrgb(h)...) for h in DIV_HIGH]
    mix(a, b, t) = Mke.RGBf(a.r + (b.r - a.r) * t, a.g + (b.g - a.g) * t, a.b + (b.b - a.b) * t)
    return map(range(0, top; length = n)) do v
        if v <= 100
            u = v / 100 * (length(lo) - 1)
            i = min(floor(Int, u), length(lo) - 2)
            mix(lo[i + 1], lo[i + 2], u - i)
        else
            mix(hi[1], hi[2], (v - 100) / (top - 100))
        end
    end
end

"Map bounds from the silhouette and the arteries."
function map_bounds(R; pad = 8.0)
    xs = vcat(R.art.x1, R.art.x2, [r[i] for r in R.sil.rows for i in 2:length(r)])
    zs = vcat(R.art.z1, R.art.z2, [r[1] for r in R.sil.rows])
    return (x0 = minimum(xs) - pad, x1 = maximum(xs) + R.sil.h + pad, z0 = minimum(zs) - pad, z1 = maximum(zs) + R.sil.h + pad)
end

"The leg at the detail stenoses, each patch by its flow as % of healthy, at rest and at maximal dilation (one row)."
function fig_maps(R, path; tall = 640)
    B = map_bounds(R)
    wide = round(Int, tall * (B.x1 - B.x0) / (B.z1 - B.z0))
    pxmm = tall / (B.z1 - B.z0)
    V = (; sil = R.sil)
    ks = R.ks
    vals = Dict((state, ds) => 100 .* R.nets[(state, ds)].Q[ks] ./ R.nets[(state, 0.0)].Q[ks] for state in (:rest, :stress) for ds in R.detail_severities)
    top = max(110.0, 10 * ceil(maximum(maximum, values(vals)) / 10))
    cmap = flow_colours(top)
    fig = Mke.Figure(backgroundcolor = SURFACE, fontsize = 15, figure_padding = 24)
    title!(fig[1, 1:7], "Blood flow in every tissue patch, as a share of its healthy flow")
    oa = sortperm(R.art.r)
    nd = length(R.detail_severities)
    for (g, (state, heading)) in enumerate(((:rest, "At rest"), (:stress, "At maximal dilation")))
        Mke.Label(fig[2, (g - 1) * nd + 1:g * nd + (g == 1 ? 0 : 0)], heading; fontsize = 18, font = :bold, color = INK, halign = :center, tellwidth = false)
        for (i, ds) in enumerate(R.detail_severities)
            ax = leg_axis(fig[3, (g - 1) * nd + i], V, B, wide, tall; title = sev_label(ds))
            ax.titlesize = 16
            ax.titlefont = :regular
            ax.titlecolor = INK_2
            Mke.linesegments!(ax, segpoints(R.art, oa); color = Mke.RGBf(hexrgb(GRIDLINE)...), linewidth = repeat(width.(R.art.r[oa], pxmm; thinnest = 0.3); inner = 2))
            v = vals[(state, ds)]
            o = sortperm(abs.(v .- 100)) # the changed patches on top
            Mke.scatter!(ax, R.pos[1, ks[o]], R.pos[3, ks[o]]; color = clamp.(v[o], 0, top), colormap = cmap, colorrange = (0, top), markersize = clamp.(1.8 .* cbrt.(R.ml[ks[o]]), 2.5, 7), strokewidth = 0)
            s1, s2 = first(R.site), last(R.site)
            Mke.scatter!(ax, [(s1.xa + s2.xb) / 2], [(s1.za + s2.zb) / 2]; marker = :hline, markersize = 16, color = INK, strokewidth = 0)
        end
    end
    Mke.Colorbar(
        fig[3, 2nd + 1]; colormap = cmap, limits = (0, top), label = "Blood flow (% of healthy)", ticks = 0:20:top,
        labelcolor = INK_2, ticklabelcolor = INK_2, tickcolor = BASELINE, spinewidth = 0, labelsize = 15, ticklabelsize = 14,
        width = 14, height = Mke.Relative(0.6),
    )
    caption!(fig[4, 1:7], "Right leg, front view. Dots: the tissue patches, at the end of each one's main artery. Black tick: the narrowing. Grey lines: the main arteries. Red: less flow than healthy; grey: unchanged; blue: more.")
    Mke.colgap!(fig.layout, 10)
    Mke.colgap!(fig.layout, nd, 34)
    Mke.colgap!(fig.layout, 2nd, 24)
    Mke.resize_to_layout!(fig)
    Mke.save(path, fig; px_per_unit = 2)
    return path
end

# ---------------------------------------------------------------------------------------------------------------------
# the page: the figures and the tables

"An HTML table: `head` cells, then `rows` of cells (strings); the first column left-aligned, the rest right."
function html_table(caption, head, rows)
    io = IOBuffer()
    print(io, "<table><caption>", _esc(caption), "</caption><thead><tr>")
    foreach(h -> print(io, "<th>", _esc(h), "</th>"), head)
    print(io, "</tr></thead><tbody>")
    for r in rows
        print(io, "<tr>")
        foreach(c -> print(io, "<td>", _esc(c), "</td>"), r)
        print(io, "</tr>")
    end
    print(io, "</tbody></table>")
    return String(take!(io))
end

"`write_figures_page` with tables (HTML strings) after the figures."
function write_stenosis_page(path, title, sub, items, tables; foot = "")
    write_figures_page(path, title, sub, items; foot)
    page = read(path, String)
    style = """
table { border-collapse: collapse; margin: 0 0 24px; background: var(--surface); border: 1px solid var(--border); border-radius: 12px; font-variant-numeric: tabular-nums; display: block; overflow-x: auto; max-width: 100%; }
caption { text-align: left; font-weight: 600; padding: 12px 14px 4px; font-size: 15px; }
th, td { padding: 6px 14px; text-align: right; white-space: nowrap; border-bottom: 1px solid var(--border); font-size: 14px; }
th { color: var(--ink2); font-weight: 600; }
th:first-child, td:first-child { text-align: left; }
tbody tr:last-child td { border-bottom: none; }
"""
    page = replace(page, "</style>" => style * "</style>"; count = 1)
    marker = "<p class=\"foot\">"
    block = "<h2 class=\"tables\">Numbers</h2>\n" * join(tables, "\n") * "\n"
    page = occursin(marker, page) ? replace(page, marker => block * marker; count = 1) : replace(page, "</main>" => block * "</main>"; count = 1)
    write(path, page)
    return path
end
