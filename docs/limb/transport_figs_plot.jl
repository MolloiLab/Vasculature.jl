# Slide figures of the leg's contrast transport at rest (transport_figs.jl draws them from the saved data). The functions
# take the data as NamedTuples and need only CairoMakie (as `Mke`), so made-up data can check the layout.
#   V: transport_view.jl's view data (vt, aif, aif_dt, art, vein, tis, sil, series, injected, leg_ml_min, npatch)
#   R: what the figures need from transport_leg.jl's result (input c on grid, outflow cs on tc, D per patch, regions)

using Printf

const SURFACE = "#fcfcfb"
const INK = "#0b0b0b"
const INK_2 = "#52514e"
const MUTED = "#898781"
const BASELINE = "#c3c2b7"
const GRIDLINE = "#e1e0d9"
const LEG_FILL = "#efeee9"
const EMPTY = "#d3d2ca"
const ARTERY = "#eb6834" # categorical slot 2
const VEIN = "#2a78d6"   # slot 1
const TISSUE = "#1baf7a" # slot 3
const TOTAL = "#52514e"  # the sum: secondary ink, a reference rather than a series of its own
const BLUES = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]

hexrgb(h) = Tuple(parse(Int, h[i:(i + 1)]; base = 16) / 255 for i in (2, 4, 6))
"A one-hue ramp from the empty grey to `hue`, so magnitude reads as more colour."
ramp(hue; n = 101) = (a = hexrgb(EMPTY); b = hexrgb(hue); [Mke.RGBf((a .+ (b .- a) .* (i / (n - 1)))...) for i in 0:(n - 1)])

"Segment end points, thin first so the wide vessels lie on top."
function segpoints(s, order)
    pts = Vector{Mke.Point2f}(undef, 2 * length(order))
    for (i, j) in enumerate(order)
        pts[2i - 1] = Mke.Point2f(s.x1[j], s.z1[j])
        pts[2i] = Mke.Point2f(s.x2[j], s.z2[j])
    end
    return pts
end

"The leg's front-view silhouette as a 0/1 image with its cell centres."
function silhouette_image(sil, bounds)
    h = sil.h
    xs = collect((bounds.x0 + h / 2):h:bounds.x1)
    zs = collect((bounds.z0 + h / 2):h:bounds.z1)
    img = zeros(Float32, length(xs), length(zs))
    for row in sil.rows
        k = clamp(round(Int, (row[1] + h / 2 - zs[1]) / h) + 1, 1, length(zs))
        for i in 2:2:length(row)
            a = clamp(round(Int, (row[i] + h / 2 - xs[1]) / h) + 1, 1, length(xs))
            b = clamp(round(Int, (row[i + 1] + h / 2 - xs[1]) / h) + 1, 1, length(xs))
            img[a:b, k] .= 1
        end
    end
    return xs, zs, img
end

function leg_bounds(V; pad = 8.0)
    xs = vcat(V.art.x1, V.art.x2, V.vein.x1, V.vein.x2, [r[i] for r in V.sil.rows for i in 2:length(r)])
    zs = vcat(V.art.z1, V.art.z2, V.vein.z1, V.vein.z2, [r[1] for r in V.sil.rows])
    return (x0 = minimum(xs) - pad, x1 = maximum(xs) + V.sil.h + pad, z0 = minimum(zs) - pad, z1 = maximum(zs) + V.sil.h + pad)
end

aif_at(V, s) = s <= 0 ? 0.0 : (x = s / V.aif_dt; i = floor(Int, x); i >= length(V.aif) - 1 ? V.aif[end] : V.aif[i + 1] + (x - i) * (V.aif[i + 2] - V.aif[i + 1]))
frame_of(V, t) = argmin(abs.(V.vt .- t))

"A leg axis at fixed size, no decorations, with the silhouette."
function leg_axis(pos, V, B, wide, tall; title = "")
    ax = Mke.Axis(pos; aspect = Mke.DataAspect(), backgroundcolor = SURFACE, title, titlesize = 22, titlecolor = INK, titlefont = :bold, width = wide, height = tall)
    Mke.hidedecorations!(ax)
    Mke.hidespines!(ax)
    xs, zs, img = silhouette_image(V.sil, B)
    Mke.heatmap!(ax, xs, zs, img; colormap = [Mke.RGBf(hexrgb(SURFACE)...), Mke.RGBf(hexrgb(LEG_FILL)...)], colorrange = (0, 1))
    Mke.limits!(ax, B.x0, B.x1, B.z0, B.z1)
    return ax
end

"Line width (figure units) for radius r (mm): the vessel's own width on the drawing, at least `thinnest`."
width(r, pxmm; thinnest = 0.45) = max(thinnest, 2r * pxmm)

"The leg at time t: tissue dots, veins, then arteries, each on its own one-hue square-root ramp."
function leg_at!(ax, V, t, pxmm, maxima)
    f = frame_of(V, t)
    tv = sqrt.(clamp.(V.tis.conc[f, :] ./ maxima.tis, 0, 1))
    Mke.scatter!(ax, V.tis.x, V.tis.z; color = tv, colormap = ramp(TISSUE), colorrange = (0, 1), markersize = clamp.(2.2 .* cbrt.(V.tis.ml), 3, 8), strokewidth = 0)
    ov = sortperm(V.vein.r)
    vv = sqrt.(clamp.(V.vein.conc[f, ov] ./ maxima.vein, 0, 1))
    Mke.linesegments!(ax, segpoints(V.vein, ov); color = repeat(vv; inner = 2), colormap = ramp(VEIN), colorrange = (0, 1), linewidth = repeat(width.(V.vein.r[ov], pxmm); inner = 2))
    oa = sortperm(V.art.r)
    av = [sqrt(clamp(aif_at(V, t - V.art.tau[j]) / maxima.art, 0, 1)) for j in oa]
    Mke.linesegments!(ax, segpoints(V.art, oa); color = repeat(av; inner = 2), colormap = ramp(ARTERY), colorrange = (0, 1), linewidth = repeat(width.(V.art.r[oa], pxmm); inner = 2))
    return nothing
end

"A horizontal colour bar for a square-root ramp: positions sqrt(value / max), labelled in values."
function sqrt_bar!(pos, hue, maxv, label)
    nice = maxv >= 10 ? [0, 1, 5, 10, 15] : maxv >= 2 ? [0, 0.5, 1, 2, 3] : maxv >= 1 ? [0, 0.25, 0.5, 1, 1.5] : [0, 0.01, 0.025, 0.05, 0.1]
    vals = filter(<=(maxv * 1.0001), nice)
    Mke.Colorbar(
        pos; colormap = ramp(hue), limits = (0, 1), vertical = false, flipaxis = false, label,
        ticks = (sqrt.(vals ./ maxv), [v == round(v) ? string(Int(v)) : string(v) for v in vals]),
        labelcolor = INK_2, ticklabelcolor = INK_2, tickcolor = BASELINE, spinewidth = 0, labelsize = 17, ticklabelsize = 15, height = 12, width = 300,
    )
end

"Snapshots of the leg at `times` (s), side by side, with their three colour scales."
function fig_snapshots(V, times, path; tall = 860)
    B = leg_bounds(V)
    wide = round(Int, tall * (B.x1 - B.x0) / (B.z1 - B.z0))
    pxmm = tall / (B.z1 - B.z0)
    maxima = (art = maximum(V.aif), vein = maximum(V.vein.conc), tis = maximum(V.tis.conc))
    fig = Mke.Figure(backgroundcolor = SURFACE, fontsize = 16, figure_padding = 24)
    for (i, t) in enumerate(times)
        ax = leg_axis(fig[1, i], V, B, wide, tall; title = @sprintf("%g s", t))
        leg_at!(ax, V, t, pxmm, maxima)
    end
    g = fig[2, 1:length(times)] = Mke.GridLayout()
    sqrt_bar!(g[1, 1], ARTERY, maxima.art, "Arteries: iodine in blood (mg/mL)")
    sqrt_bar!(g[1, 2], VEIN, maxima.vein, "Veins: iodine in blood (mg/mL)")
    sqrt_bar!(g[1, 3], TISSUE, maxima.tis, "Tissue: iodine per mL of tissue (mg/mL)")
    Mke.colgap!(g, 48)
    Mke.Label(
        fig[3, 1:length(times)],
        @sprintf("Right leg at rest, front view, times after the injection starts. Lines: the main arteries and veins, drawn at their radius. Dots: the %d tissue patches. Grey: no iodine yet; square-root colour scales.", V.npatch);
        fontsize = 15, color = INK_2, halign = :left, justification = :left, word_wrap = true, tellwidth = false,
    )
    Mke.colgap!(fig.layout, 18)
    Mke.rowgap!(fig.layout, 22)
    Mke.resize_to_layout!(fig)
    Mke.save(path, fig; px_per_unit = 2)
    return path
end

"A time axis (s, 0-120) for a chart."
function time_axis(pos, ylabel; title = "")
    return Mke.Axis(
        pos; title, titlesize = 20, titlefont = :bold, titlecolor = INK, titlealign = :left, backgroundcolor = SURFACE,
        xlabel = "Time after the injection starts (s)", ylabel, xticks = 0:20:120, xlabelcolor = INK_2, ylabelcolor = INK_2,
        xticklabelcolor = INK_2, yticklabelcolor = INK_2, xgridcolor = GRIDLINE, ygridcolor = GRIDLINE, xgridvisible = false,
        bottomspinecolor = BASELINE, leftspinecolor = BASELINE, topspinevisible = false, rightspinevisible = false,
        xtickcolor = BASELINE, ytickcolor = BASELINE, xlabelsize = 17, ylabelsize = 17, xticklabelsize = 15, yticklabelsize = 15,
        width = 760, height = 400,
    )
end

fmt_mg(v) = v >= 1000 ? replace(@sprintf("%.0f", v), r"(\d)(?=(\d{3})+$)" => s"\1,") : @sprintf("%.0f", v)

"Iodine in the leg's arteries, tissue and veins, and the total, over time."
function fig_compartments(V, path)
    fig = Mke.Figure(backgroundcolor = SURFACE, fontsize = 16, figure_padding = 24)
    Mke.Label(fig[1, 1], "Where the iodine in the leg is, at rest"; fontsize = 20, font = :bold, color = INK, halign = :left, tellwidth = false)
    ax = time_axis(fig[3, 1], "Iodine in the leg (mg)")
    t = V.vt
    S = V.series
    lines = [("Total", S.total, TOTAL), ("Tissue", S.tis, TISSUE), ("Veins", S.ven, VEIN), ("Arteries", S.art, ARTERY)]
    for (_, y, col) in lines
        Mke.lines!(ax, t, y; color = col, linewidth = 3)
    end
    Mke.Legend(
        fig[2, 1], [Mke.LineElement(color = col, linewidth = 3) for (_, _, col) in lines],
        [@sprintf("%s: %s mg at %g s", name, fmt_mg(y[end]), t[end]) for (name, y, _) in lines];
        orientation = :horizontal, nbanks = 2, framevisible = false, labelcolor = INK_2, labelsize = 15, patchsize = (24, 3), colgap = 28, rowgap = 2,
        halign = :left, tellwidth = false, backgroundcolor = :transparent, padding = (0, 0, 0, 0),
    )
    Mke.xlims!(ax, 0, 120)
    Mke.ylims!(ax, 0, nothing)
    i = argmax(S.art)
    Mke.Label(
        fig[4, 1],
        @sprintf("The arteries hold the most at %.0f s (%s mg). Arteries: the main and the patches' arteries. Tissue: arterioles, capillaries, the space between the cells, venules. Veins: the patches' and the main veins. %s mg of iodine injected; the leg takes %.0f mL/min of blood at rest.", t[i], fmt_mg(S.art[i]), fmt_mg(V.injected), V.leg_ml_min);
        fontsize = 14, color = INK_2, halign = :left, justification = :left, word_wrap = true, tellwidth = false,
    )
    Mke.rowgap!(fig.layout, 1, 6)
    Mke.resize_to_layout!(fig)
    Mke.save(path, fig; px_per_unit = 2)
    return path
end

"The leg's arterial input and its femoral venous outflow (blood, mg/mL) over time."
function fig_curves(tin, cin, tout, cout, path)
    fig = Mke.Figure(backgroundcolor = SURFACE, fontsize = 16, figure_padding = 24)
    ax = time_axis(fig[1, 1], "Iodine in blood (mg/mL)"; title = "Iodine entering and leaving the leg, at rest")
    Mke.lines!(ax, tin, cin; color = ARTERY, linewidth = 3, label = "Input: the artery at the top of the leg")
    Mke.lines!(ax, tout, cout; color = VEIN, linewidth = 3, label = "Output: the femoral vein")
    # the input's peak labelled to its right; the output's on whichever side of its curve the input does not run
    io = argmax(cout)
    j = clamp(searchsortedfirst(tin, tout[io]), 1, length(tin))
    side = cin[j] > cout[io] ? ((:center, :top), (0, -10)) : ((:center, :bottom), (0, 10))
    for (t, c, col, what, align, offset) in ((tin, cin, ARTERY, "input", (:left, :center), (10, 0)), (tout, cout, VEIN, "output", side...))
        i = argmax(c)
        Mke.scatter!(ax, [t[i]], [c[i]]; color = col, markersize = 11, strokecolor = SURFACE, strokewidth = 2)
        Mke.text!(ax, t[i], c[i]; text = @sprintf("%s peak: %.1f mg/mL at %.0f s", uppercasefirst(what), c[i], t[i]), color = INK_2, fontsize = 15, align, offset)
    end
    Mke.xlims!(ax, 0, 120)
    Mke.ylims!(ax, 0, nothing)
    Mke.axislegend(ax; position = :rc, framevisible = false, labelcolor = INK_2, labelsize = 15, patchsize = (24, 3), backgroundcolor = :transparent)
    Mke.Label(
        fig[2, 1], "Injection into the right arm: 35 mL of contrast (350 mg iodine per mL) at 5 mL/s, then 35 mL of a 30:70 contrast:saline chaser. After the first pass the iodine mixes into the whole body's blood and recirculates.";
        fontsize = 14, color = INK_2, halign = :left, justification = :left, word_wrap = true, tellwidth = false,
    )
    Mke.resize_to_layout!(fig)
    Mke.save(path, fig; px_per_unit = 2)
    return path
end

"When the iodine reaches each part of the leg: the main arteries by their arrival, each tissue patch by the arrival at its inlet."
function fig_arrival(V, Dp, medians, t_arrive, path; tall = 860, tmax = 40.0)
    B = leg_bounds(V)
    wide = round(Int, tall * (B.x1 - B.x0) / (B.z1 - B.z0))
    pxmm = tall / (B.z1 - B.z0)
    cmap = [Mke.RGBf(hexrgb(h)...) for h in BLUES[4:end]] # from step 250: the earliest vessels still show on the light leg
    fig = Mke.Figure(backgroundcolor = SURFACE, fontsize = 16, figure_padding = 24)
    Mke.Label(fig[1, 1:3], "When the iodine arrives, at rest"; fontsize = 22, font = :bold, color = INK, halign = :left, tellwidth = false)
    ax = leg_axis(fig[2, 1], V, B, wide, tall)
    Mke.scatter!(ax, V.tis.x, V.tis.z; color = clamp.(Dp, 0, tmax), colormap = cmap, colorrange = (0, tmax), markersize = clamp.(2.2 .* cbrt.(V.tis.ml), 3, 8), strokewidth = 0)
    oa = sortperm(V.art.r)
    Mke.linesegments!(ax, segpoints(V.art, oa); color = repeat(clamp.(V.art.tau[oa], 0, tmax); inner = 2), colormap = cmap, colorrange = (0, tmax), linewidth = repeat(width.(V.art.r[oa], pxmm); inner = 2))
    Mke.Colorbar(
        fig[2, 2]; colormap = cmap, limits = (0, tmax), label = "Arrival after the iodine reaches the top of the leg (s)",
        labelcolor = INK_2, ticklabelcolor = INK_2, tickcolor = BASELINE, spinewidth = 0, labelsize = 17, ticklabelsize = 15,
        width = 14, height = Mke.Relative(0.55), ticks = 0:10:tmax,
    )
    g = fig[2, 3] = Mke.GridLayout(valign = :top, tellheight = false)
    Mke.Label(g[1, 1:2], "Median arrival\nat the patches"; fontsize = 17, color = INK_2, halign = :left, justification = :left)
    for (i, (name, m)) in enumerate(medians)
        Mke.Label(g[i + 1, 1], name; fontsize = 17, color = INK_2, halign = :left)
        Mke.Label(g[i + 1, 2], @sprintf("%.0f s", m); fontsize = 17, color = INK, halign = :right)
    end
    Mke.rowgap!(g, 1, 14)
    Mke.colgap!(g, 1, 18)
    Mke.Label(
        fig[3, 1:3], @sprintf("Right leg, front view. Lines: the main arteries, by when the iodine reaches them. Dots: the tissue patches, by when it reaches the start of each one's own arteries. The iodine reaches the top of the leg %.0f s after the injection starts.", t_arrive);
        fontsize = 14, color = INK_2, halign = :left, justification = :left, word_wrap = true, tellwidth = false,
    )
    Mke.colgap!(fig.layout, 24)
    Mke.resize_to_layout!(fig)
    Mke.save(path, fig; px_per_unit = 2)
    return path
end

# ---------------------------------------------------------------------------------------------------------------------
# one page with the figures

using Base64

"Width and height (px) of a PNG, from its header."
png_size(path) = open(path) do io
    b = read(io, 24)
    return Int(ntoh(reinterpret(UInt32, b[17:20])[1])), Int(ntoh(reinterpret(UInt32, b[21:24])[1]))
end

_esc(s) = replace(s, "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")

"""
One self-contained page with the figures `items` (id, heading, alt text, png path) embedded at full resolution, each shown
at half its pixel size (they are drawn at 2×), so copying an image from the page gives the full-resolution one.
"""
function write_figures_page(path, title, sub, items; foot = "")
    io = IOBuffer()
    print(io, """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>$(_esc(title))</title>
<style>
:root { color-scheme: light; --page: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink2: #52514e; --link: #2a78d6; --border: rgba(11,11,11,0.10); }
* { box-sizing: border-box; }
html, body { margin: 0; background: var(--page); color: var(--ink); }
body { font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1160px; margin: 0 auto; padding: 24px 16px 40px; }
h1 { font-size: 22px; font-weight: 600; margin: 0 0 6px; }
.sub { color: var(--ink2); margin: 0 0 14px; }
nav { display: flex; flex-wrap: wrap; gap: 6px 18px; margin: 0 0 20px; }
nav a, .foot a { color: var(--link); text-decoration: none; }
nav a:hover, .foot a:hover { text-decoration: underline; }
figure { margin: 0 0 24px; background: var(--surface); border: 1px solid var(--border); border-radius: 12px; padding: 16px; }
figure h2 { font-size: 16px; font-weight: 600; margin: 0 0 12px; }
figure img { display: block; max-width: 100%; height: auto; margin: 0 auto; }
figcaption { color: var(--ink2); font-size: 13px; margin-top: 10px; }
.foot { color: var(--ink2); font-size: 13px; }
</style>
</head>
<body>
<main>
<h1>$(_esc(title))</h1>
<p class="sub">$(_esc(sub))</p>
<nav>""")
    for (i, (id, heading, _, _)) in enumerate(items)
        print(io, "<a href=\"#", id, "\">", i, ". ", _esc(heading), "</a>")
    end
    print(io, "</nav>\n")
    for (i, (id, heading, alt, png)) in enumerate(items)
        w, h = png_size(png)
        print(io, "<figure id=\"", id, "\"><h2>", i, ". ", _esc(heading), "</h2>")
        print(io, "<img src=\"data:image/png;base64,", base64encode(read(png)), "\" width=\"", w ÷ 2, "\" height=\"", h ÷ 2, "\" alt=\"", _esc(alt), "\">")
        print(io, "<figcaption>", _esc(basename(png)), ", ", w, " × ", h, " px</figcaption></figure>\n")
    end
    isempty(foot) || print(io, "<p class=\"foot\">", foot, "</p>\n")
    print(io, "</main>\n</body>\n</html>\n")
    write(path, String(take!(io)))
    return path
end
