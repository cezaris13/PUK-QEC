import argparse
import csv
import hashlib
import re
from pathlib import Path

import cairosvg
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from floquet import (HEX_CORNERS, QUBIT_CORNERS, hex_centers, hex_positions, memory_circuit, period,
                     qubit_kinds, round_circuit)

# Categorical slots in fixed order: one per distance, or per plaquette colour (dataviz reference palette).
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
# How each kind of qubit is drawn, as in SpinQEC's layout plots: fill and legend label. All get a dark outline.
QUBIT_KINDS = {"data": ("#333333", "data"), "main": ("white", "ancilla: syndrome"),
               "reff": ("#cfcfcf", "ancilla: reference")}


def svg_png(svg: str, png: Path) -> None:
    cairosvg.svg2png(bytestring=svg.encode(), write_to=str(png), scale=2, background_color="white")


def timeslices_png(distance: int, rounds: int, png: Path) -> None:
    svg_png(str(memory_circuit(distance, rounds).diagram("timeslice-svg")), png)


def window(points: list[complex], torus_period: complex) -> complex:
    """Top-left corner of the one-torus-period window centred on `points`."""
    xs, ys = [p.real for p in points], [p.imag for p in points]
    return complex(min(xs) - (torus_period.real - (max(xs) - min(xs))) / 2,
                   min(ys) - (torus_period.imag - (max(ys) - min(ys))) / 2)


def qubit_labels(distance: int) -> dict[complex, str]:
    """Qubit -> its label: "D<i>" for data qubit i, "S<k>" and "R<k>" for the syndrome and reference
    ancillas of edge k (its place in `edge_list`). Stim indices: D<i> is i, S<k> is N + 2k and R<k> is
    N + 2k + 1, with N the number of data qubits."""
    kinds = qubit_kinds(distance)
    n_data = list(kinds.values()).count("data")
    return {q: {"data": "D", "main": "S", "reff": "R"}[kind] + str(i if kind == "data" else (i - n_data) // 2)
            for i, (q, kind) in enumerate(kinds.items())}


def draw_layout(ax: plt.Axes, distance: int, hex_view: bool = False, numbers: bool = False) -> None:
    """The initial qubit layout: plaquettes in their colour, data qubits on the corners, and on every edge
    between them a syndrome ancilla and its reference. Brick-wall rectangles as the circuits use them, or
    regular hexagons with `hex_view`; `numbers` labels every qubit (see `qubit_labels`). Shows one torus
    period, plaquettes wrapped into it, so edges crossing the boundary still have their ancillas on a side."""
    kinds = qubit_kinds(distance)
    pos = hex_positions(distance) if hex_view else {q: q for q in kinds}
    p = period(distance)
    for c, colour in hex_centers(distance).items():
        for shift in [a * p.real + b * p.imag * 1j for a in (-1, 0, 1) for b in (-1, 0, 1)]:
            corners = [c + shift + d for d in (HEX_CORNERS if hex_view else QUBIT_CORNERS)]
            ax.fill([q.real for q in corners], [q.imag for q in corners], color=COLORS[colour], alpha=0.25,
                    edgecolor="0.4")
    for kind, (fill, label) in QUBIT_KINDS.items():
        qs = [pos[q] for q, k in kinds.items() if k == kind]
        ax.scatter([q.real for q in qs], [q.imag for q in qs], s=60 if kind == "data" else 36, facecolor=fill,
                   edgecolor="#333333", linewidths=1.2, label=f"{label} ({len(qs)})", zorder=3)
    if numbers:
        for q, text in qubit_labels(distance).items():
            ax.annotate(text, (pos[q].real, pos[q].imag), xytext=(4, 3), textcoords="offset points",
                        fontsize=6.5, color="#222222", fontweight="bold" if kinds[q] == "data" else "normal",
                        zorder=4)
    corner = window(list(pos.values()), p)
    ax.set_xlim(corner.real, corner.real + p.real)
    ax.set_ylim(corner.imag + p.imag, corner.imag)  # y grows downward, as in stim's timeslice diagrams
    ax.set_aspect("equal")
    ax.set_title(f"{len(kinds)} qubits, " + ("regular hexagons (drawing only)" if hex_view
                                            else "brick wall (the circuits' coordinates)"), fontsize=10)


def layout_figure(distance: int, hex_view: bool = False, numbers: bool = False) -> plt.Figure:
    """One view of the layout with its legend and, with `numbers`, what the labels mean."""
    fig, ax = plt.subplots(figsize=(6, 8))
    draw_layout(ax, distance, hex_view, numbers)
    handles, _ = ax.get_legend_handles_labels()
    handles += [Patch(facecolor=COLORS[c], alpha=0.25, edgecolor="0.4",
                      label=f"colour {c}: sub-round {c} checks\nthe edges linking two of these")
                for c in range(3)]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False, fontsize=9)
    if numbers:
        n_data = list(qubit_kinds(distance).values()).count("data")
        ax.text(1.04, 0.5,
                "Labels\n"
                f"D$i$  data qubit $i$, stim index $i$\n"
                f"S$k$  syndrome ancilla of edge $k$: next to its\n"
                f"       first data qubit, does the calculation,\n"
                f"       stim index {n_data} + 2$k$\n"
                f"R$k$  reference ancilla of edge $k$: next to its\n"
                f"       second data qubit, held in a known state,\n"
                f"       stim index {n_data} + 2$k$ + 1\n"
                f"$k$ is the edge's place in edge_list, so S$k$ and R$k$\n"
                f"are one spin-blockade pair.",
                transform=ax.transAxes, va="top", fontsize=9, color="#333333", linespacing=1.4)
    return fig


def layout_png(distance: int, png: Path, hex_view: bool = False, numbers: bool = False) -> None:
    fig = layout_figure(distance, hex_view, numbers)
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def color_hexes(svg: str, centers: dict[complex, int], corners: list[complex], i2pos: dict[int, complex],
                torus_period: complex) -> str:
    """Draws colour-filled plaquettes underneath a stim timeslice-svg diagram."""
    # ponytail: stim has no colour option for timeslices, so recover coordinate -> pixel from its qubit
    # dot ids. Breaks if stim changes its svg ids; the assert below says so.
    panels = [tuple(map(float, r)) for r in re.findall(
        r'id="tick_border:[^"]*" x="([-\d.]+)" y="([-\d.]+)" width="([-\d.]+)" height="([-\d.]+)"', svg)]
    # a circuit without TICKs is one borderless panel: the whole svg
    panels = panels or [(0, 0, *map(float, re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg).groups()))]
    # dot ids only carry the panel column, so find each dot's panel by position
    dots = [(i2pos[int(q)], next(k for k, (px, py, w, h) in enumerate(panels)
                                 if px <= float(x) <= px + w and py <= float(y) <= py + h), float(x), float(y))
            for q, x, y in re.findall(r'id="qubit_dot:(\d+):[^"]*" cx="([-\d.]+)" cy="([-\d.]+)"', svg)]
    assert dots and panels, "stim changed its timeslice svg format"
    xs = sorted((p.real, x) for p, t, x, _ in dots if t == 0)
    ys = sorted((p.imag, y) for p, t, _, y in dots if t == 0)
    sx = (xs[-1][1] - xs[0][1]) / (xs[-1][0] - xs[0][0])
    sy = (ys[-1][1] - ys[0][1]) / (ys[-1][0] - ys[0][0])
    origin = {t: complex(x - sx * p.real, y - sy * p.imag) for p, t, x, y in dots}
    # one torus period around the qubits, so wrapped hexes show up once on each side
    corner = window([p for p, t, _, _ in dots if t == 0], torus_period)
    left, top = corner.real, corner.imag
    shifts = [a * torus_period.real + b * torus_period.imag * 1j for a in (-1, 0, 1) for b in (-1, 0, 1)]

    # ids are page-wide once several SVGs share an HTML page (a notebook), so make them unique per diagram
    tag = hashlib.sha1(svg.encode()).hexdigest()[:10]
    out = []
    for t, o in origin.items():
        out.append(f'<clipPath id="hexclip{tag}_{t}"><rect x="{o.real + sx * left}" y="{o.imag + sy * top}" '
                   f'width="{sx * torus_period.real}" height="{sy * torus_period.imag}"/></clipPath>'
                   f'<g clip-path="url(#hexclip{tag}_{t})">')
        for h, colour in centers.items():
            for s in shifts:
                pts = " ".join(f"{o.real + sx * (h + s + d).real},{o.imag + sy * (h + s + d).imag}" for d in corners)
                out.append(f'<polygon points="{pts}" fill="{COLORS[colour]}" fill-opacity="0.3" stroke="#666"/>')
        out.append("</g>")
    return svg.replace("\n", "\n" + "\n".join(out) + "\n", 1)  # right after <svg>, i.e. below everything


def style_qubits(svg: str, kinds: list[str]) -> str:
    """Stim draws every qubit as the same small dot: redraw each as its kind (`kinds[i]` for qubit i)."""
    def dot(m: re.Match[str]) -> str:
        fill = QUBIT_KINDS[kinds[int(m[1])]][0]
        return (f'<circle id="qubit_dot:{m[1]}:{m[2]}" cx="{m[3]}" cy="{m[4]}" r="5" fill="{fill}" '
                f'stroke="#333333" stroke-width="1.5"/>')
    styled, n = re.subn(r'<circle id="qubit_dot:(\d+):([^"]*)" cx="([-\d.]+)" cy="([-\d.]+)" r="2" '
                        r'stroke="none" fill="black"/>', dot, svg)
    assert n, "stim changed its timeslice svg format"
    return styled


def round_svg(distance: int, r: int, hex_view: bool) -> str:
    """Timeslice view of sub-round r, plaquettes coloured underneath; brick-wall or regular hexagons."""
    circuit = round_circuit(distance, r, hex_view)
    i2pos = {i: complex(*xy) for i, xy in circuit.get_final_qubit_coordinates().items()}
    svg = color_hexes(str(circuit.diagram("timeslice-svg")), hex_centers(distance),
                      HEX_CORNERS if hex_view else QUBIT_CORNERS, i2pos, period(distance))
    return style_qubits(svg, list(qubit_kinds(distance).values()))


def honeycomb_pngs(distance: int, out: Path, numbers: bool = False) -> None:
    """layout.png and layout_hex.png (qubits labelled with `numbers`), then round{r}.png and
    round{r}_hex.png for each of the 3 sub-rounds."""
    layout_png(distance, out / "layout.png", False, numbers)
    layout_png(distance, out / "layout_hex.png", True, numbers)
    for r in range(6):
        svg_png(round_svg(distance, r, False), out / f"round{r}.png")
        svg_png(round_svg(distance, r, True), out / f"round{r}_hex.png")


def draw_ler(ax: plt.Axes, runs: list[list[dict]]) -> None:
    """One log-log line of LER vs p per run (a run = one distance's rows); the crossing is the threshold."""
    runs = sorted(runs, key=lambda rows: int(rows[0]["distance"]))
    if len(runs) > len(COLORS):
        raise ValueError(f"{len(runs)} distances, only {len(COLORS)} distinguishable colours; plot fewer")
    for color, rows in zip(COLORS, runs):
        label = f"d = {rows[0]['distance']}"
        rows = [r for r in rows if int(r["errors"])]  # zero errors has no place on a log axis
        ax.loglog([float(r["p"]) for r in rows], [float(r["ler"]) for r in rows], "o-",
                  color=color, lw=1.5, ms=5, label=label)
    ax.set(xlabel="physical error rate p", ylabel="logical error rate",
           title=f"{runs[0][0]['noise']} noise, η = {float(runs[0][0]['eta']):g}")
    ax.grid(which="major", color="#e4e4e0", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)


def ler_png(csv_paths: list[Path], png: Path) -> None:
    runs = []
    for path in csv_paths:
        with open(path, newline="") as f:
            runs.append(list(csv.DictReader(f)))
    fig, ax = plt.subplots(figsize=(6, 4.5))
    draw_ler(ax, runs)
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="figure", required=True)
    ts = sub.add_parser("timeslices")
    ts.add_argument("--distance", type=int, default=3)
    ts.add_argument("--rounds", type=int, default=3)
    ts.add_argument("--png", type=Path, required=True)
    hc = sub.add_parser("honeycomb", help="qubit layout and coloured sub-round timeslices, into --out")
    hc.add_argument("--distance", type=int, default=1)
    hc.add_argument("--out", type=Path, required=True)
    hc.add_argument("--numbers", action="store_true",
                    help="label every qubit in the layouts: Di data, Sk / Rk syndrome / reference ancilla")
    ler = sub.add_parser("ler")
    ler.add_argument("csvs", type=Path, nargs="+")
    ler.add_argument("--png", type=Path, required=True)
    args = parser.parse_args()

    if args.figure == "honeycomb":
        args.out.mkdir(parents=True, exist_ok=True)
        honeycomb_pngs(args.distance, args.out, args.numbers)
        print(f"wrote {args.out}/")
        return
    args.png.parent.mkdir(parents=True, exist_ok=True)
    if args.figure == "timeslices":
        timeslices_png(args.distance, args.rounds, args.png)
    else:
        ler_png(args.csvs, args.png)
    print(f"wrote {args.png}")


if __name__ == "__main__":
    main()
