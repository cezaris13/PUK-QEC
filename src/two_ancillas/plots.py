import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from drawing import COLORS, QUBIT_KINDS, color_hexes, drawing_dir, label_qubits, style_qubits, svg_png, window
from floquet import HEX_CORNERS, QUBIT_CORNERS, hex_centers, memory_circuit, period
from pairs import positions, qubit_kinds, round_circuit


def layout_figure(distance: int, hex_view: bool = False, numbers: bool = False) -> plt.Figure:
    """One view of the layout with its legend and, with `numbers`, what the labels mean."""
    fig, ax = plt.subplots(figsize=(6, 8))
    _draw_layout(ax, distance, hex_view, numbers)
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


def round_svg(distance: int, r: int, hex_view: bool, numbers: bool = False, code: str = "css") -> str:
    """Timeslice view of sub-round r, plaquettes coloured underneath; brick-wall or regular hexagons.
    `numbers` labels every qubit as in the layouts (see `_qubit_labels`); `code` "x3z3" draws the X3Z3 step."""
    circuit = round_circuit(distance, r, hex_view, code)
    i2pos = {i: complex(*xy) for i, xy in circuit.get_final_qubit_coordinates().items()}
    svg = color_hexes(str(circuit.diagram("timeslice-svg")), hex_centers(distance),
                      HEX_CORNERS if hex_view else QUBIT_CORNERS, i2pos, period(distance))
    svg = style_qubits(svg, list(qubit_kinds(distance).values()))
    return label_qubits(svg, list(_qubit_labels(distance).values())) if numbers else svg


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
    # Below this line the logical qubit beats a bare physical one.
    ax.axline((1e-3, 1e-3), (1e-2, 1e-2), color="0.5", ls="--", lw=1, label="LER = p")
    ax.set(xlabel="physical error rate p", ylabel="logical error rate",
           title=f"{runs[0][0]['noise']} noise, η = {float(runs[0][0]['eta']):g}")
    ax.grid(which="major", color="#e4e4e0", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)


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
    ler = sub.add_parser("ler")
    ler.add_argument("csvs", type=Path, nargs="+")
    ler.add_argument("--png", type=Path, required=True)
    args = parser.parse_args()

    if args.figure == "honeycomb":
        _honeycomb_pngs(args.distance, args.out)
        print(f"wrote {args.out}/")
        return
    args.png.parent.mkdir(parents=True, exist_ok=True)
    if args.figure == "timeslices":
        _timeslices_png(args.distance, args.rounds, args.png)
    else:
        _ler_png(args.csvs, args.png)
    print(f"wrote {args.png}")


def _timeslices_png(distance: int, rounds: int, png: Path) -> None:
    svg_png(str(memory_circuit(distance, rounds, round_circuit).diagram("timeslice-svg")), png)


def _qubit_labels(distance: int) -> dict[complex, str]:
    """Qubit -> its label: "D<i>" for data qubit i, "S<k>" and "R<k>" for the syndrome and reference
    ancillas of edge k (its place in `edge_list`). Stim indices: D<i> is i, S<k> is N + 2k and R<k> is
    N + 2k + 1, with N the number of data qubits."""
    kinds = qubit_kinds(distance)
    n_data = list(kinds.values()).count("data")
    return {q: {"data": "D", "main": "S", "reff": "R"}[kind] + str(i if kind == "data" else (i - n_data) // 2)
            for i, (q, kind) in enumerate(kinds.items())}


def _draw_layout(ax: plt.Axes, distance: int, hex_view: bool = False, numbers: bool = False) -> None:
    """The initial qubit layout: plaquettes in their colour, data qubits on the corners, and on every edge
    between them a syndrome ancilla and its reference. Brick-wall rectangles as the circuits use them, or
    regular hexagons with `hex_view`; `numbers` labels every qubit (see `_qubit_labels`). Shows one torus
    period, plaquettes wrapped into it, so edges crossing the boundary still have their ancillas on a side."""
    kinds = qubit_kinds(distance)
    pos = positions(distance) if hex_view else {q: q for q in kinds}
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
        for q, text in _qubit_labels(distance).items():
            ax.annotate(text, (pos[q].real, pos[q].imag), xytext=(4, 3), textcoords="offset points",
                        fontsize=6.5, color="#222222", fontweight="bold" if kinds[q] == "data" else "normal",
                        zorder=4)
    corner = window(list(pos.values()), p)
    ax.set_xlim(corner.real, corner.real + p.real)
    ax.set_ylim(corner.imag + p.imag, corner.imag)  # y grows downward, as in stim's timeslice diagrams
    ax.set_aspect("equal")
    ax.set_title(f"{len(kinds)} qubits, " + ("regular hexagons (drawing only)" if hex_view
                                            else "brick wall (the circuits' coordinates)"), fontsize=10)


def _layout_png(distance: int, png: Path, hex_view: bool = False, numbers: bool = False) -> None:
    fig = layout_figure(distance, hex_view, numbers)
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _honeycomb_pngs(distance: int, out: Path) -> None:
    """layout.png and round{r}.png for each of the 6 steps, in every `drawing_dir` of `out`."""
    for hex_view in (False, True):
        for numbers in (False, True):
            folder = drawing_dir(out, hex_view, numbers)
            _layout_png(distance, folder / "layout.png", hex_view, numbers)
            for r in range(6):
                svg_png(round_svg(distance, r, hex_view, numbers), folder / f"round{r}.png")


def _ler_png(csv_paths: list[Path], png: Path) -> None:
    runs = []
    for path in csv_paths:
        with open(path, newline="") as f:
            runs.append(list(csv.DictReader(f)))
    fig, ax = plt.subplots(figsize=(6, 4.5))
    draw_ler(ax, runs)
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
