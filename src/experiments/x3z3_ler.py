"""The X3Z3 Floquet code against the CSS Floquet code, for every readout scheme, under biased noise: logical error
rate vs p at each eta, as Setiawan & McLauchlan (arXiv:2411.04974) for these circuits. x3z3_bias.py draws threshold
and p_L vs eta from the same sweep.

    .venv/bin/python src/experiments/x3z3_ler.py --noise spin --distances 3 5 --etas 1 1000
    .venv/bin/python src/experiments/x3z3_ler.py --plot-only --out <out>    # redraw from <out>.csv

Every (code, scheme, distance, eta, p) runs the Z and the X memory experiment (floquet.memory_circuit, d rounds,
noise.NOISE_MODELS[--noise] at (p, eta)), decoded by matching with the exact noise (memory.logical_errors).
Their logical error rates combine as the paper's Eq. (9), p_L = 1 - (1 - p_Z)(1 - p_X): either logical failing.
The 95% Clopper-Pearson interval of each basis combines the same way. The p grid is --p-min..--p-max.

    <out>_ler_eta<eta>.png  p_L vs p, one panel per scheme, both codes: colour = distance, X3Z3 filled, CSS hollow
                            and faded, shaded 95% interval.

Each finished point is appended to <out>.csv at once, and a rerun with the same settings skips the points already
there, so a cut-short run picks up where it stopped. x3z3_bias.py with the same settings writes the same CSV, so
each reuses the other's points. Points with no logical errors are left out of the plots.
"""
import argparse
import csv
import os
from functools import lru_cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from floquet import CODES, memory_circuit
from memory import logical_errors, run_points
from noise import NOISE_MODELS, add_noise
from physical_to_logical import DISTANCE_COLOUR, SCHEMES, STYLES, clopper_pearson

ROOT = Path(__file__).resolve().parents[2]
# As physical_to_logical.py and footprint.py: colour = distance (DISTANCE_COLOUR), line style and marker = scheme
# (STYLES). The code is the marker fill: X3Z3 filled, CSS hollow and faded.
CODE_LABEL = {"x3z3": "X$^3$Z$^3$", "css": "CSS"}
CODE_ALPHA = {"x3z3": 1.0, "css": 0.55}


def main() -> None:
    args = sweep_parser(__doc__).parse_args()
    rows = sweep(args)
    for eta in sorted({float(r["eta"]) for r in rows}):
        png = Path(f"{args.out}_ler_eta{eta:g}.png")
        _ler_figure(rows, eta, png)
        print(f"wrote {png}")


def sweep_parser(doc: str) -> argparse.ArgumentParser:
    """The settings of the (code, scheme, distance, eta, p) sweep, shared with x3z3_bias.py."""
    parser = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--codes", nargs="+", choices=CODES, default=list(CODES))
    parser.add_argument("--schemes", nargs="+", choices=SCHEMES, default=list(SCHEMES))
    parser.add_argument("--noise", choices=NOISE_MODELS, default="spin")
    parser.add_argument("--distances", type=int, nargs="+", default=[3, 5])
    parser.add_argument("--etas", type=float, nargs="+", default=[0.5, 1, 3, 10, 30, 100, 1000])
    parser.add_argument("--p-min", type=float, default=2e-4)
    parser.add_argument("--p-max", type=float, default=1.5e-2)
    parser.add_argument("--num", type=int, default=12, help="p values, log-spaced")
    parser.add_argument("--shots", type=int, default=20_000, help="most shots per point and basis")
    parser.add_argument("--max-errors", type=int, default=500, help="stop a point at this many logical errors")
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--out", type=Path, help="writes <out>.csv and <out>_*.png")
    parser.add_argument("--plot-only", action="store_true", help="redraw from <out>.csv")
    return parser


def sweep(args, extra_ps: set[float] = frozenset()) -> list[dict]:
    """Every point of the grid `args` sets, plus `extra_ps`, in <args.out>.csv (args.out defaults from the
    settings); simulates nothing with --plot-only. The CSV's rows, as text."""
    if args.out is None:
        args.out = ROOT / "results" / "single_ancilla" / "x3z3" / (
            f"{args.noise}_d{'-'.join(map(str, args.distances))}_p{args.p_min:g}-{args.p_max:g}x{args.num}"
            f"_shots{args.shots}")
    csv_path = args.out.with_suffix(".csv")
    if not args.plot_only:
        ps = sorted(set(np.geomspace(args.p_min, args.p_max, args.num).tolist()) | set(extra_ps))
        tasks = [dict(code=c, scheme=s, distance=d, eta=eta, p=p, basis=b, noise=args.noise, shots_max=args.shots,
                      max_errors=args.max_errors)
                 for c in args.codes for s in args.schemes for d in args.distances for eta in args.etas
                 for p in ps for b in "ZX"]
        # slowest first (largest distance, then highest p), so the last workers aren't left with one big job
        run_points(csv_path, sorted(tasks, key=lambda t: (-t["distance"], -t["p"])), _point, args.workers)
    return list(csv.DictReader(csv_path.open()))


def combined(rows: list[dict]) -> dict:
    """{(code, scheme, distance, eta, p): (p_L, lower, upper, qubits)}: the Z and X memories as Eq. (9),
    each end of the interval combined the same way; points where neither basis saw an error are left out."""
    by = {}
    for r in rows:
        by.setdefault((r["code"], r["scheme"], int(r["distance"]), float(r["eta"]), float(r["p"])), {})[r["basis"]] = r
    out = {}
    for k, bases in by.items():
        if len(bases) < 2 or not any(int(b["errors"]) for b in bases.values()):
            continue
        z, x = (_interval(int(bases[b]["errors"]), int(bases[b]["shots"])) for b in "ZX")
        out[k] = (*(1 - (1 - a) * (1 - b) for a, b in zip(z, x)), int(bases["Z"]["qubits"]))
    return out


def style_axes(ax) -> None:
    """physical_to_logical.py's plot styling."""
    ax.grid(which="major", color="#e4e4e0", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)


def blocks(ax, titled: list[tuple[str, list[Line2D]]], top: float = 1.0) -> None:
    """Legend blocks stacked right of `ax`, as physical_to_logical.legend_blocks: (title, handles), empty ones
    skipped."""
    y = top
    for title, handles in titled:
        if handles:
            ax.add_artist(ax.legend(handles=handles, title=title, fontsize=7, title_fontsize=8, loc="upper left",
                                    bbox_to_anchor=(1.02, y), frameon=False))
            y -= 0.08 + 0.055 * len(handles)


def save(fig, png: Path) -> None:
    """Save cropped to its content: bbox_inches="tight" alone leaves out legends added with add_artist, and then
    the suptitle too once extra artists are given."""
    extra = [a for ax in fig.axes for a in ax.artists] + [t for t in [fig._suptitle] if t is not None]
    fig.savefig(png, dpi=200, bbox_inches="tight", bbox_extra_artists=extra)
    plt.close(fig)


@lru_cache(maxsize=None)
def _circuit(code: str, scheme: str, distance: int, basis: str):
    """The noiseless memory circuit, d rounds: once per worker."""
    return memory_circuit(distance, distance, SCHEMES[scheme], basis, code)


def _point(task: dict) -> dict:
    c = _circuit(task["code"], task["scheme"], task["distance"], task["basis"])
    errors, shots = logical_errors(add_noise(c, NOISE_MODELS[task["noise"]](task["p"], task["eta"])),
                                   task["shots_max"], task["max_errors"])
    return dict(task, qubits=c.num_qubits, shots=shots, errors=errors)


def _interval(errors: int, shots: int) -> tuple[float, float, float]:
    """(rate, lower, upper), 95% Clopper-Pearson."""
    return (errors / shots, *clopper_pearson(errors, shots))


def _mark(code: str, colour: str) -> dict:
    """Marker fill and line opacity for `code`: X3Z3 filled, CSS hollow and faded."""
    return dict(color=colour, mfc=colour if code == "x3z3" else "white", mec=colour, alpha=CODE_ALPHA[code])


def _ler_figure(rows: list[dict], eta: float, png: Path) -> None:
    """physical_to_logical.py's LER vs p, one panel per scheme with both codes: colour = distance, the scheme's
    line style and marker, X3Z3 filled, CSS hollow and faded, shaded 95% Clopper-Pearson band."""
    pl = combined(rows)
    schemes = [s for s in SCHEMES if any(r["scheme"] == s for r in rows)]
    distances = sorted({int(r["distance"]) for r in rows})
    fig, axes = plt.subplots(1, len(schemes), figsize=(4 * len(schemes) + 1.5, 4.2), sharey=True, squeeze=False)
    fig.subplots_adjust(wspace=0.08, right=0.88)
    for ax, scheme in zip(axes[0], schemes):
        ls, marker = STYLES[scheme]
        for code in CODES:
            for d in distances:
                pts = sorted((k[4], v) for k, v in pl.items() if k[:4] == (code, scheme, d, eta))
                if not pts:
                    continue
                ps, vs = [p for p, _ in pts], [v for _, v in pts]
                colour = DISTANCE_COLOUR.get(d, "0.3")
                ax.plot(ps, [v[0] for v in vs], ls=ls, marker=marker, ms=4, lw=1.5, **_mark(code, colour))
                ax.fill_between(ps, [v[1] for v in vs], [v[2] for v in vs], color=colour,
                                alpha=0.2 * CODE_ALPHA[code], lw=0)
        ax.axline((1e-3, 1e-3), (1e-2, 1e-2), color="0.5", ls=":", lw=1)
        ax.loglog()
        ax.set(xlabel="physical error rate p", title=scheme)
        style_axes(ax)
    axes[0][0].set_ylabel("logical failure probability $p_L$ (d rounds)")
    codes = [c for c in CODES if any(r["code"] == c for r in rows)]
    blocks(axes[0][-1], [
        ("Schemes", [Line2D([], [], marker=STYLES[s][1], ls=STYLES[s][0], color="0.3", label=s) for s in schemes]),
        ("Codes", [Line2D([], [], marker="o", ls="none", label=CODE_LABEL[c], **_mark(c, "0.3")) for c in codes]),
        ("Distances", [Line2D([], [], marker="s", ls="none", color=DISTANCE_COLOUR.get(d, "0.3"), label=str(d))
                       for d in distances])])
    fig.suptitle(f"{rows[0]['noise']} noise, η = {eta:g}; dotted: LER = p")
    save(fig, png)


if __name__ == "__main__":
    main()
