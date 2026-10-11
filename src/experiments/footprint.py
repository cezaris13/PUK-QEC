"""How many physical qubits each readout scheme needs at one physical error rate: logical error rate vs qubit
count, one line per scheme, extrapolated to a target logical error rate.

    .venv/bin/python src/experiments/footprint.py --p 1e-3 --distances 3 5 --shots 10000000 --target 1e-6 \
        --schemes method_a method_b method_c --out results/single_ancilla/footprint
    .venv/bin/python src/experiments/footprint.py --codes css x3z3 --eta 1000   # X3Z3 vs CSS Floquet code

Each point runs the Z and the X memory experiment (d rounds at distance d) at p, decoded by matching, and combines
them as x3z3_ler.py does, p_L = 1 - (1 - p_Z)(1 - p_X): either logical failing. Both are needed: under biased noise a
CSS Z memory barely sees the dominant Z errors, so alone it would make CSS look far better than it is. The LER falls
exponentially in d, so a straight line through log(LER) vs d (points with at least one logical error) gives
the d where it reaches --target; the qubit count there is the scheme's qubits per d^2 times d^2. Points with no
errors are left out: at low p the largest distances need a lot of shots to show any.

Plotted on square-root-log axes (as the N2E3N2 paper's Figure 5a), where each fit is a straight line since the
qubit count grows as d^2; its slope is how fast the scheme suppresses errors, given as
Lambda = LER(d) / LER(d + 2). Each point is appended to <out>.csv as it finishes (writing.results_csv.run_points),
so a rerun picks up where a run stopped; --plot-only redraws the PNG from <out>.csv. With more than one of
--codes, X3Z3 markers are filled and CSS ones hollow, one fit per (code, scheme). --sources warns if that code
changed since <out>.csv was written.
"""
import argparse
import math
import os
from functools import lru_cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from floquet import CODES, memory_circuit
from memory import logical_errors
from noise import NOISE_MODELS, add_noise
from physical_to_logical import DISTANCE_COLOUR, SCHEMES, STYLES, clopper_pearson, legend_blocks
from writing.results_csv import read_rows, run_points


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--schemes", nargs="+", choices=SCHEMES, default=list(SCHEMES))
    parser.add_argument("--codes", nargs="+", choices=CODES, default=["css"], help="css and/or x3z3")
    parser.add_argument("--distances", type=int, nargs="+", default=[3, 5, 7])
    parser.add_argument("--noise", choices=NOISE_MODELS, default="spin")
    parser.add_argument("--eta", type=float, default=10.0)
    parser.add_argument("--p", type=float, default=1e-3)
    parser.add_argument("--shots", type=int, default=10_000_000, help="most shots per point")
    parser.add_argument("--max-errors", type=int, default=100, help="stop a point at this many logical errors")
    parser.add_argument("--target", type=float, default=1e-6, help="logical error rate to extrapolate to")
    parser.add_argument("--out", type=Path, default=Path("results/single_ancilla/footprint"), help="writes <out>.csv and <out>.png")
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--plot-only", action="store_true", help="redraw <out>.png from <out>.csv")
    parser.add_argument("--sources", nargs="*", default=[], help="warn if any is newer than <out>.csv")
    args = parser.parse_args()

    if args.plot_only:
        rows = read_rows(args.out.with_suffix(".csv"))
    else:
        tasks = [dict(scheme=name, code=code, distance=d, basis=b, rounds=d, noise=args.noise, eta=args.eta, p=args.p,
                      shots_max=args.shots, max_errors=args.max_errors)
                 for code in args.codes for name in args.schemes for d in args.distances for b in "ZX"]
        rows = run_points(args.out.with_suffix(".csv"), tasks, _point, args.workers, sources=args.sources)
    rows = _combined(rows)
    for code, name in dict.fromkeys((r["code"], r["scheme"]) for r in rows):
        fit = _qubits_needed([r for r in rows if (r["code"], r["scheme"]) == (code, name)], args.target)
        print(f"{code:5} {name:12} " + (f"{fit[2]:,.0f} qubits for LER {args.target:g}" if fit else
                                         "no fit: under two distances with errors (more shots or larger p), or the "
                                         "LER rises with d (p above threshold: lower p)"))
    _plot(rows, args.target, args.out.with_suffix(".png"))
    print(f"wrote {args.out}.csv, {args.out}.png")


@lru_cache(maxsize=None)
def _circuit(scheme: str, distance: int, code: str, basis: str):
    """The noiseless memory circuit, d rounds: once per worker."""
    return memory_circuit(distance, distance, SCHEMES[scheme], basis, code)


def _point(task: dict) -> dict:
    """One (scheme, code, distance, basis) point: its qubit count and LER at p, from up to shots_max shots,
    stopping at max_errors logical errors."""
    c = _circuit(task["scheme"], task["distance"], task["code"], task["basis"])
    errors, ran = logical_errors(add_noise(c, NOISE_MODELS[task["noise"]](task["p"], task["eta"])), task["shots_max"],
                                 task["max_errors"])
    return dict(task, qubits=c.num_qubits, shots=ran, errors=errors, ler=errors / ran)


def _combined(rows: list[dict]) -> list[dict]:
    """One row per (code, scheme, distance) with both bases: ler = 1 - (1 - p_Z)(1 - p_X), and lo / hi its 95%
    Clopper-Pearson interval, each basis's ends combined the same way (x3z3_ler.combined)."""
    either = lambda a, b: 1 - (1 - a) * (1 - b)
    by = {}
    for r in rows:
        by.setdefault((r["code"], r["scheme"], r["distance"]), {})[r["basis"]] = r
    out = []
    for bases in by.values():
        if set(bases) != {"Z", "X"}:
            continue
        z, x = bases["Z"], bases["X"]
        (z_lo, z_hi), (x_lo, x_hi) = clopper_pearson(z["errors"], z["shots"]), clopper_pearson(x["errors"], x["shots"])
        out.append(dict(z, basis="ZX", errors=z["errors"] + x["errors"], ler=either(z["ler"], x["ler"]),
                        lo=either(z_lo, x_lo), hi=either(z_hi, x_hi)))
    return out


def _qubits_needed(rows: list[dict], target: float) -> tuple[float, float, float] | None:
    """(slope, intercept, qubits) of the fit log(LER) = intercept + slope * d for one (code, scheme)'s rows, with the
    qubit count where it crosses `target`; None without two distances that saw errors, or if the LER doesn't
    fall with d (p above threshold)."""
    pts = [(r["distance"], math.log(r["ler"])) for r in rows if r["errors"]]
    if len({d for d, _ in pts}) < 2:
        return None
    slope, intercept = np.polyfit(*zip(*pts), 1)
    if slope >= 0:
        return None
    d = (math.log(target) - intercept) / slope
    per_d2 = rows[0]["qubits"] / rows[0]["distance"] ** 2  # every scheme's qubit count is a multiple of d^2
    return slope, intercept, per_d2 * d ** 2


def _plot(rows: list[dict], target: float, png: Path) -> None:
    """LER vs qubits on square-root-log axes, where each scheme's fit is a straight line: points coloured by
    distance with 95% Clopper-Pearson bars, marker = scheme, grey fit in the scheme's line style out to the
    target (star), labelled with the qubits needed and Lambda = LER(d) / LER(d + 2)."""
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    fig.subplots_adjust(left=0.10, right=0.70)
    schemes = list(dict.fromkeys(r["scheme"] for r in rows))
    codes = list(dict.fromkeys(r["code"] for r in rows))
    fits = []
    for code, name in dict.fromkeys((r["code"], r["scheme"]) for r in rows):
        mine = [r for r in rows if (r["code"], r["scheme"]) == (code, name)]
        ls, marker = STYLES[name]
        hollow = len(codes) > 1 and code == "css"
        for r in mine:
            if r["errors"]:
                lo, hi = r["lo"], r["hi"]
                colour = DISTANCE_COLOUR[r["distance"]]
                ax.errorbar(r["qubits"], r["ler"], yerr=[[r["ler"] - lo], [hi - r["ler"]]], marker=marker,
                            color=colour, mfc="white" if hollow else colour, ms=6, lw=1, capsize=2, ls="none")
        fit = _qubits_needed(mine, target)
        if fit:
            slope, intercept, needed = fit
            per_d2 = mine[0]["qubits"] / mine[0]["distance"] ** 2
            ds = np.linspace(min(r["distance"] for r in mine), math.sqrt(needed / per_d2), 50)
            grey = "0.7" if hollow else "0.35"
            ax.plot(per_d2 * ds ** 2, np.exp(intercept + slope * ds), ls=ls, color=grey, lw=1)
            ax.plot(needed, target, marker="*", color=grey, ms=9)
            label = f"{code}, {name}" if len(codes) > 1 else name
            fits.append(Line2D([], [], ls=ls, color=grey,
                               label=f"{label}: {needed:,.0f} qubits, Λ = {math.exp(-2 * slope):.1f}"))
    ax.axhline(target, color="0.5", ls=":", lw=1)
    ax.set_xscale("function", functions=(np.sqrt, np.square))
    # autoscale's margin can push the left edge below 0, where sqrt is NaN and the whole axis vanishes
    ax.set_xlim(left=0.8 * min(r["qubits"] for r in rows))
    ax.xaxis.set_major_locator(MaxNLocator(6))  # sqrt spacing crowds the labels at the right end
    ax.set_yscale("log")
    ax.set(xlabel="Total Physical Qubits", ylabel="Logical Error Rate, Z or X memory (d rounds)")
    ax.set_title(f"{rows[0]['noise']} noise, η = {rows[0]['eta']:g}, p = {rows[0]['p']:g};"
                 f" dotted: target {target:g}", fontsize=10)
    ax.grid(which="both", alpha=0.3, lw=0.5)
    legend_blocks(ax, schemes, sorted({r["distance"] for r in rows}), lines=False)
    if fits:  # none when every (code, scheme) is above threshold
        second = ax.get_legend()
        ax.legend(handles=fits, title="Linefits", fontsize=7, title_fontsize=8, loc="upper left",
                  bbox_to_anchor=(1.02, 0.2), frameon=False)
        ax.add_artist(second)
    fig.savefig(png, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
