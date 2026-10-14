"""Memory experiment for every readout scheme, decoded the same way: logical error rate vs p.

    .venv/bin/python src/experiments/physical_to_logical.py --distances 3 5 --noise spin --eta 10 --shots 10000
    .venv/bin/python src/experiments/physical_to_logical.py --schemes method_a method_b method_c

"pairs" is src/two_ancillas/pairs.py's syndrome + reference ancilla pair on every edge; the rest have one ancilla per
edge: "method_a", "method_b" and "method_c" are hexes.pdf's methods A, B and C (src/single_ancilla/<name>/). All run through
floquet.memory_circuit, so they share every detector and the observable, and only their step circuits (and so
their noise) differ. Decoding is memory.count_logical_errors: Stim's detector error model, decomposed into
a graph, matched by PyMatching (docs/decoding.pdf). One CSV row per (scheme, distance, p), appended as each
finishes (writing.results_csv.run_points), so rerunning the same command after a crash only runs the points
still missing; then a plot. --sources warns if that code changed since <out>.csv was written.
"""
import argparse
import os
from functools import lru_cache, partial
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.stats import beta

import floquet
import pairs
from memory import logical_errors
from noise import NOISE_MODELS, add_noise
from writing.results_csv import read_rows, run_points

import method_a
import method_b
import method_c

SCHEMES = {"pairs": pairs.round_circuit, "method_a": method_a.round_circuit, "method_b": method_b.round_circuit,
           "method_c": method_c.round_circuit, "pairs_bell": partial(pairs.round_circuit, bell=True)}
# line style and marker per scheme; colour is the distance
STYLES = {"pairs": ("-", "o"), "method_a": ("--", "D"), "method_b": ("-.", "^"), "method_c": (":", "v"), "pairs_bell": ((0, (6, 2, 1, 2, 1, 2)), "s")}
# colour per distance: the N2E3N2 paper's Figure 4 green / blue / red for the odd distances used here
DISTANCE_COLOUR = {3: "#2ca02c", 5: "#1f77b4", 7: "#d62728", 9: "#9467bd",
                   4: "#8fce8f", 6: "#8fb4e8", 8: "#f0928f", 10: "#c5b0d5"}


def clopper_pearson(errors: int, shots: int) -> tuple[float, float]:
    """95% Clopper-Pearson interval of errors / shots."""
    low = beta.ppf(0.025, errors, shots - errors + 1) if errors else 0.0
    high = beta.ppf(0.975, errors + 1, shots - errors) if errors < shots else 1.0
    return low, high


def combined(rows: list[dict]) -> list[dict]:
    """One row per (code, scheme, distance, p) with both memory bases: ler = 1 - (1 - p_Z)(1 - p_X), either logical
    failing, and lo / hi its 95% Clopper-Pearson interval, each basis's ends combined the same way. Rows from
    before the code / basis columns (CSS Z memory only) pass through as they are."""
    either = lambda a, b: 1 - (1 - a) * (1 - b)
    by, out = {}, []
    for r in rows:
        if "basis" not in r:
            lo, hi = clopper_pearson(r["errors"], r["shots"])
            out.append(dict(r, code="css", basis="Z", lo=lo, hi=hi))
            continue
        by.setdefault((r["code"], r["scheme"], r["distance"], r["p"]), {})[r["basis"]] = r
    for bases in by.values():
        if set(bases) != {"Z", "X"}:
            continue
        z, x = bases["Z"], bases["X"]
        (z_lo, z_hi), (x_lo, x_hi) = clopper_pearson(z["errors"], z["shots"]), clopper_pearson(x["errors"], x["shots"])
        out.append(dict(z, basis="ZX", errors=z["errors"] + x["errors"], ler=either(z["ler"], x["ler"]),
                        lo=either(z_lo, x_lo), hi=either(z_hi, x_hi)))
    return out


def legend_blocks(ax, schemes: list[str], distances: list[int], lines: bool) -> None:
    """Schemes (marker, and line style if `lines`) and Distances (colour) as two legend blocks right of `ax`."""
    first = ax.legend(handles=[Line2D([], [], marker=STYLES[n][1], ls=STYLES[n][0] if lines else "none",
                                      color="0.3", label=n) for n in schemes],
                      title="Schemes", fontsize=7, title_fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1),
                      frameon=False)
    ax.add_artist(first)
    ax.legend(handles=[Line2D([], [], marker="s", ls="none", color=DISTANCE_COLOUR[d], label=str(d))
                       for d in distances], title="Distances", fontsize=7, title_fontsize=8, loc="upper left",
              bbox_to_anchor=(1.02, 0.45), frameon=False)


def dem_stats(distance: int, noise: str, eta: float, p: float, schemes: list[str]) -> dict:
    """What the decoder gets from each scheme: circuit size, how many DEM mechanisms are hyperedges (3+
    detectors, split up by decomposition) and the shortest undetected logical error (graphlike, an upper
    bound on the circuit distance)."""
    out = {}
    for name in schemes:
        step = SCHEMES[name]
        circuit = floquet.memory_circuit(distance, distance, step)
        noisy = add_noise(circuit, NOISE_MODELS[noise](p, eta))
        dem = noisy.detector_error_model(approximate_disjoint_errors=True)  # raises if a detector is random
        sizes = [sum(t.is_relative_detector_id() for t in e.targets_copy()) for e in dem.flattened()
                 if e.type == "error"]
        out[name] = dict(qubits=circuit.num_qubits, ticks=circuit.num_ticks, detectors=circuit.num_detectors,
                         mechanisms=len(sizes), hyperedges=sum(s > 2 for s in sizes),
                         distance=len(noisy.shortest_graphlike_error()))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--distances", type=int, nargs="+", default=[3, 5])
    parser.add_argument("--schemes", nargs="+", choices=SCHEMES, default=list(SCHEMES))
    parser.add_argument("--codes", nargs="+", choices=floquet.CODES, default=["css"],
                        help="css and/or x3z3; each runs the Z and the X memory, combined in the plot")
    parser.add_argument("--noise", choices=NOISE_MODELS, default="spin")
    parser.add_argument("--eta", type=float, default=10.0)
    parser.add_argument("--p-min", type=float, default=1e-4)
    parser.add_argument("--p-max", type=float, default=1e-2)
    parser.add_argument("--num", type=int, default=9, help="p values, log-spaced")
    parser.add_argument("--shots", type=int, default=10_000, help="most shots per point")
    parser.add_argument("--max-errors", type=int, help="stop a point at this many logical errors")
    parser.add_argument("--out", type=Path, help="writes <out>.csv and <out>.png; default named from the settings "
                                                     "in results/physical_to_logical/")
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--plot-only", action="store_true", help="redraw <out>.png from <out>.csv")
    parser.add_argument("--sources", nargs="*", default=[], help="warn if any is newer than <out>.csv")
    args = parser.parse_args()
    if args.out is None:  # named like the Makefile's files
        args.out = Path("results/physical_to_logical") / (
            f"{args.noise}_eta{args.eta:g}_p{args.p_min:g}-{args.p_max:g}x{args.num}_shots{args.shots}"
            f"_{'-'.join(args.schemes)}_physical_to_logical")

    if args.plot_only:
        _plot(read_rows(Path(f"{args.out}.csv")), Path(f"{args.out}.png"))
        return

    for name, s in dem_stats(args.distances[0], args.noise, args.eta, 1e-3, args.schemes).items():
        print(f"{name:12} d={args.distances[0]}: " + ", ".join(f"{k} {v}" for k, v in s.items()))

    tasks = [dict(scheme=name, code=code, basis=b, distance=d, rounds=d, noise=args.noise, eta=args.eta, p=float(p), shots_max=args.shots,
                  # only when set, so CSVs written without it still resume
                  **({"max_errors": args.max_errors} if args.max_errors else {}))
             for code in args.codes for name in args.schemes for d in args.distances for b in "ZX"
             for p in np.geomspace(args.p_min, args.p_max, args.num)]
    rows = run_points(Path(f"{args.out}.csv"), tasks, _point, args.workers, sources=args.sources)
    _plot(rows, Path(f"{args.out}.png"))
    print(f"wrote {args.out}.csv, {args.out}.png")


@lru_cache(maxsize=None)
def _circuit(scheme: str, distance: int, code: str, basis: str) -> "stim.Circuit":
    """The noiseless memory circuit, d rounds: once per worker."""
    return floquet.memory_circuit(distance, distance, SCHEMES[scheme], basis, code)


def _point(task: dict) -> dict:
    """One (scheme, code, basis, distance, p) point of the sweep: its logical errors in up to shots_max shots, stopping at
    max_errors if set."""
    noisy = add_noise(_circuit(task["scheme"], task["distance"], task["code"], task["basis"]), NOISE_MODELS[task["noise"]](task["p"], task["eta"]))
    errors, ran = logical_errors(noisy, task["shots_max"], task.get("max_errors"))
    return dict(task, shots=ran, errors=errors, ler=errors / ran)


def _plot(rows: list[dict], png: Path) -> None:
    """LER (Z and X memory combined, `combined`) vs p: colour = distance, line style and marker = scheme; with both
    codes, CSS hollow and lighter. Legend blocks for only what is in `rows`. A point with no logical errors is a
    bar at its 95% upper bound with an arrow down: the LER is below there, not zero."""
    rows = combined(rows)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    fig.subplots_adjust(right=0.75)
    distances = sorted({r["distance"] for r in rows})
    schemes = [n for n in SCHEMES if any(r["scheme"] == n for r in rows)]  # fixed order, not the CSV's
    codes = [c for c in floquet.CODES if any(r["code"] == c for r in rows)]
    for k, d in enumerate(distances):
        color = DISTANCE_COLOUR[d]
        nudge = 1.06 ** (k - (len(distances) - 1) / 2)  # same shots give the same bound: spread the distances out
        for code in codes:
            hollow = len(codes) > 1 and code == "css"
            alpha = 0.6 if hollow else 1
            for scheme in schemes:
                ls, marker = STYLES[scheme]
                mine = sorted((r for r in rows if (r["distance"], r["code"], r["scheme"]) == (d, code, scheme)),
                              key=lambda r: r["p"])
                seen = [r for r in mine if r["errors"]]
                ax.loglog([r["p"] for r in seen], [r["ler"] for r in seen], marker=marker, ls=ls, color=color,
                          mfc="white" if hollow else color, alpha=alpha, lw=1.5, ms=4)
                for r in mine:
                    if not r["errors"]:
                        x = r["p"] * nudge
                        ax.plot(x, r["hi"], marker="_", color=color, ms=7, ls="none", alpha=alpha)
                        ax.plot(x, r["hi"] / 4, alpha=0)  # annotate doesn't widen the axis: keep the arrow inside it
                        ax.annotate("", xy=(x, r["hi"] / 4), xytext=(x, r["hi"]),
                                    arrowprops=dict(arrowstyle="->", color=color, lw=0.8, alpha=alpha))
    lo, hi = min(r["p"] for r in rows), max(r["p"] for r in rows)
    ax.axline((lo, lo), (hi, hi), color="0.5", ls=":", lw=1)  # through the data's own p: fixed points widen the axis
    memory = "Z or X memory" if any(r["basis"] == "ZX" for r in rows) else "Z memory"
    ax.set(xlabel="physical error rate p", ylabel=f"logical error rate, {memory} (d rounds)",
           title=f"{rows[0]['noise']} noise, η = {rows[0]['eta']:g}")
    ax.grid(which="major", color="#e4e4e0", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    legend_blocks(ax, schemes, distances, lines=True)
    second = ax.get_legend()
    extra = [Line2D([], [], marker="o", color="0.3", ls="none", label="X3Z3"),
             Line2D([], [], marker="o", mfc="white", color="0.3", alpha=0.6, ls="none", label="CSS")] if len(codes) > 1 else []
    ax.legend(handles=extra + [Line2D([], [], color="0.5", ls=":", label="LER = p"),
                               Line2D([], [], marker="$\\downarrow$", color="0.3", ls="none",
                                      label="0 errors: 95% upper bound")],
              fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 0.2), frameon=False)
    ax.add_artist(second)
    fig.savefig(png, dpi=200)
    plt.close(fig)

if __name__ == "__main__":
    main()
