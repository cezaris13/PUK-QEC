"""Memory experiment for every readout scheme, decoded the same way: logical error rate vs p.

    .venv/bin/python src/single_ancilla/decode.py --distances 3 5 --noise spin --eta 10 --shots 10000
    .venv/bin/python src/single_ancilla/decode.py --schemes method_a method_b method_c

"pairs" is src/two_ancillas/pairs.py's syndrome + reference ancilla pair on every edge; the rest have one ancilla per
edge: "method_a", "method_b" and "method_c" are hexes.pdf's methods A, B and C (src/single_ancilla/<name>/). All run through
floquet.memory_circuit, so they share every detector and the observable, and only their step circuits (and so
their noise) differ. Decoding is memory.count_logical_errors: Stim's detector error model, decomposed into
a graph, matched by PyMatching (docs/decoding.pdf). Writes one CSV row per (scheme, distance, p) and a plot.
Each row is also appended to <out>.partial.jsonl as it finishes, so rerunning the same command after a crash
only runs the points still missing.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.stats import beta
from tqdm import tqdm

SRC = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(SRC / "shared"), str(SRC / "two_ancillas")]
import floquet
import pairs
from memory import logical_errors
from noise import NOISE_MODELS, add_noise

for name in ("method_a", "method_b", "method_c"):
    sys.path.insert(0, str(Path(__file__).resolve().parent / name))
import method_a
import method_b
import method_c

SCHEMES = {"pairs": pairs.round_circuit, "method_a": method_a.round_circuit, "method_b": method_b.round_circuit,
           "method_c": method_c.round_circuit}
# line style and marker per scheme; colour is the distance
STYLES = {"pairs": ("-", "o"), "method_a": ("--", "D"), "method_b": ("-.", "^"), "method_c": (":", "v")}
# colour per distance: the N2E3N2 paper's Figure 4 green / blue / red for the odd distances used here
DISTANCE_COLOUR = {3: "#2ca02c", 5: "#1f77b4", 7: "#d62728", 9: "#9467bd",
                   4: "#8fce8f", 6: "#8fb4e8", 8: "#f0928f", 10: "#c5b0d5"}


def clopper_pearson(errors: int, shots: int) -> tuple[float, float]:
    """95% Clopper-Pearson interval of errors / shots."""
    low = beta.ppf(0.025, errors, shots - errors + 1) if errors else 0.0
    high = beta.ppf(0.975, errors + 1, shots - errors) if errors < shots else 1.0
    return low, high


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


def read_rows(csv_path: Path) -> list[dict]:
    """A CSV written here, numbers back as int or float."""
    def number(text):
        for kind in (int, float):
            try:
                return kind(text)
            except ValueError:
                pass
        return text
    with open(csv_path, newline="") as f:
        return [{k: number(v) for k, v in row.items()} for row in csv.DictReader(f)]


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


def plot(rows: list[dict], png: Path) -> None:
    """LER vs p: colour = distance, line style = scheme. Zero-error points are left out (log axis)."""
    fig, ax = plt.subplots(figsize=(6, 4.5))
    distances = sorted({r["distance"] for r in rows})
    for d in distances:
        color = DISTANCE_COLOUR[d]
        for scheme, (ls, marker) in STYLES.items():
            pts = [r for r in rows if r["distance"] == d and r["scheme"] == scheme and r["errors"]]
            ax.loglog([r["p"] for r in pts], [r["ler"] for r in pts], marker=marker, ls=ls, color=color, lw=1.5, ms=4,
                      label=f"d = {d}, {scheme}")
    ax.axline((1e-3, 1e-3), (1e-2, 1e-2), color="0.5", ls=":", lw=1, label="LER = p")
    ax.set(xlabel="physical error rate p", ylabel="logical error rate",
           title=f"{rows[0]['noise']} noise, η = {rows[0]['eta']:g}")
    ax.grid(which="major", color="#e4e4e0", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--distances", type=int, nargs="+", default=[3, 5])
    parser.add_argument("--schemes", nargs="+", choices=SCHEMES, default=list(SCHEMES))
    parser.add_argument("--noise", choices=NOISE_MODELS, default="spin")
    parser.add_argument("--eta", type=float, default=10.0)
    parser.add_argument("--p-min", type=float, default=1e-4)
    parser.add_argument("--p-max", type=float, default=1e-2)
    parser.add_argument("--num", type=int, default=9, help="p values, log-spaced")
    parser.add_argument("--shots", type=int, default=10_000)
    parser.add_argument("--out", type=Path, default=Path("results/single_ancilla/decode"), help="writes <out>.csv and <out>.png")
    parser.add_argument("--plot-only", action="store_true", help="redraw <out>.png from <out>.csv")
    args = parser.parse_args()

    if args.plot_only:
        plot(read_rows(args.out.with_suffix(".csv")), args.out.with_suffix(".png"))
        return

    for name, s in dem_stats(args.distances[0], args.noise, args.eta, 1e-3, args.schemes).items():
        print(f"{name:12} d={args.distances[0]}: " + ", ".join(f"{k} {v}" for k, v in s.items()))

    ps = list(np.geomspace(args.p_min, args.p_max, args.num))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    partial = args.out.with_suffix(".partial.jsonl")
    done = {}
    if partial.exists():
        for line in partial.read_text().splitlines():
            row = json.loads(line)
            done[(row["scheme"], row["distance"], row["p"])] = row
    points = [(name, d, p) for name in args.schemes for d in args.distances for p in ps]
    with open(partial, "a") as f, tqdm(total=len(points), initial=sum(k in done for k in points)) as bar:
        for name, d, p in points:
            if (name, d, p) in done:
                continue
            bar.set_description(f"{name} d={d}")
            circuit = floquet.memory_circuit(d, d, SCHEMES[name])
            errors, ran = logical_errors(add_noise(circuit, NOISE_MODELS[args.noise](p, args.eta)), args.shots)
            done[(name, d, p)] = row = dict(scheme=name, distance=d, rounds=d, noise=args.noise, eta=args.eta,
                                            p=float(p), shots=ran, errors=errors, ler=errors / ran)
            f.write(json.dumps(row) + "\n")
            f.flush()
            bar.update()
    rows = [done[k] for k in points]
    with open(args.out.with_suffix(".csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    plot(rows, args.out.with_suffix(".png"))
    partial.unlink()
    print(f"wrote {args.out}.csv, {args.out}.png")


if __name__ == "__main__":
    main()
