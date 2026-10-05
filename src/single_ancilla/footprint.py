"""How many physical qubits each readout scheme needs at one physical error rate: logical error rate vs qubit
count, one line per scheme, extrapolated to a target logical error rate.

    .venv/bin/python src/single_ancilla/footprint.py --p 1e-3 --distances 3 5 --shots 10000000 --target 1e-6 \
        --schemes method_a method_b method_c --out results/single_ancilla/footprint

Each point is decode.py's memory experiment (d rounds at distance d) at p, decoded by matching. The LER falls
exponentially in d, so a straight line through log(LER) vs d (points with at least one logical error) gives
the d where it reaches --target; the qubit count there is the scheme's qubits per d^2 times d^2. Points with no
errors are left out: at low p the largest distances need a lot of shots to show any.
"""
import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from decode import SCHEMES, STYLES
from floquet import memory_circuit
from memory import sweep
from noise import NOISE_MODELS
from drawing import COLORS


def measure(schemes: list[str], distances: list[int], noise: str, eta: float, p: float, shots: int,
            max_errors: int) -> list[dict]:
    """One row per (scheme, distance): its qubit count and the LER at p, from up to `shots` shots, stopping at
    `max_errors` logical errors."""
    return [dict(scheme=name, qubits=memory_circuit(d, 1, SCHEMES[name]).num_qubits, **row)
            for name in schemes for d in distances
            for row in sweep(d, d, noise, eta, [p], shots, SCHEMES[name], max_errors)]


def qubits_needed(rows: list[dict], target: float) -> tuple[float, float, float] | None:
    """(slope, intercept, qubits) of the fit log(LER) = intercept + slope * d for one scheme's rows, with the
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


def plot(rows: list[dict], target: float, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for color, name in zip(COLORS, dict.fromkeys(r["scheme"] for r in rows)):
        mine = [r for r in rows if r["scheme"] == name]
        seen = [r for r in mine if r["errors"]]
        ls, marker = STYLES[name]
        fit = qubits_needed(mine, target)
        label = name + (f": {fit[2]:,.0f} qubits" if fit else ": no fit")
        ax.loglog([r["qubits"] for r in seen], [r["ler"] for r in seen], marker=marker, ls="none", color=color,
                  ms=6, label=label)
        if fit:
            slope, intercept, needed = fit
            per_d2 = mine[0]["qubits"] / mine[0]["distance"] ** 2
            ds = np.linspace(min(r["distance"] for r in mine), math.sqrt(needed / per_d2), 50)
            ax.loglog(per_d2 * ds ** 2, np.exp(intercept + slope * ds), ls=ls, color=color, lw=1.2)
            ax.plot(needed, target, marker="*", color=color, ms=10)
    ax.axhline(target, color="0.5", ls=":", lw=1, label=f"target LER {target:g}")
    ax.set(xlabel="physical qubits (data + ancillas)", ylabel="logical error rate (d rounds)",
           title=f"{rows[0]['noise']} noise, η = {rows[0]['eta']:g}, p = {rows[0]['p']:g}")
    ax.grid(which="major", color="#e4e4e0", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--schemes", nargs="+", choices=SCHEMES, default=list(SCHEMES))
    parser.add_argument("--distances", type=int, nargs="+", default=[3, 5, 7])
    parser.add_argument("--noise", choices=NOISE_MODELS, default="spin")
    parser.add_argument("--eta", type=float, default=10.0)
    parser.add_argument("--p", type=float, default=1e-3)
    parser.add_argument("--shots", type=int, default=10_000_000, help="most shots per point")
    parser.add_argument("--max-errors", type=int, default=100, help="stop a point at this many logical errors")
    parser.add_argument("--target", type=float, default=1e-6, help="logical error rate to extrapolate to")
    parser.add_argument("--out", type=Path, default=Path("results/single_ancilla/footprint"), help="writes <out>.csv and <out>.png")
    args = parser.parse_args()

    rows = measure(args.schemes, args.distances, args.noise, args.eta, args.p, args.shots, args.max_errors)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out.with_suffix(".csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    for name in args.schemes:
        fit = qubits_needed([r for r in rows if r["scheme"] == name], args.target)
        print(f"{name:12} " + (f"{fit[2]:,.0f} qubits for LER {args.target:g}" if fit else
                                "not enough distances with errors to fit (more shots or larger p)"))
    plot(rows, args.target, args.out.with_suffix(".png"))
    print(f"wrote {args.out}.csv, {args.out}.png")


if __name__ == "__main__":
    main()
