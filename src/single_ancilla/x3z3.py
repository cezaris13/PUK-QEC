"""The X3Z3 Floquet code against the CSS Floquet code, for every readout scheme, under biased noise: the figures
of Setiawan & McLauchlan (arXiv:2411.04974) for these circuits.

    .venv/bin/python src/single_ancilla/x3z3.py --noise spin --distances 3 5 --etas 0.5 1 10 100 1000
    .venv/bin/python src/single_ancilla/x3z3.py --plot-only --out <out>    # redraw from <out>.csv

Every (code, scheme, distance, eta, p) runs the Z and the X memory experiment (floquet.memory_circuit, d rounds,
noise.NOISE_MODELS[--noise] at (p, eta)), decoded by matching with the exact noise (memory.logical_errors).
Their logical error rates combine as the paper's Eq. (9), p_L = 1 - (1 - p_Z)(1 - p_X): either logical failing.
The 95% Clopper-Pearson interval of each basis combines the same way. The p grid is --p-min..--p-max plus --p-fixed.

    <out>_bias.png             (a) threshold vs eta: where the largest distance's p_L curve crosses the smallest's,
                               the lower of the Z and X memories', as paper Fig. 6; (b) p_L vs eta at --p-fixed for
                               the largest distance, as Fig. 7. Colour = scheme, solid X3Z3, dashed CSS.
    <out>_ler_eta<eta>.png     p_L vs p, one panel per scheme, both codes: colour = distance, solid X3Z3, dashed
                               CSS, shaded 95% interval.
    <out>_qubits_eta<eta>.png  p_L at --p-fixed vs total physical qubits on square-root / log axes, a line
                               log p_L = a + b sqrt(qubits) per (scheme, code): marker = scheme, filled X3Z3,
                               hollow CSS, colour = distance.

Each finished point is appended to <out>.csv at once, and a rerun with the same settings skips the points already
there, so a cut-short run picks up where it stopped. Points with no logical errors are left out of the plots.
"""
import argparse
import os
from functools import lru_cache
from pathlib import Path

import numpy as np

from decode import SCHEMES
from floquet import CODES, memory_circuit
from memory import logical_errors, run_points
from noise import NOISE_MODELS, add_noise
from x3z3_plots import plot

ROOT = Path(__file__).resolve().parents[2]


@lru_cache(maxsize=None)
def circuit(code: str, scheme: str, distance: int, basis: str):
    """The noiseless memory circuit, d rounds: once per worker."""
    return memory_circuit(distance, distance, SCHEMES[scheme], basis, code)


def point(task: dict) -> dict:
    c = circuit(task["code"], task["scheme"], task["distance"], task["basis"])
    errors, shots = logical_errors(add_noise(c, NOISE_MODELS[task["noise"]](task["p"], task["eta"])),
                                   task["shots_max"], task["max_errors"])
    return dict(task, qubits=c.num_qubits, shots=shots, errors=errors)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--codes", nargs="+", choices=CODES, default=list(CODES))
    parser.add_argument("--schemes", nargs="+", choices=SCHEMES, default=list(SCHEMES))
    parser.add_argument("--noise", choices=NOISE_MODELS, default="spin")
    parser.add_argument("--distances", type=int, nargs="+", default=[3, 5])
    parser.add_argument("--etas", type=float, nargs="+", default=[0.5, 1, 3, 10, 30, 100, 1000])
    parser.add_argument("--p-min", type=float, default=2e-4)
    parser.add_argument("--p-max", type=float, default=1.5e-2)
    parser.add_argument("--num", type=int, default=12, help="p values, log-spaced")
    parser.add_argument("--p-fixed", type=float, default=1e-3, help="p of the qubit-count figure and of (b)")
    parser.add_argument("--shots", type=int, default=20_000, help="most shots per point and basis")
    parser.add_argument("--max-errors", type=int, default=500, help="stop a point at this many logical errors")
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--out", type=Path, help="writes <out>.csv and <out>_*.png")
    parser.add_argument("--plot-only", action="store_true", help="redraw from <out>.csv")
    args = parser.parse_args()

    if args.out is None:
        args.out = ROOT / "results" / "single_ancilla" / "x3z3" / (
            f"{args.noise}_d{'-'.join(map(str, args.distances))}_p{args.p_min:g}-{args.p_max:g}x{args.num}"
            f"_shots{args.shots}")
    csv_path = args.out.with_suffix(".csv")
    if not args.plot_only:
        ps = sorted(set(np.geomspace(args.p_min, args.p_max, args.num).tolist()) | {args.p_fixed})
        tasks = [dict(code=c, scheme=s, distance=d, eta=eta, p=p, basis=b, noise=args.noise, shots_max=args.shots,
                      max_errors=args.max_errors)
                 for c in args.codes for s in args.schemes for d in args.distances for eta in args.etas
                 for p in ps for b in "ZX"]
        # slowest first (largest distance, then highest p), so the last workers aren't left with one big job
        run_points(csv_path, sorted(tasks, key=lambda t: (-t["distance"], -t["p"])), point, args.workers)
    for png in plot(csv_path, args.out, args.p_fixed):
        print(f"wrote {png}")


if __name__ == "__main__":
    main()
