"""Memory experiment: logical error rate vs physical p, as CSV. `make plots` runs it per distance.

    .venv/bin/python src/memory.py --distance 3 --rounds 3 --noise spin --eta 10 --csv results/spin_eta10_d3.csv

The noiseless circuit comes from `floquet.memory_circuit`; `noise.add_noise` applies the model at every p.
"""
import argparse
import csv
from pathlib import Path

import numpy as np
import pymatching
import stim
from tqdm import tqdm

from floquet import memory_circuit
from noise import NOISE_MODELS, add_noise


def count_logical_errors(circuit: stim.Circuit, shots: int) -> int:
    """Shots, out of `shots`, where MWPM gets any logical observable wrong."""
    dem = circuit.detector_error_model(decompose_errors=True, approximate_disjoint_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(dem)
    detectors, observables = circuit.compile_detector_sampler().sample(shots, separate_observables=True)
    return int(np.any(matcher.decode_batch(detectors) != observables, axis=1).sum())


def sweep(distance: int, rounds: int, noise: str, eta: float, ps: list[float], shots: int) -> list[dict]:
    """One row per p: the noise model at that p on `memory_circuit(distance, rounds)`, decoded."""
    circuit = memory_circuit(distance, rounds)
    rows = []
    for p in tqdm(ps, desc=f"d={distance}"):
        errors = count_logical_errors(add_noise(circuit, NOISE_MODELS[noise](p, eta)), shots)
        rows.append(dict(distance=distance, rounds=rounds, noise=noise, eta=eta, p=p, shots=shots,
                         errors=errors, ler=errors / shots))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--distance", type=int, default=3)
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--noise", choices=NOISE_MODELS, default="spin")
    parser.add_argument("--eta", type=float, default=10.0)
    parser.add_argument("--p-min", type=float, default=1e-4)
    parser.add_argument("--p-max", type=float, default=1e-2)
    parser.add_argument("--num", type=int, default=9, help="p values, log-spaced")
    parser.add_argument("--shots", type=int, default=10_000)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()

    ps = list(np.geomspace(args.p_min, args.p_max, args.num))
    rows = sweep(args.distance, args.rounds, args.noise, args.eta, ps, args.shots)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with open(args.csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {args.csv}")


if __name__ == "__main__":
    main()
