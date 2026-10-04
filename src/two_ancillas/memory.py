"""Memory experiment: logical error rate vs physical p, as CSV. `make plots` runs it per distance.

    .venv/bin/python src/two_ancillas/memory.py --distance 3 --rounds 3 --noise spin --eta 10 --csv results/two_ancillas/spin_eta10_d3.csv

The noiseless circuit comes from `floquet.memory_circuit`; `noise.add_noise` applies the model at every p.
"""
import argparse
import csv
from pathlib import Path

import numpy as np
import pymatching
import stim
from tqdm import tqdm

from floquet import memory_circuit, round_circuit
from noise import NOISE_MODELS, add_noise


def logical_errors(circuit: stim.Circuit, shots: int, max_errors: int | None = None,
                   batch: int = 100_000) -> tuple[int, int]:
    """(errors, shots run): shots where MWPM gets any logical observable wrong, sampled `batch` at a time so
    millions of shots fit in memory, stopping early once `max_errors` are seen."""
    dem = circuit.detector_error_model(decompose_errors=True, approximate_disjoint_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(dem)
    sampler = circuit.compile_detector_sampler()
    errors = done = 0
    while done < shots and (max_errors is None or errors < max_errors):
        n = min(batch, shots - done)
        detectors, observables = sampler.sample(n, separate_observables=True)
        errors += int(np.any(matcher.decode_batch(detectors) != observables, axis=1).sum())
        done += n
    return errors, done


def count_logical_errors(circuit: stim.Circuit, shots: int) -> int:
    """Shots, out of `shots`, where MWPM gets any logical observable wrong."""
    return logical_errors(circuit, shots)[0]


def sweep(distance: int, rounds: int, noise: str, eta: float, ps: list[float], shots: int,
          step=round_circuit, max_errors: int | None = None) -> list[dict]:
    """One row per p: the noise model at that p on `memory_circuit(distance, rounds, step)`, decoded; up to
    `shots` shots, fewer once `max_errors` logical errors are in (the row says how many ran)."""
    circuit = memory_circuit(distance, rounds, step)
    rows = []
    for p in tqdm(ps, desc=f"d={distance}"):
        errors, ran = logical_errors(add_noise(circuit, NOISE_MODELS[noise](p, eta)), shots, max_errors)
        rows.append(dict(distance=distance, rounds=rounds, noise=noise, eta=eta, p=p, shots=ran,
                         errors=errors, ler=errors / ran))
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
