"""Memory experiment: logical error rate vs physical p, as CSV. `make plots` runs it per distance.
Each p is appended to the CSV as it finishes, and a rerun skips the p already there, so a cut-short run picks up
where it stopped (writing.results_csv.run_points, which every result script here uses).

    .venv/bin/python src/shared/memory.py --distance 3 --rounds 3 --noise spin --eta 10 --csv results/two_ancillas/spin_eta10_d3.csv

Run as a script it sweeps the two-ancilla scheme (src/two_ancillas/pairs.py). The noiseless circuit comes from
`floquet.memory_circuit`; `noise.add_noise` applies the model at every p.
"""
import argparse
from functools import lru_cache
from pathlib import Path

import numpy as np
import pymatching
import stim
from tqdm import tqdm

from floquet import memory_circuit
from noise import NOISE_MODELS, add_noise, uniform_noise_model
from writing.results_csv import run_points


def uniform_matcher(circuit: stim.Circuit, p: float = 0.01) -> pymatching.Matching:
    """The decoder of IBM's QEC-with-spin-qubits, for any circuit: matching weights from one uniform error rate
    `p` (theirs is 0.01) instead of the real noise, so it knows nothing of biases or ratios between error
    sources. Theirs is a matching graph written by hand for the heavy-hex lattice; this one is Stim's from the
    noiseless `circuit` under noise.uniform_noise_model(p), built once and reused at every noise point."""
    dem = add_noise(circuit, uniform_noise_model(p)).detector_error_model(decompose_errors=True,
                                                                            approximate_disjoint_errors=True)
    return pymatching.Matching.from_detector_error_model(dem)


def logical_errors(circuit: stim.Circuit, shots: int, max_errors: int | None = None,
                   batch: int = 100_000, matcher: pymatching.Matching | None = None) -> tuple[int, int]:
    """(errors, shots run): shots where MWPM gets any logical observable wrong, sampled `batch` at a time so
    millions of shots fit in memory, stopping early once `max_errors` are seen. Decodes with `matcher`, or by
    default one built from `circuit`'s own error model, which knows the noise exactly."""
    if matcher is None:
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
          step, max_errors: int | None = None, code: str = "css") -> list[dict]:
    """One row per p, in memory (notebooks): the noise model at that p on `memory_circuit(distance, rounds, step,
    code=code)`, decoded; up to `shots` shots, fewer once `max_errors` logical errors are in."""
    circuit = memory_circuit(distance, rounds, step, code=code)
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
    parser.add_argument("--csv", type=Path, required=True, help="appended to point by point; a rerun resumes it")
    parser.add_argument("--sources", nargs="*", default=[], help="warn if any is newer than the CSV")
    args = parser.parse_args()

    tasks = [dict(distance=args.distance, rounds=args.rounds, noise=args.noise, eta=args.eta, p=float(p),
                  shots_max=args.shots) for p in np.geomspace(args.p_min, args.p_max, args.num)]
    run_points(args.csv, tasks, _pairs_point, desc=f"d={args.distance}", sources=args.sources)
    print(f"wrote {args.csv}")


@lru_cache(maxsize=None)
def _pairs_circuit(distance: int, rounds: int) -> stim.Circuit:
    from pairs import round_circuit  # ponytail: the CLI runs the two-ancilla scheme only; physical_to_logical.py runs every scheme
    return memory_circuit(distance, rounds, round_circuit)


def _pairs_point(task: dict) -> dict:
    """`main`'s point: the pairs memory experiment at the task's (distance, rounds, noise, eta, p)."""
    noisy = add_noise(_pairs_circuit(task["distance"], task["rounds"]), NOISE_MODELS[task["noise"]](task["p"], task["eta"]))
    errors, ran = logical_errors(noisy, task["shots_max"])
    return dict(task, shots=ran, errors=errors, ler=errors / ran)


if __name__ == "__main__":
    main()
