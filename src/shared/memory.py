"""Memory experiment: logical error rate vs physical p, as CSV. `make plots` runs it per distance.
Each p is appended to the CSV as it finishes, and a rerun skips the p already there, so a cut-short run picks up
where it stopped (run_points, which every result script here uses).

    .venv/bin/python src/shared/memory.py --distance 3 --rounds 3 --noise spin --eta 10 --csv results/two_ancillas/spin_eta10_d3.csv

Run as a script it sweeps the two-ancilla scheme (src/two_ancillas/pairs.py). The noiseless circuit comes from
`floquet.memory_circuit`; `noise.add_noise` applies the model at every p.
"""
import argparse
import csv
import os
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path
from typing import Callable

import numpy as np
import pymatching
import stim
from tqdm import tqdm

from floquet import memory_circuit
from noise import NOISE_MODELS, add_noise, uniform_noise_model


def read_rows(csv_path: Path) -> list[dict]:
    """A CSV written here, numbers back as int or float; [] if there is none yet."""
    def number(text):
        for kind in (int, float):
            try:
                return kind(text)
            except ValueError:
                pass
        return text
    if not Path(csv_path).exists():
        return []
    with open(csv_path, newline="") as f:
        return [{k: number(v) for k, v in row.items()} for row in csv.DictReader(f)]


def point_key(row: dict, fields) -> tuple:
    """A point's identity: its input `fields`, as text, so a row read back from the CSV matches its task."""
    return tuple(str(row[f]) for f in fields)


def append_row(csv_path: Path, row: dict) -> None:
    """Append one row to `csv_path`, with the header first if the file is new; flushed at once, so a crash loses
    nothing that finished. Every row of a file must have the same columns, written in the file's order."""
    new = not Path(csv_path).exists() or Path(csv_path).stat().st_size == 0
    fields = list(row) if new else next(csv.reader(open(csv_path, newline="")))  # the file's own column order
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        if new:
            writer.writeheader()
        writer.writerow(row)


def warn_stale(csv_path: Path, sources=()) -> None:
    """Warn if any of `sources` (the code that simulates `csv_path`'s points) is newer than it: its points may be
    stale, and only deleting it starts over."""
    # ponytail: warn, never delete: git sets checkout times, so a pull would make every result look stale
    csv_path = Path(csv_path)
    newer = [str(f) for f in sources if csv_path.exists() and Path(f).stat().st_mtime > csv_path.stat().st_mtime]
    if newer:
        print(f"note: {', '.join(newer)} changed since {csv_path} was written; its points are reused. "
              f"Delete it to simulate them again.", flush=True)


def run_points(csv_path: Path, tasks: list[dict], compute: Callable[[dict], dict], workers: int = os.cpu_count(),
               desc: str | None = None, sources=()) -> list[dict]:
    """Every scheme's data collection: the row of each task (a dict of the point's inputs), in task order.
    Rows already in `csv_path` are read back; the rest are `compute(task)`d (a top-level function returning
    the task plus its results) in a Pool of `workers` and appended as each finishes. So a run cut short, or
    rerun with more points, only computes what is missing; `warn_stale` if any of `sources` changed since."""
    if not tasks:
        return []
    warn_stale(csv_path, sources)
    fields = list(tasks[0])
    Path(csv_path).parent.mkdir(parents=True, exist_ok=True)
    done = {point_key(r, fields): r for r in read_rows(csv_path)}
    todo = [t for t in tasks if point_key(t, fields) not in done]
    with Pool(workers) as pool, tqdm(total=len(tasks), initial=len(tasks) - len(todo), desc=desc) as bar:
        for row in pool.imap_unordered(compute, todo):
            append_row(csv_path, row)
            done[point_key(row, fields)] = row
            bar.update()
    return [done[point_key(t, fields)] for t in tasks]


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
