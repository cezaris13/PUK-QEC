"""Memory experiment: logical error rate vs physical p, as CSV. `make plots` runs it per distance.
Each p is appended to the CSV as it finishes, and a rerun skips the p already there, so a cut-short run picks up
where it stopped (saved_rows / append_row, which every result script here uses).

    .venv/bin/python src/shared/memory.py --distance 3 --rounds 3 --noise spin --eta 10 --csv results/two_ancillas/spin_eta10_d3.csv

Run as a script it sweeps the two-ancilla scheme (src/two_ancillas/pairs.py). The noiseless circuit comes from
`floquet.memory_circuit`; `noise.add_noise` applies the model at every p.
"""
import argparse
import csv
import io
import sys
from pathlib import Path

import numpy as np
import pymatching
import stim
from tqdm import tqdm

from floquet import memory_circuit
from noise import NOISE_MODELS, add_noise, uniform_noise_model


def _value(text: str):
    for kind in (int, float):
        try:
            return kind(text)
        except ValueError:
            pass
    return text


def saved_rows(csv_path: Path) -> list[dict]:
    """The rows already in `csv_path`, numbers back as int or float; [] if it doesn't exist yet."""
    if not csv_path.exists():
        return []
    with open(csv_path, newline="") as f:
        return [{k: _value(v) for k, v in row.items()} for row in csv.DictReader(f)]


def point_key(row: dict, fields: tuple) -> tuple:
    """`row`'s values at `fields`, the same whether they come from a CSV (text) or from code (numbers)."""
    def norm(v):
        try:
            return repr(float(v))
        except (TypeError, ValueError):
            return str(v)
    return tuple(norm(row[k]) for k in fields)


def start_csv(csv_path: Path, fields, sources=()) -> None:
    """Create `csv_path` with its header, unless it already exists. Warns if any of `sources` (the code that
    simulates its points) is newer than it: its points may be stale, and only deleting it starts over."""
    # ponytail: warn, never delete: git sets checkout times, so a pull would make every result look stale
    newer = [str(f) for f in sources if csv_path.exists() and Path(f).stat().st_mtime > csv_path.stat().st_mtime]
    if newer:
        print(f"note: {', '.join(newer)} changed since {csv_path} was written; its points are reused. "
              f"Delete it to simulate them again.", flush=True)
    if not csv_path.exists():
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        csv_path.write_text(",".join(fields) + "\n")


def append_row(csv_path: Path, fields, row: dict) -> None:
    """Append `row` to `csv_path` (header from `start_csv`) and flush, so a crash loses at most this point."""
    line = io.StringIO()
    csv.DictWriter(line, fieldnames=fields, extrasaction="ignore", lineterminator="\n").writerow(row)
    # ponytail: one short append per row, so pool workers writing the same file don't interleave lines
    with open(csv_path, "a") as f:
        f.write(line.getvalue())


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
    """One row per p: the noise model at that p on `memory_circuit(distance, rounds, step, code=code)`, decoded;
    up to `shots` shots, fewer once `max_errors` logical errors are in (the row says how many ran)."""
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
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--sources", nargs="*", default=[], help="warn if any is newer than the CSV")
    args = parser.parse_args()

    # ponytail: the CLI runs the two-ancilla scheme only; decode.py runs every scheme
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "two_ancillas"))
    from pairs import round_circuit

    fields = ("distance", "rounds", "noise", "eta", "p", "shots_max", "shots", "errors", "ler")
    start_csv(args.csv, fields, args.sources)
    done = {point_key(r, ("p",)) for r in saved_rows(args.csv)}
    todo = [float(p) for p in np.geomspace(args.p_min, args.p_max, args.num) if point_key({"p": p}, ("p",)) not in done]
    print(f"{args.num - len(todo)} of {args.num} p already in {args.csv}", flush=True)
    circuit = memory_circuit(args.distance, args.rounds, round_circuit)
    for p in tqdm(todo, desc=f"d={args.distance}"):
        errors, ran = logical_errors(add_noise(circuit, NOISE_MODELS[args.noise](p, args.eta)), args.shots)
        append_row(args.csv, fields, dict(distance=args.distance, rounds=args.rounds, noise=args.noise, eta=args.eta,
                                          p=p, shots_max=args.shots, shots=ran, errors=errors, ler=errors / ran))
    print(f"wrote {args.csv}")


if __name__ == "__main__":
    main()
