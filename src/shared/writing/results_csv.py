"""The result CSVs every experiment writes: one row per simulated point, appended as soon as it is done, so a run
cut short (or rerun with more points) only computes what is missing (`run_points`)."""
import csv
import os
from multiprocessing import Pool
from pathlib import Path
from typing import Callable

from tqdm import tqdm


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


def write_rows(csv_path: Path, fields, rows: list[dict]) -> None:
    """Overwrite `csv_path` with `rows` under the header `fields`; keys not in `fields` are left out."""
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
