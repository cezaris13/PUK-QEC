import argparse
import csv
from pathlib import Path

import cairosvg
import matplotlib.pyplot as plt

from floquet import memory_circuit

# Categorical slots in fixed order, one per distance (dataviz reference palette, light mode).
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]


def timeslices_png(distance: int, rounds: int, png: Path) -> None:
    svg = str(memory_circuit(distance, rounds).diagram("timeslice-svg"))
    cairosvg.svg2png(bytestring=svg.encode(), write_to=str(png), scale=2, background_color="white")


def draw_ler(ax: plt.Axes, runs: list[list[dict]]) -> None:
    """One log-log line of LER vs p per run (a run = one distance's rows); the crossing is the threshold."""
    runs = sorted(runs, key=lambda rows: int(rows[0]["distance"]))
    if len(runs) > len(COLORS):
        raise ValueError(f"{len(runs)} distances, only {len(COLORS)} distinguishable colours; plot fewer")
    for color, rows in zip(COLORS, runs):
        label = f"d = {rows[0]['distance']}"
        rows = [r for r in rows if int(r["errors"])]  # zero errors has no place on a log axis
        ax.loglog([float(r["p"]) for r in rows], [float(r["ler"]) for r in rows], "o-",
                  color=color, lw=1.5, ms=5, label=label)
    ax.set(xlabel="physical error rate p", ylabel="logical error rate",
           title=f"{runs[0][0]['noise']} noise, η = {float(runs[0][0]['eta']):g}")
    ax.grid(which="major", color="#e4e4e0", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)


def ler_png(csv_paths: list[Path], png: Path) -> None:
    runs = []
    for path in csv_paths:
        with open(path, newline="") as f:
            runs.append(list(csv.DictReader(f)))
    fig, ax = plt.subplots(figsize=(6, 4.5))
    draw_ler(ax, runs)
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="figure", required=True)
    ts = sub.add_parser("timeslices")
    ts.add_argument("--distance", type=int, default=3)
    ts.add_argument("--rounds", type=int, default=3)
    ts.add_argument("--png", type=Path, required=True)
    ler = sub.add_parser("ler")
    ler.add_argument("csvs", type=Path, nargs="+")
    ler.add_argument("--png", type=Path, required=True)
    args = parser.parse_args()

    args.png.parent.mkdir(parents=True, exist_ok=True)
    if args.figure == "timeslices":
        timeslices_png(args.distance, args.rounds, args.png)
    else:
        ler_png(args.csvs, args.png)
    print(f"wrote {args.png}")


if __name__ == "__main__":
    main()
