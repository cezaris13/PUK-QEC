"""Threshold surface and threshold-vs-bias plots for any readout scheme (pairs, method_a/b/c), ported from IBM's
QEC-with-spin-qubits (https://github.com/IBM/QEC-with-spin-qubits, Hetenyi & Wootton, "Tailoring quantum error
correction to spin qubits", arXiv:2306.17786): their figures (b) and (c), for these circuits.

    .venv/bin/python src/experiments/thresholds.py --scheme method_b --distances 3 5 --nphi 6
    .venv/bin/python src/experiments/thresholds.py --scheme method_b --plot-only   # redraw from the JSON
    .venv/bin/python src/experiments/thresholds.py --scheme method_b --code x3z3   # the X3Z3 Floquet code

Noise (noise.spin_qubit_noise_model) has three independent sources, gates p_G, idling p_T and readout p_R. A
direction (theta, phi) and a size p set them, as in IBM's plot_utils.LogFail_of_d_p:

    p_G = p cos(theta) cos(phi)   split as 1-qubit 2 p_G / (1 + eta_G), 2-qubit 2 p_G eta_G / (1 + eta_G)
    p_T = p cos(theta) sin(phi)   split as relaxation p_T / (1 + eta_T), dephasing p_T eta_T / (1 + eta_T)
    p_R = p sin(theta)

so eta_G is the 2-qubit to 1-qubit gate error ratio and eta_T the dephasing to relaxation idle ratio. Along
each direction the threshold p_th is where the logical error rates of the distances cross
(`_threshold_from_log_fail`, IBM's Threshold_from_LogFail).

(b) Over a grid of directions spanning the octant, the threshold points (p_G, p_T, p_R) * p_th / p form a
surface, coloured by its length p_th. (c) Along pure gate noise the threshold vs eta_G, along pure idling vs
eta_T; their ranges are the arrows on the axes of (b).

Decoding is IBM's idea (memory.uniform_matcher): weights from one uniform error rate, blind to biases. Where
IBM hand-wrote a graph for the heavy-hex lattice, Stim builds it here for whichever circuit, so it works for
every scheme. `--decoder dem` decodes with the exact noise instead. Every direction runs the Z and the X memory
experiment, and its threshold is the lower of the two, as IBM's: a Z memory never sees dephasing, an X memory
never relaxation-like X errors. Where the curves don't cross inside the scanned window, the window moves 2.5x
towards the threshold, up to --retries times. Without --pg/--pt/--pr the axis thresholds that aim the grid come
from a coarse scan first.

Writes results/two_ancillas/thresholds/ (pairs) or results/single_ancilla/<scheme>/thresholds/: <name>.csv with
every simulated point, one row per (direction, basis, distance, p) appended as soon as it is done, and
<name>.png, and <name>_thresholds.csv, one row per direction with its threshold (what the viewer behind
`make serve` reads; derived, like the PNG). Each scan reads its points from the CSV before simulating any, so rerunning the same command after a crash
only runs the points still missing, and a rerun with other settings (more NBIAS, say) reuses every point they
share; --plot-only replays the scans from the CSV alone, or draws <name>_thresholds.csv if it has no CSV.
--sources warns if that code changed since <name>.csv was written.
"""
# Parts of this file are adapted from IBM's QEC-with-spin-qubits, (C) Copyright IBM 2023, licensed under the
# Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0): LogFail_of_d_p,
# Threshold_from_LogFail and plot_3d_threshold of plot_utils.py, the direction grid of generate_3d_plot.py, the
# bias lists of generate_Gbias_plot.py / generate_Tbias_plot.py, and Arrow3D and the bias figure of
# Threshold_surfaces_demo.ipynb. Modified: rewritten for this repo's circuits, noise and decoder.
import argparse
import os
from multiprocessing import Lock
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, colors
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from scipy.interpolate import griddata
from tqdm import tqdm

from physical_to_logical import SCHEMES
from floquet import CODES, memory_circuit
from memory import logical_errors, uniform_matcher
from noise import add_noise, spin_qubit_noise_model
from writing.results_csv import append_row, point_key, read_rows, warn_stale, write_rows

ROOT = Path(__file__).resolve().parents[2]
# <name>_thresholds.csv: one row per direction; group is axis (the three axis thresholds aiming the grid, in
# the order G, T, R), surface, g_bias or t_bias
SUMMARY_FIELDS = ("group", "theta", "phi", "eta_g", "eta_t", "p_th", "p_th_error", "p_th_Z", "p_th_X", "scheme",
                  "code", "decoder", "distances")
AXES = ((0, 0), (0, np.pi / 2), (np.pi / 2, 0))  # (theta, phi) of pure gate, idling and readout noise
# One simulated point: what it is (POINT) and its result (shots, errors). A worker's view of <out>.csv, set by `_init`.
POINT = ("scheme", "code", "decoder", "basis", "distance", "p", "theta", "phi", "eta_g", "eta_t", "shots_max", "batch",
         "max_fail")
_csv: Path | None = None
_lock = None
_done: dict = {}


class Arrow3D(FancyArrowPatch):
    """An arrow between two 3D points (Threshold_surfaces_demo.ipynb)."""

    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs, ys, zs = proj3d.proj_transform(*self._verts3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scheme", choices=SCHEMES, default="pairs")
    parser.add_argument("--code", choices=CODES, default="css", help="css, or x3z3: the X3Z3 Floquet code")
    parser.add_argument("--distances", type=int, nargs="+", default=[3, 5])
    parser.add_argument("--decoder", choices=["uniform", "dem"], default="uniform",
                        help="uniform: IBM's bias-blind weights (default); dem: weights from the exact noise")
    parser.add_argument("--eta-g", type=float, default=1.0, help="2-qubit / 1-qubit gate error, for (b)")
    parser.add_argument("--eta-t", type=float, default=20.0, help="dephasing / relaxation idling, for (b)")
    parser.add_argument("--pg", type=float, help="threshold along pure gate noise, to aim the grid")
    parser.add_argument("--pt", type=float, help="threshold along pure idling")
    parser.add_argument("--pr", type=float, help="threshold along pure readout")
    parser.add_argument("--nphi", type=int, default=6, help="grid of (b): nphi (nphi + 1) / 2 directions")
    parser.add_argument("--nbias", type=int, default=9, help="eta values from 0.01 to 100 for (c)")
    parser.add_argument("--num-p", type=int, default=10, help="error rates per direction")
    parser.add_argument("--bias-num-p", type=int, help="error rates per direction of (c), default --num-p")
    parser.add_argument("--delpth", type=float, default=0.5, help="(b) scans p_guess * (1 +- delpth)")
    parser.add_argument("--bias-delpth", type=float, default=1.0, help="(c) scans p_axis * (1 +- this)")
    parser.add_argument("--shots", type=int, default=20_000, help="most shots per point")
    parser.add_argument("--batch", type=int, default=2_000)
    parser.add_argument("--max-fail", type=int, default=2_000, help="stop a point at this many logical errors")
    parser.add_argument("--max-fail-rate", type=float, default=0.45, help="stop raising p past this LER")
    parser.add_argument("--retries", type=int, default=3, help="window moves when the curves don't cross")
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--out", type=Path, help="writes <out>.csv (every point), <out>_thresholds.csv and <out>.png")
    parser.add_argument("--plot-only", action="store_true",
                        help="redraw <out>.png from <out>.csv, simulating nothing (or from <out>_thresholds.csv)")
    parser.add_argument("--bias-only", action="store_true",
                        help="run only (c) and replace it in <out>_thresholds.csv, keeping its surface and axes")
    parser.add_argument("--sources", nargs="*", default=[], help="warn if any is newer than <out>.csv")
    args = parser.parse_args()

    if args.out is None:
        folder = (ROOT / "results" / "two_ancillas" if args.scheme == "pairs"
                  else ROOT / "results" / "single_ancilla" / args.scheme) / "thresholds"
        args.out = folder / (("" if args.code == "css" else f"{args.code}_") + f"{args.scheme}_d{'-'.join(map(str, args.distances))}_{args.decoder}"
                             f"_etaG{args.eta_g:g}_etaT{args.eta_t:g}_shots{args.shots}_nphi{args.nphi}")
    # not with_suffix: a name like ..._etaG0.5_... has a dot of its own
    csv_path, summary, png = (args.out.parent / (args.out.name + end) for end in (".csv", "_thresholds.csv", ".png"))
    saved = _data_from_summary(read_rows(summary)) if summary.exists() else None
    if args.plot_only and saved and not csv_path.exists():
        data = saved  # a run from before the CSV: nothing to replay
    elif args.plot_only:
        if saved:  # aim at the saved axes: a run converted from before the CSV never kept its coarse scan
            args.pg, args.pt, args.pr = (a if a is not None else b for a, b in zip((args.pg, args.pt, args.pr),
                                                                                  saved["axes"]))
        try:
            with Pool(args.workers, initializer=_init, initargs=(csv_path, Lock())) as pool:
                data = _run(args, pool)
        except KeyError as missing:  # a point the replay needs isn't in the CSV
            if not saved:
                raise
            print(f"replay incomplete ({missing}); drawing {summary.name} as saved")
            data = saved
    else:
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        warn_stale(csv_path, args.sources)
        old = saved if args.bias_only else None
        if old:  # aim (c) at the axis thresholds the saved run found, so no coarse scan
            args.pg, args.pt, args.pr = (a if a is not None else b for a, b in zip((args.pg, args.pt, args.pr),
                                                                                  old["axes"]))
        with Pool(args.workers, initializer=_init, initargs=(csv_path, Lock())) as pool:
            data = _run(args, pool, args.bias_only)
        if args.bias_only:  # keep the saved surface, replace (c)
            new = data
            data = old or dict(new, surface=[])
            data.update(g_bias=new["g_bias"], t_bias=new["t_bias"])
        _write_summary(data, summary)
    for key, name in (("g_bias", "eta_g"), ("t_bias", "eta_t")):
        print(f"{key}: " + ", ".join(f"{r[name]:g}: {100 * r['p_th']:.2f}%" for r in data[key]))
    _plot(data, png)
    print(f"wrote {csv_path}, {summary}, {png}")


def _noise_at(p: float, theta: float, phi: float, eta_g: float, eta_t: float):
    """The noise model at size p in direction (theta, phi), as IBM's LogFail_of_d_p sets its error rates."""
    # max(0, .): on the octant's walls cos/sin of pi/2 can round to -1e-17, which Stim refuses as a probability
    p_g, p_t, p_r = (max(0.0, x) for x in (p * np.cos(theta) * np.cos(phi), p * np.cos(theta) * np.sin(phi),
                                            p * np.sin(theta)))
    return spin_qubit_noise_model(p_g1=2 * p_g / (1 + eta_g), p_g2=2 * p_g * eta_g / (1 + eta_g),
                                  p_t1=p_t / (1 + eta_t), p_t2=p_t * eta_t / (1 + eta_t), p_r=p_r)


@lru_cache(maxsize=None)
def _circuit_and_matcher(scheme: str, distance: int, decoder: str, basis: str, code: str = "css"):
    """The noiseless memory circuit (d rounds) and, for the uniform decoder, its fixed matcher: once per worker."""
    circuit = memory_circuit(distance, distance, SCHEMES[scheme], basis, code)
    return circuit, uniform_matcher(circuit) if decoder == "uniform" else None


def _init(csv_path: Path, lock) -> None:
    """Pool initializer: every worker reads the points already in `csv_path`, and appends new ones under `lock`."""
    global _csv, _lock, _done
    _csv, _lock = csv_path, lock
    _done = {point_key(r, POINT): r for r in read_rows(csv_path)}


def _log_fail(task: dict) -> list:
    """[[d, [[p, LER], ...]], ...] for one direction: IBM's LogFail_of_d_p. Each distance runs its error rates
    in order and stops once the LER passes max_fail_rate; each point runs batches of `batch` shots until
    `shots` or `max_fail` logical errors, or is read back from the CSV if it is there already (with
    task["cache_only"], it has to be)."""
    out = []
    for d in task["distances"]:
        rows = []
        for p in task["errors"]:
            if rows and rows[-1][1] >= task["max_fail_rate"]:
                break
            point = dict(scheme=task["scheme"], code=task["code"], decoder=task["decoder"], basis=task["basis"],
                         distance=d, p=float(p), theta=float(task["theta"]), phi=float(task["phi"]),
                         eta_g=float(task["eta_g"]), eta_t=float(task["eta_t"]), shots_max=task["shots"],
                         batch=task["batch"], max_fail=task["max_fail"])
            row = _done.get(point_key(point, POINT))
            if row is None:
                if task.get("cache_only"):
                    raise KeyError(f"{point} is not in {_csv}: run without --plot-only")
                circuit, matcher = _circuit_and_matcher(task["scheme"], d, task["decoder"], task["basis"], task["code"])
                noisy = add_noise(circuit, _noise_at(p, task["theta"], task["phi"], task["eta_g"], task["eta_t"]))
                errors, ran = logical_errors(noisy, task["shots"], task["max_fail"], task["batch"], matcher)
                row = dict(point, shots=ran, errors=errors)
                with _lock:
                    append_row(_csv, row)
            rows.append([float(p), row["errors"] / row["shots"]])
        out.append([d, rows])
    return out


# A p only counts towards the crossing when every distance saw at least this many logical errors there
# (LER >= MIN_ERRORS / shots): below that, counting noise alone puts d = 5 above d = 3 now and then, and IBM's
# rule took the lowest such p as the threshold, which gave the dips in the surfaces and bias plots.
MIN_ERRORS = 10


def _threshold_from_log_fail(log_fail_d_p: list, min_ler: float = 0.0, cutoff: float = 0.495) -> tuple[float, float]:
    """(p_th, error), (0, 0) if the curves never cross: IBM's Threshold_from_LogFail. Below threshold the
    LER falls with d at every p, above it rises; p_th sits between the highest p of the first kind and the
    lowest of the second, nearer the one whose curves are closer together. Only p where every distance's LER
    is at least `min_ler` (and above 0) count."""
    n = min(next((i for i, (_, ler) in enumerate(rows) if ler >= cutoff), len(rows) - 1)
            for _, rows in log_fail_d_p)
    errors = [log_fail_d_p[-1][1][i][0] for i in range(n)]
    below, above, spread = [], [], []
    for i, p in enumerate(errors):
        lers = [rows[i][1] for _, rows in log_fail_d_p]
        spread.append(max(lers) - min(lers))
        if min(lers) > 0 and min(lers) >= min_ler:  # too few failures to tell the distances apart
            if lers[::-1] == sorted(lers):
                below.append(p)
            if lers == sorted(lers):
                above.append(p)
    if not below or not above:
        return 0.0, 0.0
    lo, hi = max(below), min(above)
    s_lo, s_hi = spread[errors.index(lo)], spread[errors.index(hi)]
    p_th = (s_hi * lo + s_lo * hi) / (s_lo + s_hi) if s_lo + s_hi else (lo + hi) / 2
    return p_th, max(p_th - lo, hi - p_th)


def _scan(task: dict) -> list:
    """`_log_fail`, and if the distances don't cross, again with the window moved 2.5x towards the threshold:
    down if the smallest distance already fails often at the lowest p, up otherwise; up to task["retries"]
    times. Error rates above 0.3 are dropped (some channels would go past probability 1)."""
    errors = task["errors"]
    for _ in range(task["retries"] + 1):
        errors = [p for p in errors if p <= 0.3]
        lf = _log_fail(dict(task, errors=errors))
        if not errors or _threshold_from_log_fail(lf, MIN_ERRORS / task["shots"])[0] > 0:
            break
        factor = 1 / 2.5 if lf[0][1][0][1] >= 0.1 else 2.5
        errors = [p * factor for p in errors]
    return lf


def _combine(lf_z: list, lf_x: list, min_ler: float) -> dict:
    """The lower of the Z and X memories' thresholds (one alone if the other found none), as IBM's."""
    (z, z_err), (x, x_err) = _threshold_from_log_fail(lf_z, min_ler), _threshold_from_log_fail(lf_x, min_ler)
    p_th, err = min(((t, e) for t, e in ((z, z_err), (x, x_err)) if t > 0), default=(0.0, 0.0))
    return dict(p_th=p_th, p_th_error=err, p_th_Z=z, p_th_X=x)


def _indexed(item: tuple) -> tuple:
    i, task = item
    return i, _scan(task)


def _run_tasks(pool: Pool, tasks: list[dict]) -> list:
    """`_scan` of every task, in order; each point it simulates lands in the CSV as soon as it is done."""
    out = [None] * len(tasks)
    with tqdm(total=len(tasks), unit="dir", desc=tasks[0]["scheme"] if tasks else None) as bar:
        for i, lf in pool.imap_unordered(_indexed, enumerate(tasks)):
            out[i] = lf
            bar.set_postfix_str(f"last: {tasks[i]['basis']} memory, threshold "
                                f"{100 * _threshold_from_log_fail(lf, MIN_ERRORS / tasks[i]['shots'])[0]:.3f}%")
            bar.update()
    return out


def _direction(p_g: float, p_t: float, p_r: float) -> tuple[float, float, float]:
    """(theta, phi, length) of the point (p_G, p_T, p_R)."""
    return np.arctan2(p_r, np.hypot(p_g, p_t)), np.arctan2(p_t, p_g), float(np.sqrt(p_g ** 2 + p_t ** 2 + p_r ** 2))


def _run_both(pool: Pool, tasks: list[dict]) -> list[dict]:
    """Each task in the Z and the X basis, `_combine`d."""
    lfs = _run_tasks(pool, [dict(task, basis=b) for task in tasks for b in "ZX"])
    return [_combine(z, x, MIN_ERRORS / task["shots"]) for task, z, x in zip(tasks, lfs[::2], lfs[1::2])]


def _run(args, pool: Pool, bias_only: bool = False) -> dict:
    base = dict(scheme=args.scheme, code=args.code, distances=args.distances, shots=args.shots, batch=args.batch,
                max_fail=args.max_fail, max_fail_rate=args.max_fail_rate, decoder=args.decoder,
                retries=args.retries, cache_only=args.plot_only)
    eta = dict(eta_g=args.eta_g, eta_t=args.eta_t)
    axes = [args.pg, args.pt, args.pr]
    if None in axes:  # coarse scan along each axis for where to aim the grid
        coarse = list(np.geomspace(1e-3, 0.3, 16))
        tasks = [dict(base, **eta, theta=t, phi=f, errors=coarse) for t, f in AXES]
        print("coarse scan along the three axes", flush=True)
        found = [r["p_th"] for r in _run_both(pool, tasks)]
        axes = [a if a is not None else (f or 0.05) for a, f in zip(axes, found)]
        print("axis thresholds p_G, p_T, p_R:", ", ".join(f"{a:.4f}" for a in axes), flush=True)
    p_g_max, p_t_max, p_r_max = axes

    # (b): IBM's grid, (1-s-t) p_G + s p_T + t p_R over s + t <= 1, each scanned around its linear guess
    grid = [(max(0.0, 1 - s - t) * p_g_max, s * p_t_max, t * p_r_max)
            for s in np.linspace(0, 1, args.nphi) for t in np.linspace(0, 1, args.nphi) if s + t <= 1 + 1e-9]
    spread = np.linspace(1 - args.delpth, 1 + args.delpth, args.num_p)
    surface = []
    for p_g, p_t, p_r in grid:
        theta, phi, guess = _direction(p_g, p_t, p_r)
        surface.append(dict(base, **eta, theta=theta, phi=phi, errors=list(guess * spread)))
    # (c): pure gate noise vs eta_G, pure idling vs eta_T, from about 0 to twice the axis threshold
    biases = list(np.logspace(-2, 2, args.nbias))
    wide = np.linspace(1 - args.bias_delpth, 1 + args.bias_delpth, args.bias_num_p or args.num_p)[1:]
    g_bias = [dict(base, eta_g=b, eta_t=args.eta_t, theta=0, phi=0, errors=list(p_g_max * wide)) for b in biases]
    t_bias = [dict(base, eta_g=args.eta_g, eta_t=b, theta=0, phi=np.pi / 2, errors=list(p_t_max * wide))
              for b in biases]

    groups = {"g_bias": g_bias, "t_bias": t_bias} if bias_only else {"surface": surface, "g_bias": g_bias, "t_bias": t_bias}
    print("bias directions" if bias_only else "surface and bias directions", flush=True)
    results = iter(_run_both(pool, [task for tasks in groups.values() for task in tasks]))
    data = dict(settings={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}, axes=axes)
    for key, tasks in groups.items():
        data[key] = []
        for task, result in zip(tasks, results):
            data[key].append(dict(theta=task["theta"], phi=task["phi"], eta_g=task["eta_g"], eta_t=task["eta_t"],
                                  **result))
    return data


def _write_summary(data: dict, path: Path) -> None:
    """`_run`'s data as <name>_thresholds.csv: the three axis thresholds, then a row per direction of each group."""
    s = data["settings"]
    common = dict(scheme=s["scheme"], code=s["code"], decoder=s["decoder"],
                  distances=" ".join(map(str, s["distances"])))
    rows = [dict(common, group="axis", theta=t, phi=f, eta_g=s["eta_g"], eta_t=s["eta_t"], p_th=a)
            for (t, f), a in zip(AXES, data["axes"])]
    # dict(r, ...): rows kept from a saved summary (--bias-only) already carry group and the common columns
    rows += [dict(r, group=g, **common) for g in ("surface", "g_bias", "t_bias") for r in data.get(g, [])]
    write_rows(path, SUMMARY_FIELDS, rows)


def _data_from_summary(rows: list[dict]) -> dict:
    """The rows of <name>_thresholds.csv back as `_run`'s data."""
    axis = [r for r in rows if r["group"] == "axis"]
    first = axis[0]
    settings = dict(scheme=first["scheme"], code=first["code"], decoder=first["decoder"],
                    distances=[int(d) for d in str(first["distances"]).split()], eta_g=first["eta_g"],
                    eta_t=first["eta_t"])
    return dict(settings=settings, axes=[r["p_th"] for r in axis],
                **{g: [r for r in rows if r["group"] == g] for g in ("surface", "g_bias", "t_bias")})


def _plot_surface(ax, rows: list, cmap: str = "copper", alpha: float = 0.5) -> None:
    """IBM's plot_3d_threshold for one logical: threshold points in (p_G, p_T, p_R) in %, the surface through
    them interpolated, both coloured by p_th."""
    pts = [(r["p_th"] * np.cos(r["theta"]) * np.cos(r["phi"]), r["p_th"] * np.cos(r["theta"]) * np.sin(r["phi"]),
            r["p_th"] * np.sin(r["theta"])) for r in rows if r["p_th"] > 0]
    if not pts:  # e.g. a --bias-only run with no surface saved yet
        return
    p_g, p_t, p_r = (100 * np.array(c) for c in zip(*pts))
    p_th = np.sqrt(p_g ** 2 + p_t ** 2 + p_r ** 2)
    norm = colors.Normalize(p_th.min(), p_th.max(), clip=True)
    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    if len(pts) >= 3:
        xi, yi = np.linspace(0, p_g.max(), 200), np.linspace(0, p_t.max(), 200)
        zi = griddata((p_g, p_t), p_r, (xi[None, :], yi[:, None]), method="linear")
        xg, yg = np.meshgrid(xi, yi)
        face = mappable.to_rgba(np.sqrt(xg ** 2 + yg ** 2 + np.nan_to_num(zi) ** 2))
        face[np.isnan(zi)] = (1, 1, 1, 0)
        ax.plot_surface(xg, yg, zi, rstride=1, cstride=1, facecolors=face, alpha=alpha, antialiased=False,
                        linewidth=0, shade=False)
    ax.scatter(p_g, p_t, p_r, c=p_th, cmap=cmap, norm=norm, marker=".")
    cbar = plt.colorbar(mappable, ax=ax, shrink=0.6, pad=0.1)
    cbar.set_label(r"$p_{th}$ [%]", rotation=270, labelpad=13)
    ax.set(xlabel=r"$p_G$ [%]", ylabel=r"$p_T$ [%]", zlabel=r"$p_R$ [%]")
    ax.set_xlim(0)
    ax.set_ylim(0)
    ax.set_zlim(0)
    ax.view_init(elev=20, azim=20)


def _plot(data: dict, png: Path) -> None:
    """(b) the threshold surface with the bias ranges of (c) as arrows on its axes, (c) p_T^th vs eta_T and
    p_G^th vs eta_G (IBM's bias figure), dotted lines at the eta the surface used."""
    s = data["settings"]
    fig = plt.figure(figsize=(11, 4.5))
    ax_b = fig.add_subplot(1, 2, 1, projection="3d")
    _plot_surface(ax_b, data["surface"])
    g = [(r["eta_g"], 100 * r["p_th"]) for r in data["g_bias"] if r["p_th"] > 0]
    t = [(r["eta_t"], 100 * r["p_th"]) for r in data["t_bias"] if r["p_th"] > 0]
    for pts, xyz, color in ((g, lambda lo, hi: ([lo, hi], [0, 0], [0, 0]), "b"),
                            (t, lambda lo, hi: ([0, 0], [lo, hi], [0, 0]), "g")):
        if pts:
            lo, hi = min(p for _, p in pts), max(p for _, p in pts)
            ax_b.add_artist(Arrow3D(*xyz(lo, hi), mutation_scale=10, lw=2, arrowstyle="<|-|>", color=color))
    code = {"css": "", "x3z3": "X$^3$Z$^3$ "}[s.get("code", "css")]
    ax_b.set_title(f"(b) {code}{s['scheme']}, d = {', '.join(map(str, s['distances']))}, "
                   f"$\\eta_G$ = {s['eta_g']:g}, $\\eta_T$ = {s['eta_t']:g}", fontsize=10)

    ax_c = fig.add_subplot(1, 2, 2)
    ax_g = ax_c.twinx()
    handles = []
    if t:
        handles += ax_c.plot(*zip(*t), "o-", c="g", label=r"$p_T^{th}(\eta_T)$")
    if g:
        handles += ax_g.plot(*zip(*g), "o-", c="b", label=r"$p_G^{th}(\eta_G)$")
    ax_c.semilogx()
    ax_c.axvline(s["eta_t"], ls="dotted", c="g")
    ax_c.axvline(s["eta_g"], ls="dotted", c="b")
    ax_c.set_ylim(0, 1.2 * max([p for _, p in t] or [1]))
    ax_g.set_ylim(0, 1.2 * max([p for _, p in g] or [1]))
    ax_c.tick_params(axis="y", colors="g")
    ax_g.tick_params(axis="y", colors="b")
    ax_c.set_xlabel(r"error bias $\eta$")
    ax_c.set_ylabel(r"$p_T^{th}$ [%]", c="g")
    ax_g.set_ylabel(r"$p_G^{th}$ [%]", c="b")
    ax_c.legend(handles=handles, loc="lower left")
    ax_c.set_title("(c)", fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
