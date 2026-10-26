"""The X3Z3 Floquet code against the CSS Floquet code, for every readout scheme, vs noise bias eta: Setiawan &
McLauchlan's (arXiv:2411.04974) Figs. 6 and 7 for these circuits. The sweep is x3z3_ler.py's (see its docstring),
plus --p-fixed in the p grid; with the same settings both write and reuse the same CSV.

    .venv/bin/python src/experiments/x3z3_bias.py --noise spin --distances 3 5 --etas 0.5 1 10 100 1000
    .venv/bin/python src/experiments/x3z3_bias.py --plot-only --out <out>    # redraw from <out>.csv

    <out>_bias.png  (a) threshold vs eta: where the largest distance's p_L curve crosses the smallest's, the lower
                    of the Z and X memories', as paper Fig. 6; (b), (c), ... p_L vs eta at
                    --p-fixed, one panel per distance on a shared axis, as Fig. 7. Both with 95% bands from the shot noise. Colour = scheme, bright and filled X3Z3, pale and hollow CSS.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import floquet
from floquet import CODES, code_distance_png, distance_label, rounds_label
from physical_to_logical import SCHEMES, STYLES
from x3z3_ler import CODE_LABEL, blocks, combined, save, style_axes, sweep, sweep_parser

# With no distance to colour, colour is the scheme: bright for X3Z3, pale for CSS (physical_to_logical.py's bright /
# pale distance pairs).
SCHEME_COLOUR = {"pairs": ("#1f77b4", "#8fb4e8"), "method_a": ("#d62728", "#f0928f"),
                 "method_b": ("#2ca02c", "#8fce8f"), "method_c": ("#9467bd", "#c5b0d5"),
                 "pairs_bell": ("#ff7f0e", "#ffbb78")}  # (X3Z3, CSS)


def main() -> None:
    parser = sweep_parser(__doc__)
    parser.add_argument("--p-fixed", type=float, default=1e-3, help="p of (b), added to the p grid")
    parser.add_argument("--largest-only", action="store_true",
                        help="p_L panel for the largest distance only, into <out>_bias_largest.png")
    args = parser.parse_args()
    rows = sweep(args, {args.p_fixed})
    floquet.LABEL_CODE_DISTANCE = args.code_distance
    png = code_distance_png(Path(f"{args.out}_bias{'_largest' if args.largest_only else ''}.png"))
    _bias_figure(rows, args.p_fixed, png, args.largest_only)
    print(f"wrote {png}")


def _threshold(rows: list[dict], code: str, scheme: str, eta: float) -> tuple[float, float, float]:
    """(p_th, lo, hi) of the lower of the Z and X memories' thresholds: the p where log(p_L) of the largest
    distance first rises above the smallest's, interpolated in log p; nan if they never cross inside the grid.
    lo / hi are its 95% interval (+-1.96 sigma): each point's binomial sigma(log p_L) ~ sqrt((1 - p_L) / errors),
    carried through the interpolation, to match the Clopper-Pearson bands of (b)."""
    found = []
    for basis in "ZX":
        mine = [r for r in rows if (r["code"], r["scheme"], r["basis"]) == (code, scheme, basis)
                and float(r["eta"]) == eta and int(r["errors"])]
        ds = sorted({int(r["distance"]) for r in mine})
        if len(ds) < 2:
            continue
        ler = {(int(r["distance"]), float(r["p"])): (int(r["errors"]), int(r["shots"])) for r in mine}
        ps = sorted(p for p in {p for _, p in ler} if (ds[0], p) in ler and (ds[-1], p) in ler)
        log = lambda d, p: np.log(ler[d, p][0] / ler[d, p][1])
        var = lambda d, p: (1 - ler[d, p][0] / ler[d, p][1]) / ler[d, p][0]
        gap = [log(ds[-1], p) - log(ds[0], p) for p in ps]
        for i in range(len(ps) - 1):
            if gap[i] < 0 <= gap[i + 1]:
                g0, g1, span = gap[i], gap[i + 1], np.log(ps[i + 1]) - np.log(ps[i])
                t = g0 / (g0 - g1)
                v0, v1 = (var(ds[-1], ps[k]) + var(ds[0], ps[k]) for k in (i, i + 1))
                sigma = 1.96 * span * np.hypot(g1 * np.sqrt(v0), g0 * np.sqrt(v1)) / (g0 - g1) ** 2  # d(t)/d(gap)
                found.append(float(np.log(ps[i]) + t * span) + np.array([0, -sigma, sigma]))
                break
    if not found:
        return (float("nan"),) * 3
    return tuple(float(x) for x in np.exp(min(found, key=lambda f: f[0])))


def _combos(rows: list[dict]) -> list[tuple[str, str]]:
    return [(c, s) for c in CODES for s in SCHEMES if any((r["code"], r["scheme"]) == (c, s) for r in rows)]


def _bias_figure(rows: list[dict], p_fixed: float, png: Path, largest_only: bool = False) -> None:
    """Paper Figs. 6 and 7: (a) threshold vs eta, then p_L at p_fixed vs eta, one panel per distance with a shared
    y axis. Line style, marker and colour = scheme; X3Z3 bright and filled, CSS pale and hollow."""
    pl = combined(rows)
    etas = sorted({float(r["eta"]) for r in rows})
    distances = sorted({int(r["distance"]) for r in rows})[-1 if largest_only else 0:]
    fig, (a, *bs) = plt.subplots(1, 1 + len(distances), figsize=(5.5 * (1 + len(distances)), 4.4))
    fig.subplots_adjust(wspace=0.3, right=0.88)
    for b in bs[1:]:
        b.sharey(bs[0])
    for code, scheme in _combos(rows):
        ls, marker = STYLES[scheme]
        colour = SCHEME_COLOUR[scheme][code == "css"]
        style = dict(ls=ls, marker=marker, ms=5, lw=1.5, color=colour, mec=colour,
                     mfc=colour if code == "x3z3" else "white")
        th = np.array([_threshold(rows, code, scheme, eta) for eta in etas])
        a.plot(etas, th[:, 0], **style)
        a.fill_between(etas, th[:, 1], th[:, 2], color=colour, alpha=0.12, lw=0)
        for d, b in zip(distances, bs):
            pts = [(eta, pl[k]) for eta in etas if (k := (code, scheme, d, eta, p_fixed)) in pl]
            if pts:
                b.plot([e for e, _ in pts], [v[0] for _, v in pts], **style)
                b.fill_between([e for e, _ in pts], [v[1] for _, v in pts], [v[2] for _, v in pts],
                               color=colour, alpha=0.12, lw=0)
    a.set(xscale="log", xlabel="noise bias η", ylabel="threshold $p_{th}$", title="(a) threshold")
    for i, (d, b) in enumerate(zip(distances, bs)):
        b.set(xscale="log", yscale="log", xlabel="noise bias η", title=f"({'bcdefg'[i]}) d = {distance_label(d)}, {rounds_label(str(d))}, p = {p_fixed:g}")
    bs[0].set_ylabel("logical failure probability $p_L$")
    for ax in (a, *bs):
        style_axes(ax)
    present = _combos(rows)
    blocks(bs[-1], [("Schemes", [Line2D([], [], marker=STYLES[s][1], ls=STYLES[s][0], color=SCHEME_COLOUR[s][0], label=s)
                           for s in dict.fromkeys(s for _, s in present)]),
              ("Codes", [Line2D([], [], marker="o", ls="-", color=g, mec=g, mfc=g if c == "x3z3" else "white",
                                label=f"{CODE_LABEL[c]} ({'bright' if c == 'x3z3' else 'pale'})")
                         for c, g in (("x3z3", "0.25"), ("css", "0.7")) if any(c == x for x, _ in present)])])
    fig.suptitle(f"{rows[0]['noise']} noise")
    save(fig, png)


if __name__ == "__main__":
    main()
