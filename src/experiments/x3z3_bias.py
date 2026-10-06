"""The X3Z3 Floquet code against the CSS Floquet code, for every readout scheme, vs noise bias eta: Setiawan &
McLauchlan's (arXiv:2411.04974) Figs. 6 and 7 for these circuits. The sweep is x3z3_ler.py's (see its docstring),
plus --p-fixed in the p grid; with the same settings both write and reuse the same CSV.

    .venv/bin/python src/experiments/x3z3_bias.py --noise spin --distances 3 5 --etas 0.5 1 10 100 1000
    .venv/bin/python src/experiments/x3z3_bias.py --plot-only --out <out>    # redraw from <out>.csv

    <out>_bias.png  (a) threshold vs eta: where the largest distance's p_L curve crosses the smallest's, the lower
                    of the Z and X memories', as paper Fig. 6; (b) p_L vs eta at --p-fixed for the largest distance,
                    as Fig. 7. Colour = scheme, bright and filled X3Z3, pale and hollow CSS.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from floquet import CODES
from physical_to_logical import SCHEMES, STYLES
from x3z3_ler import CODE_LABEL, blocks, combined, save, style_axes, sweep, sweep_parser

# With no distance to colour, colour is the scheme: bright for X3Z3, pale for CSS (physical_to_logical.py's bright /
# pale distance pairs).
SCHEME_COLOUR = {"pairs": ("#1f77b4", "#8fb4e8"), "method_a": ("#d62728", "#f0928f"),
                 "method_b": ("#2ca02c", "#8fce8f"), "method_c": ("#9467bd", "#c5b0d5")}  # (X3Z3, CSS)


def main() -> None:
    parser = sweep_parser(__doc__)
    parser.add_argument("--p-fixed", type=float, default=1e-3, help="p of (b), added to the p grid")
    args = parser.parse_args()
    rows = sweep(args, {args.p_fixed})
    png = Path(f"{args.out}_bias.png")
    _bias_figure(rows, args.p_fixed, png)
    print(f"wrote {png}")


def _threshold(rows: list[dict], code: str, scheme: str, eta: float) -> float:
    """Lower of the Z and X memories' thresholds: the p where log(p_L) of the largest distance first rises
    above the smallest's, interpolated in log p; nan if they never cross inside the grid."""
    found = []
    for basis in "ZX":
        mine = [r for r in rows if (r["code"], r["scheme"], r["basis"]) == (code, scheme, basis)
                and float(r["eta"]) == eta and int(r["errors"])]
        ds = sorted({int(r["distance"]) for r in mine})
        if len(ds) < 2:
            continue
        ler = {(int(r["distance"]), float(r["p"])): int(r["errors"]) / int(r["shots"]) for r in mine}
        ps = sorted(p for p in {p for _, p in ler} if (ds[0], p) in ler and (ds[-1], p) in ler)
        gap = [np.log(ler[ds[-1], p]) - np.log(ler[ds[0], p]) for p in ps]
        for i in range(len(ps) - 1):
            if gap[i] < 0 <= gap[i + 1]:
                t = gap[i] / (gap[i] - gap[i + 1])
                found.append(float(np.exp(np.log(ps[i]) + t * (np.log(ps[i + 1]) - np.log(ps[i])))))
                break
    return min(found, default=float("nan"))


def _combos(rows: list[dict]) -> list[tuple[str, str]]:
    return [(c, s) for c in CODES for s in SCHEMES if any((r["code"], r["scheme"]) == (c, s) for r in rows)]


def _bias_figure(rows: list[dict], p_fixed: float, png: Path) -> None:
    """Paper Figs. 6 and 7: (a) threshold vs eta, (b) p_L at p_fixed vs eta, largest distance. Line style, marker
    and colour = scheme; X3Z3 bright and filled, CSS pale and hollow."""
    pl = combined(rows)
    etas = sorted({float(r["eta"]) for r in rows})
    d_max = max(int(r["distance"]) for r in rows)
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.4))
    fig.subplots_adjust(wspace=0.3, right=0.85)
    for code, scheme in _combos(rows):
        ls, marker = STYLES[scheme]
        colour = SCHEME_COLOUR[scheme][code == "css"]
        style = dict(ls=ls, marker=marker, ms=5, lw=1.5, color=colour, mec=colour,
                     mfc=colour if code == "x3z3" else "white")
        a.plot(etas, [_threshold(rows, code, scheme, eta) for eta in etas], **style)
        pts = [(eta, pl[k]) for eta in etas if (k := (code, scheme, d_max, eta, p_fixed)) in pl]
        if pts:
            b.plot([e for e, _ in pts], [v[0] for _, v in pts], **style)
            b.fill_between([e for e, _ in pts], [v[1] for _, v in pts], [v[2] for _, v in pts],
                           color=colour, alpha=0.12, lw=0)
    a.set(xscale="log", xlabel="noise bias η", ylabel="threshold $p_{th}$", title="(a) threshold")
    b.set(xscale="log", yscale="log", xlabel="noise bias η", ylabel="logical failure probability $p_L$",
          title=f"(b) d = {d_max}, {d_max} rounds, p = {p_fixed:g}")
    for ax in (a, b):
        style_axes(ax)
    present = _combos(rows)
    blocks(b, [("Schemes", [Line2D([], [], marker=STYLES[s][1], ls=STYLES[s][0], color=SCHEME_COLOUR[s][0], label=s)
                           for s in dict.fromkeys(s for _, s in present)]),
              ("Codes", [Line2D([], [], marker="o", ls="-", color=g, mec=g, mfc=g if c == "x3z3" else "white",
                                label=f"{CODE_LABEL[c]} ({'bright' if c == 'x3z3' else 'pale'})")
                         for c, g in (("x3z3", "0.25"), ("css", "0.7")) if any(c == x for x, _ in present)])])
    fig.suptitle(f"{rows[0]['noise']} noise")
    save(fig, png)


if __name__ == "__main__":
    main()
