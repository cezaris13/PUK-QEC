"""The figures of x3z3.py (see its docstring), drawn from its CSV: x3z3.py --plot-only redraws them.

    <out>_bias.png             threshold and p_L at --p-fixed vs eta (paper Figs. 6, 7)
    <out>_ler_eta<eta>.png     p_L vs p, one panel per scheme, both codes
    <out>_qubits_eta<eta>.png  p_L at --p-fixed vs total physical qubits, with fits
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from decode import DISTANCE_COLOUR, SCHEMES, STYLES, clopper_pearson  # also puts src/shared on the path
from floquet import CODES

# As decode.py and footprint.py: colour = distance (DISTANCE_COLOUR), line style and marker = scheme (STYLES), fits in
# grey. The code is the marker fill: X3Z3 filled, CSS hollow and faded. Where there is no distance to colour (the bias
# figure), colour is the scheme, bright for X3Z3 and pale for CSS (decode.py's bright / pale distance pairs).
CODE_LABEL = {"x3z3": "X$^3$Z$^3$", "css": "CSS"}
SCHEME_COLOUR = {"pairs": ("#1f77b4", "#8fb4e8"), "method_a": ("#d62728", "#f0928f"),
                 "method_b": ("#2ca02c", "#8fce8f"), "method_c": ("#9467bd", "#c5b0d5")}  # (X3Z3, CSS)
CODE_ALPHA = {"x3z3": 1.0, "css": 0.55}
FIT_GREY = {"x3z3": "0.35", "css": "0.7"}

def interval(errors: int, shots: int) -> tuple[float, float, float]:
    """(rate, lower, upper), 95% Clopper-Pearson."""
    return (errors / shots, *clopper_pearson(errors, shots))


def combined(rows: list[dict]) -> dict:
    """{(code, scheme, distance, eta, p): (p_L, lower, upper, qubits)}: the Z and X memories as Eq. (9),
    each end of the interval combined the same way; points where neither basis saw an error are left out."""
    by = {}
    for r in rows:
        by.setdefault((r["code"], r["scheme"], int(r["distance"]), float(r["eta"]), float(r["p"])), {})[r["basis"]] = r
    out = {}
    for k, bases in by.items():
        if len(bases) < 2 or not any(int(b["errors"]) for b in bases.values()):
            continue
        z, x = (interval(int(bases[b]["errors"]), int(bases[b]["shots"])) for b in "ZX")
        out[k] = (*(1 - (1 - a) * (1 - b) for a, b in zip(z, x)), int(bases["Z"]["qubits"]))
    return out


def threshold(rows: list[dict], code: str, scheme: str, eta: float) -> float:
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


def combos(rows: list[dict]) -> list[tuple[str, str]]:
    return [(c, s) for c in CODES for s in SCHEMES if any((r["code"], r["scheme"]) == (c, s) for r in rows)]


def mark(code: str, colour: str) -> dict:
    """Marker fill and line opacity for `code`: X3Z3 filled, CSS hollow and faded."""
    return dict(color=colour, mfc=colour if code == "x3z3" else "white", mec=colour, alpha=CODE_ALPHA[code])


def decode_axes(ax) -> None:
    """decode.py's plot styling."""
    ax.grid(which="major", color="#e4e4e0", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)


def blocks(ax, titled: list[tuple[str, list[Line2D]]], top: float = 1.0) -> None:
    """Legend blocks stacked right of `ax`, as decode.legend_blocks: (title, handles), empty ones skipped."""
    y = top
    for title, handles in titled:
        if handles:
            ax.add_artist(ax.legend(handles=handles, title=title, fontsize=7, title_fontsize=8, loc="upper left",
                                    bbox_to_anchor=(1.02, y), frameon=False))
            y -= 0.08 + 0.055 * len(handles)


def scheme_handles(schemes: list[str], lines: bool) -> list[Line2D]:
    return [Line2D([], [], marker=STYLES[s][1], ls=STYLES[s][0] if lines else "none", color="0.3", label=s)
            for s in schemes]


def code_handles(codes: list[str]) -> list[Line2D]:
    return [Line2D([], [], marker="o", ls="none", label=CODE_LABEL[c], **mark(c, "0.3")) for c in codes]


def distance_handles(distances: list[int]) -> list[Line2D]:
    return [Line2D([], [], marker="s", ls="none", color=DISTANCE_COLOUR.get(d, "0.3"), label=str(d)) for d in distances]


def save(fig, png: Path) -> None:
    """Save cropped to its content: bbox_inches="tight" alone leaves out legends added with add_artist, and then
    the suptitle too once extra artists are given."""
    extra = [a for ax in fig.axes for a in ax.artists] + [t for t in [fig._suptitle] if t is not None]
    fig.savefig(png, dpi=200, bbox_inches="tight", bbox_extra_artists=extra)
    plt.close(fig)


def bias_figure(rows: list[dict], p_fixed: float, png: Path) -> None:
    """Paper Figs. 6 and 7: (a) threshold vs eta, (b) p_L at p_fixed vs eta, largest distance. Line style, marker
    and colour = scheme; X3Z3 bright and filled, CSS pale and hollow."""
    pl = combined(rows)
    etas = sorted({float(r["eta"]) for r in rows})
    d_max = max(int(r["distance"]) for r in rows)
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.4))
    fig.subplots_adjust(wspace=0.3, right=0.85)
    for code, scheme in combos(rows):
        ls, marker = STYLES[scheme]
        colour = SCHEME_COLOUR[scheme][code == "css"]
        style = dict(ls=ls, marker=marker, ms=5, lw=1.5, color=colour, mec=colour,
                     mfc=colour if code == "x3z3" else "white")
        a.plot(etas, [threshold(rows, code, scheme, eta) for eta in etas], **style)
        pts = [(eta, pl[k]) for eta in etas if (k := (code, scheme, d_max, eta, p_fixed)) in pl]
        if pts:
            b.plot([e for e, _ in pts], [v[0] for _, v in pts], **style)
            b.fill_between([e for e, _ in pts], [v[1] for _, v in pts], [v[2] for _, v in pts],
                           color=colour, alpha=0.12, lw=0)
    a.set(xscale="log", xlabel="noise bias η", ylabel="threshold $p_{th}$", title="(a) threshold")
    b.set(xscale="log", yscale="log", xlabel="noise bias η", ylabel="logical failure probability $p_L$",
          title=f"(b) d = {d_max}, {d_max} rounds, p = {p_fixed:g}")
    for ax in (a, b):
        decode_axes(ax)
    present = combos(rows)
    blocks(b, [("Schemes", [Line2D([], [], marker=STYLES[s][1], ls=STYLES[s][0], color=SCHEME_COLOUR[s][0], label=s)
                            for s in dict.fromkeys(s for _, s in present)]),
               ("Codes", [Line2D([], [], marker="o", ls="-", color=g, mec=g, mfc=g if c == "x3z3" else "white",
                                 label=f"{CODE_LABEL[c]} ({'bright' if c == 'x3z3' else 'pale'})")
                          for c, g in (("x3z3", "0.25"), ("css", "0.7")) if any(c == x for x, _ in present)])])
    fig.suptitle(f"{rows[0]['noise']} noise")
    save(fig, png)


def ler_figure(rows: list[dict], eta: float, png: Path) -> None:
    """decode.py's LER vs p, one panel per scheme with both codes: colour = distance, the scheme's line style and
    marker, X3Z3 filled, CSS hollow and faded, shaded 95% Clopper-Pearson band."""
    pl = combined(rows)
    schemes = [s for s in SCHEMES if any(r["scheme"] == s for r in rows)]
    distances = sorted({int(r["distance"]) for r in rows})
    fig, axes = plt.subplots(1, len(schemes), figsize=(4 * len(schemes) + 1.5, 4.2), sharey=True, squeeze=False)
    fig.subplots_adjust(wspace=0.08, right=0.88)
    for ax, scheme in zip(axes[0], schemes):
        ls, marker = STYLES[scheme]
        for code in CODES:
            for d in distances:
                pts = sorted((k[4], v) for k, v in pl.items() if k[:4] == (code, scheme, d, eta))
                if not pts:
                    continue
                ps, vs = [p for p, _ in pts], [v for _, v in pts]
                colour = DISTANCE_COLOUR.get(d, "0.3")
                ax.plot(ps, [v[0] for v in vs], ls=ls, marker=marker, ms=4, lw=1.5, **mark(code, colour))
                ax.fill_between(ps, [v[1] for v in vs], [v[2] for v in vs], color=colour,
                                alpha=0.2 * CODE_ALPHA[code], lw=0)
        ax.axline((1e-3, 1e-3), (1e-2, 1e-2), color="0.5", ls=":", lw=1)
        ax.loglog()
        ax.set(xlabel="physical error rate p", title=scheme)
        decode_axes(ax)
    axes[0][0].set_ylabel("logical failure probability $p_L$ (d rounds)")
    blocks(axes[0][-1], [("Schemes", scheme_handles(schemes, lines=True)),
                         ("Codes", code_handles([c for c in CODES if any(r["code"] == c for r in rows)])),
                         ("Distances", distance_handles(distances))])
    fig.suptitle(f"{rows[0]['noise']} noise, η = {eta:g}; dotted: LER = p")
    save(fig, png)


def qubits_figure(rows: list[dict], eta: float, p_fixed: float, png: Path) -> None:
    """footprint.py's figure: p_L at p_fixed vs total physical qubits, square-root x / log y axes, colour =
    distance with 95% Clopper-Pearson bars, marker = scheme (X3Z3 filled, CSS hollow), and for each (scheme, code)
    with two distances the fit log p_L = a + b sqrt(qubits) in grey (lighter for CSS), the scheme's line style."""
    pl = combined(rows)
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    fig.subplots_adjust(left=0.10, right=0.70)
    fits = []
    for code, scheme in combos(rows):
        ls, marker = STYLES[scheme]
        pts = sorted((v[3], v, k[2]) for k, v in pl.items() if k[:2] == (code, scheme) and k[3:] == (eta, p_fixed))
        for n, (c, lo, hi, _), d in pts:
            ax.errorbar(n, c, yerr=[[c - lo], [hi - c]], marker=marker, ms=6, lw=1, capsize=2, ls="none",
                        **mark(code, DISTANCE_COLOUR.get(d, "0.3")))
        if len(pts) >= 2:
            xs, ys = [n for n, *_ in pts], [v[0] for _, v, _ in pts]
            slope, intercept = np.polyfit(np.sqrt(xs), np.log(ys), 1)
            grid = np.linspace(min(xs) * 0.8, max(xs) * 1.6, 100)
            ax.plot(grid, np.exp(slope * np.sqrt(grid) + intercept), ls=ls, color=FIT_GREY[code], lw=1)
            fits.append(Line2D([], [], ls=ls, color=FIT_GREY[code], label=f"{CODE_LABEL[code]}, {scheme}"))
    ax.set_xscale("function", functions=(np.sqrt, np.square))
    ax.set_yscale("log")
    ax.set(xlabel="Total Physical Qubits", ylabel="Logical Error Rate (d rounds)")
    ax.set_title(f"{rows[0]['noise']} noise, η = {eta:g}, p = {p_fixed:g}", fontsize=10)
    ax.grid(which="both", alpha=0.3, lw=0.5)
    present = combos(rows)
    blocks(ax, [("Schemes", scheme_handles(list(dict.fromkeys(s for _, s in present)), lines=False)),
                ("Codes", code_handles(list(dict.fromkeys(c for c, _ in present)))),
                ("Linefits", fits),
                ("Distances", distance_handles(sorted({int(r["distance"]) for r in rows})))])
    save(fig, png)


def plot(csv_path: Path, out: Path, p_fixed: float) -> list[Path]:
    rows = list(csv.DictReader(csv_path.open()))
    written = [Path(f"{out}_bias.png")]
    bias_figure(rows, p_fixed, written[0])
    for eta in sorted({float(r["eta"]) for r in rows}):
        written += [Path(f"{out}_ler_eta{eta:g}.png"), Path(f"{out}_qubits_eta{eta:g}.png")]
        ler_figure(rows, eta, written[-2])
        qubits_figure(rows, eta, p_fixed, written[-1])
    return written
