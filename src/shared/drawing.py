"""Drawing helpers both readout schemes use: the colour palette, PNG export, and the stim timeslice-svg
post-processing (coloured plaquettes underneath, qubits styled by kind and labelled)."""
import hashlib
import re
from pathlib import Path

import cairosvg

# Categorical slots in fixed order: one per distance, or per plaquette colour (dataviz reference palette).
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
# How each kind of qubit is drawn, as in SpinQEC's layout plots: fill and legend label. All get a dark outline.
QUBIT_KINDS = {"data": ("#333333", "data"), "main": ("white", "ancilla: syndrome"),
               "reff": ("#cfcfcf", "ancilla: reference")}


def svg_png(svg: str, png: Path) -> None:
    """2x scale, less when that would pass cairo's 32767 px limit on either side."""
    w, h = map(float, re.search(r'viewBox="[\d.]+ [\d.]+ ([\d.]+) ([\d.]+)"', svg).groups())
    cairosvg.svg2png(bytestring=svg.encode(), write_to=str(png), scale=min(2, 32000 / max(w, h)),
                     background_color="white")


def window(points: list[complex], torus_period: complex) -> complex:
    """Top-left corner of the one-torus-period window centred on `points`."""
    xs, ys = [p.real for p in points], [p.imag for p in points]
    return complex(min(xs) - (torus_period.real - (max(xs) - min(xs))) / 2,
                   min(ys) - (torus_period.imag - (max(ys) - min(ys))) / 2)


def color_hexes(svg: str, centers: dict[complex, int], corners: list[complex], i2pos: dict[int, complex],
                torus_period: complex) -> str:
    """Draws colour-filled plaquettes underneath a stim timeslice-svg diagram."""
    # ponytail: stim has no colour option for timeslices, so recover coordinate -> pixel from its qubit
    # dot ids. Breaks if stim changes its svg ids; the assert below says so.
    panels = [tuple(map(float, r)) for r in re.findall(
        r'id="tick_border:[^"]*" x="([-\d.]+)" y="([-\d.]+)" width="([-\d.]+)" height="([-\d.]+)"', svg)]
    # a circuit without TICKs is one borderless panel: the whole svg
    panels = panels or [(0, 0, *map(float, re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg).groups()))]
    # dot ids only carry the panel column, so find each dot's panel by position
    dots = [(i2pos[int(q)], next(k for k, (px, py, w, h) in enumerate(panels)
                                 if px <= float(x) <= px + w and py <= float(y) <= py + h), float(x), float(y))
            for q, x, y in re.findall(r'id="qubit_dot:(\d+):[^"]*" cx="([-\d.]+)" cy="([-\d.]+)"', svg)]
    assert dots and panels, "stim changed its timeslice svg format"
    xs = sorted((p.real, x) for p, t, x, _ in dots if t == 0)
    ys = sorted((p.imag, y) for p, t, _, y in dots if t == 0)
    sx = (xs[-1][1] - xs[0][1]) / (xs[-1][0] - xs[0][0])
    sy = (ys[-1][1] - ys[0][1]) / (ys[-1][0] - ys[0][0])
    origin = {t: complex(x - sx * p.real, y - sy * p.imag) for p, t, x, y in dots}
    # one torus period around the qubits, so wrapped hexes show up once on each side
    corner = window([p for p, t, _, _ in dots if t == 0], torus_period)
    left, top = corner.real, corner.imag
    shifts = [a * torus_period.real + b * torus_period.imag * 1j for a in (-1, 0, 1) for b in (-1, 0, 1)]

    # ids are page-wide once several SVGs share an HTML page (a notebook), so make them unique per diagram
    tag = hashlib.sha1(svg.encode()).hexdigest()[:10]
    out = []
    for t, o in origin.items():
        out.append(f'<clipPath id="hexclip{tag}_{t}"><rect x="{o.real + sx * left}" y="{o.imag + sy * top}" '
                   f'width="{sx * torus_period.real}" height="{sy * torus_period.imag}"/></clipPath>'
                   f'<g clip-path="url(#hexclip{tag}_{t})">')
        for h, colour in centers.items():
            for s in shifts:
                pts = " ".join(f"{o.real + sx * (h + s + d).real},{o.imag + sy * (h + s + d).imag}" for d in corners)
                out.append(f'<polygon points="{pts}" fill="{COLORS[colour]}" fill-opacity="0.3" stroke="#666"/>')
        out.append("</g>")
    return svg.replace("\n", "\n" + "\n".join(out) + "\n", 1)  # right after <svg>, i.e. below everything


def style_qubits(svg: str, kinds: list[str]) -> str:
    """Stim draws every qubit as the same small dot: redraw each as its kind (`kinds[i]` for qubit i)."""
    def dot(m: re.Match[str]) -> str:
        fill = QUBIT_KINDS[kinds[int(m[1])]][0]
        return (f'<circle id="qubit_dot:{m[1]}:{m[2]}" cx="{m[3]}" cy="{m[4]}" r="5" fill="{fill}" '
                f'stroke="#333333" stroke-width="1.5"/>')
    styled, n = re.subn(r'<circle id="qubit_dot:(\d+):([^"]*)" cx="([-\d.]+)" cy="([-\d.]+)" r="2" '
                        r'stroke="none" fill="black"/>', dot, svg)
    assert n, "stim changed its timeslice svg format"
    return styled


def label_qubits(svg: str, labels: list[str]) -> str:
    """Writes `labels[i]` beside every dot of qubit i in a timeslice diagram (data bold, as in the layouts).
    They go in right after the dots, so gates, drawn later, stay on top."""
    def text(m: re.Match[str]) -> str:
        label = labels[int(m[1])]
        return (f'{m[0]}<text x="{float(m[2]) + 6}" y="{float(m[3]) - 5}" font-size="9" font-family="sans-serif" '
                f'fill="#222222" font-weight="{"bold" if label.startswith("D") else "normal"}">{label}</text>')
    labelled, n = re.subn(r'<circle id="qubit_dot:(\d+):[^"]*" cx="([-\d.]+)" cy="([-\d.]+)"[^>]*/>', text, svg)
    assert n, "stim changed its timeslice svg format"
    return labelled


def drawing_dir(out: Path, hex_view: bool, numbers: bool) -> Path:
    """Where a drawing goes: out/rectangle or out/hex (brick-wall or regular-hexagon view), then numbered or
    plain. Made if missing."""
    folder = out / ("hex" if hex_view else "rectangle") / ("numbered" if numbers else "plain")
    folder.mkdir(parents=True, exist_ok=True)
    return folder
