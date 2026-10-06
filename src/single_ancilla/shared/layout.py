"""The one-ancilla lattice every single-ancilla method uses: one spin ancilla at the midpoint of each edge,
where it is drawn, and drawing each tick's qubits with their own kinds."""
import re

from drawing import label_qubits, style_qubits
import floquet
from floquet import hex_positions, qubits, torus


def edge_list(distance: int) -> list[tuple[int, complex, complex, complex]]:
    """Every edge once, in floquet.edge_list order: (colour, data qubit, data qubit, ancilla), the ancilla at the
    edge's midpoint."""
    return [(e.colour, *e.data, torus((e.ends[0] + e.ends[1]) / 2, distance)) for e in floquet.edge_list(distance)]


def qubit_kinds(distance: int) -> dict[complex, str]:
    """Qubit -> "data" or "main" (the ancilla; drawn as plots.py draws a syndrome ancilla), in index order:
    data qubits first, then the ancilla of each edge in `edge_list` order."""
    kinds = {q: "data" for q in qubits(distance)}
    edges = edge_list(distance)
    kinds.update({anc: "main" for *_, anc in edges})
    assert len(kinds) == len(qubits(distance)) + len(edges), "two ancillas on one spot"
    return kinds


def positions(distance: int, hex_view: bool) -> dict[complex, complex]:
    """Qubit -> where it is drawn: its brick-wall coordinates, or on regular hexagons (ancillas still halfway)."""
    if not hex_view:
        return {q: q for q in qubit_kinds(distance)}
    pos = hex_positions(distance)
    for e in floquet.edge_list(distance):
        pos[torus((e.ends[0] + e.ends[1]) / 2, distance)] = torus((e.hex_ends[0] + e.hex_ends[1]) / 2, distance)
    return pos


def style_ticks(svg: str, kinds_per_tick: list[list[str]], labels_per_tick: list[list[str]] | None = None) -> str:
    """plots.style_qubits, and plots.label_qubits given labels, with their own kinds and labels for each tick:
    qubit i in tick t drawn as kinds_per_tick[t][i]."""
    # ponytail: stim's dot ids carry the panel column, not the tick, so find each dot's panel by position
    panels = {int(t): tuple(map(float, r)) for t, *r in re.findall(
        r'id="tick_border:(\d+):[^"]*" x="([-\d.]+)" y="([-\d.]+)" width="([-\d.]+)" height="([-\d.]+)"', svg)}

    def dot(m: re.Match[str]) -> str:
        x, y = float(m[1]), float(m[2])
        t = next(t for t, (px, py, w, h) in panels.items() if px <= x <= px + w and py <= y <= py + h)
        styled = style_qubits(m[0], kinds_per_tick[t])
        return label_qubits(styled, labels_per_tick[t]) if labels_per_tick else styled
    return re.sub(r'<circle id="qubit_dot:\d+:[^"]*" cx="([-\d.]+)" cy="([-\d.]+)"[^>]*/>', dot, svg)
