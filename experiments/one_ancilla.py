import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import stim
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from floquet import (EDGE_TYPES, HEX_CORNERS, QUBIT_CORNERS, edge_ends, hex_centers, hex_positions, period,
                     qubits, sorted_complex, torus)
from plots import COLORS, QUBIT_KINDS, color_hexes, label_qubits, style_qubits, svg_png, window


def edge_list(distance: int) -> list[tuple[int, complex, complex, complex]]:
    """Every edge once: (colour, data qubit, data qubit, ancilla), the ancilla at the edge's midpoint."""
    out = []
    for h, colour in hex_centers(distance).items():
        for e in EDGE_TYPES:
            a, b = edge_ends(h, e, QUBIT_CORNERS)
            out.append((colour, torus(a, distance), torus(b, distance), torus((a + b) / 2, distance)))
    return out


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
    data = hex_positions(distance)
    pos = {q: data[q] for q in qubits(distance)}
    for h in hex_centers(distance):
        for e in EDGE_TYPES:
            a, b = edge_ends(h, e, QUBIT_CORNERS)
            ha, hb = edge_ends(h, e, HEX_CORNERS)
            pos[torus((a + b) / 2, distance)] = torus((ha + hb) / 2, distance)
    return pos


def ring_swaps(distance: int, colour: int) -> list[tuple[complex, complex]]:
    """(data, ancilla) around every hex of `colour`: each corner with the ancilla before it clockwise (as
    drawn, y down), so SWAPping every pair rotates the data one step counterclockwise around their hex, and
    the same SWAPs again rotate them back clockwise."""
    ancillas = {a for *_, a in edge_list(distance)}
    pairs = []
    for h, c in hex_centers(distance).items():
        if c == colour:
            for i, corner in enumerate(QUBIT_CORNERS):
                anc = torus(h + (QUBIT_CORNERS[i - 1] + corner) / 2, distance)  # mid of the edge before it
                assert anc in ancillas
                pairs.append((torus(h + corner, distance), anc))
    return pairs


def syndrome_checks(distance: int, colour: int) -> list[tuple[int, complex, complex, complex]]:
    """(k, u, v, a) for every other edge u-v around each hex of `colour`: right to bottom-right, bottom-left to
    left, top-left to top-right; sorted by k. After the counterclockwise swaps the ancilla on corner u is the
    syndrome, with u's data moved in on one side and v's on the other, and it reads out against a, the
    ancilla halfway along edge k: the `colour` edge at u, the syndrome's third neighbour. The checked edges
    are one whole colour class, so every data qubit is in exactly one check and every a in at most one."""
    at = {}  # corner -> (k, ancilla) of the `colour` edge ending there
    for k, (c, d0, d1, a) in enumerate(edge_list(distance)):
        if c == colour:
            at[d0] = at[d1] = (k, a)
    out = []
    for h, c in hex_centers(distance).items():
        if c == colour:
            ring = [torus(h + corner, distance) for corner in QUBIT_CORNERS]
            out += [(*at[ring[i]], ring[i], ring[i + 1]) for i in (0, 2, 4)]
    return sorted(((k, u, v, a) for k, a, u, v in out), key=lambda check: check[0])


def corner_check(pauli: str, d0: int, d1: int, syndrome: int, reference: int) -> list[list[tuple[str, list[int]]]]:
    """Measure P_d0 P_d1 onto `syndrome` from both of its sides, then read it out against `reference` by MZZ:
    floquet.parity_check's gates, without its SWAP since both data are next to the syndrome. Z: the syndrome
    starts in |0> and both data CX onto it. X and Y: it starts in |+>, controls a P on each, and H turns its
    X into the Z that MZZ reads. The reference starts in |0>. Layers, for the caller to TICK between."""
    if pauli == "Z":
        return [[("R", [syndrome, reference])], [("CX", [d0, syndrome])], [("CX", [d1, syndrome])],
                [("MZZ", [syndrome, reference])]]
    return [[("RX", [syndrome]), ("R", [reference])], [(f"C{pauli}", [syndrome, d0])],
            [(f"C{pauli}", [syndrome, d1])], [("H", [syndrome])], [("MZZ", [syndrome, reference])]]


def rotated(r: int) -> int:
    """The hexes step r rotates. Rotating the colour c hexes checks the colour c - 1 edges, so the colour r % 3
    edges step r measures (bX rZ gX bZ rX gZ, as in floquet.py) come from rotating colour (r + 1) % 3:
    blue edges from red hexes, red edges from green hexes, green edges from blue hexes."""
    return (r + 1) % 3


def round_circuit(distance: int, r: int, hex_view: bool = False) -> stim.Circuit:
    """Step r, a TICK after every layer:
    1. SWAP the data qubits and ancillas around every `rotated(r)` hex, rotating the data counterclockwise
       (see `ring_swaps`), then an empty TICK so diagrams show where everything is.
    2. On every other corner of those hexes, the syndrome measures the XX (even r) or ZZ (odd r) parity of
       the two data next to it, which together are every colour r % 3 edge, read out against the ancilla of
       the corner's `rotated(r)` edge (see `syndrome_checks`, `corner_check`). Record j is the j-th check's
       parity. The other corners idle.
    3. The same SWAPs again, rotating everything back clockwise."""
    pos = positions(distance, hex_view)
    q2i = {q: i for i, q in enumerate(qubit_kinds(distance))}
    circuit = stim.Circuit()
    for q, i in q2i.items():
        circuit.append("QUBIT_COORDS", [i], [pos[q].real, pos[q].imag])
    pairs = ring_swaps(distance, rotated(r))
    swaps = [q2i[q] for pair in pairs for q in pair]
    moved = dict(pairs)  # data qubit's corner -> the ancilla spot its state is on after the swaps
    circuit.append("SWAP", swaps)
    circuit.append("TICK")
    circuit.append("TICK")
    checks = [corner_check("XZ"[r % 2], q2i[moved[u]], q2i[moved[v]], q2i[u], q2i[a])
              for _, u, v, a in syndrome_checks(distance, rotated(r))]
    for layer in zip(*checks):  # layer k of every check at once
        for part in layer:
            for name, targets in part:
                circuit.append(name, targets)
        circuit.append("TICK")
    circuit.append("SWAP", swaps)
    circuit.append("TICK")
    return circuit


def roles(distance: int, r: int, swapped: bool) -> tuple[list[str], list[str]]:
    """(kinds, labels) by stim index, before step r's swaps or after them. Before: D<i> data qubit i on its
    corner, A<k> the ancilla of edge k (as in the numbered layout). After: data and ancillas traded places
    with their labels, then each check's syndrome S<k> (`main`) on its corner, and the ancilla A<k> of edge k
    drawn as its reference (`reff`), as `round_circuit` uses them (see `syndrome_checks`)."""
    kinds = qubit_kinds(distance)
    q2i = {q: i for i, q in enumerate(kinds)}
    n_data = len(qubits(distance))
    kind = list(kinds.values())
    label = [f"D{i}" if k == "data" else f"A{i - n_data}" for i, k in enumerate(kind)]
    if swapped:
        for d, a in ring_swaps(distance, rotated(r)):
            i, j = q2i[d], q2i[a]
            kind[i], kind[j], label[i], label[j] = kind[j], kind[i], label[j], label[i]
        for k, u, _, a in syndrome_checks(distance, rotated(r)):
            kind[q2i[u]], label[q2i[u]] = "main", f"S{k}"
            kind[q2i[a]] = "reff"
    return kind, label


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


def check(distance: int) -> None:
    """Each colour's swaps pair every data qubit exactly once, each with a different ancilla, and its checks
    are one whole colour class of real edges, each with its own reference; every step reads its j-th check's
    parity into record j and leaves every data qubit back where it started."""
    colours = {frozenset((d0, d1)): c for c, d0, d1, _ in edge_list(distance)}
    for colour in range(3):
        pairs = ring_swaps(distance, colour)
        assert sorted_complex(d for d, _ in pairs) == qubits(distance), colour
        assert len({a for _, a in pairs}) == len(pairs), colour
        checks = syndrome_checks(distance, colour)
        assert sorted_complex(q for _, u, v, _ in checks for q in (u, v)) == qubits(distance), colour
        assert len({colours[frozenset((u, v))] for _, u, v, _ in checks}) == 1, colour
        assert len({a for *_, a in checks}) == len(checks), colour
    q2i = {q: i for i, q in enumerate(qubit_kinds(distance))}
    n_data = len(qubits(distance))
    for r in range(6):
        circuit = round_circuit(distance, r)
        checks = syndrome_checks(distance, rotated(r))
        assert {colours[frozenset((u, v))] for _, u, v, _ in checks} == {r % 3}, r  # the schedule's colour
        # the diagram labels agree with the circuit: MZZ reads S<k>, A<k> of edge k; every data qubit once
        _, label = roles(distance, r, True)
        mzz = [t.value for inst in circuit if inst.name == "MZZ" for t in inst.targets_copy()]
        assert [label[i] for i in mzz] == [f"{x}{k}" for k, *_ in checks for x in "SA"], r
        assert sorted(lab for lab in label if lab[0] == "D") == sorted(f"D{i}" for i in range(n_data)), r
        for j, (_, u, v, _) in enumerate(checks):
            assert frozenset((u, v)) in colours, (r, j)  # neighbouring data
            parity = stim.PauliString(circuit.num_qubits)
            parity[q2i[u]] = parity[q2i[v]] = "XZ"[r % 2]
            assert circuit.has_flow(stim.Flow(input=parity, measurements=[j])), (r, j)  # record j reads it
        for d in qubits(distance):  # one qubit at a time: MZZ copies a pair's parity onto its ancillas too
            single = stim.PauliString(circuit.num_qubits)
            single[q2i[d]] = "XZ"[r % 2]
            assert circuit.has_flow(stim.Flow(input=single, output=single)), (r, d)  # back in place


def layout_figure(distance: int, hex_view: bool = False, numbers: bool = False) -> plt.Figure:
    """Plaquettes in their colour, data qubits on the corners, one ancilla on every edge; `numbers` labels
    each qubit with its stim index (D<i> data qubit i, A<k> ancilla of edge k). One torus period is shown."""
    kinds = qubit_kinds(distance)
    pos = positions(distance, hex_view)
    p = period(distance)
    fig, ax = plt.subplots(figsize=(6, 8))
    for c, colour in hex_centers(distance).items():
        for shift in [a * p.real + b * p.imag * 1j for a in (-1, 0, 1) for b in (-1, 0, 1)]:
            corners = [c + shift + d for d in (HEX_CORNERS if hex_view else QUBIT_CORNERS)]
            ax.fill([q.real for q in corners], [q.imag for q in corners], color=COLORS[colour], alpha=0.25,
                    edgecolor="0.4")
    for kind, label in (("data", "data"), ("main", "ancilla")):
        qs = [pos[q] for q, k in kinds.items() if k == kind]
        ax.scatter([q.real for q in qs], [q.imag for q in qs], s=60 if kind == "data" else 36,
                   facecolor=QUBIT_KINDS[kind][0], edgecolor="#333333", linewidths=1.2,
                   label=f"{label} ({len(qs)})", zorder=3)
    n_data = list(kinds.values()).count("data")
    if numbers:
        for i, (q, kind) in enumerate(kinds.items()):
            ax.annotate(f"D{i}" if kind == "data" else f"A{i - n_data}", (pos[q].real, pos[q].imag),
                        xytext=(4, 3), textcoords="offset points", fontsize=6.5, color="#222222",
                        fontweight="bold" if kind == "data" else "normal", zorder=4)
        ax.text(1.04, 0.5, f"Labels\nD$i$  data qubit $i$, stim index $i$\n"
                           f"A$k$  ancilla of edge $k$, stim index {n_data} + $k$",
                transform=ax.transAxes, va="top", fontsize=9, color="#333333", linespacing=1.4)
    corner = window(list(pos.values()), p)
    ax.set_xlim(corner.real, corner.real + p.real)
    ax.set_ylim(corner.imag + p.imag, corner.imag)  # y grows downward, as in stim's timeslice diagrams
    ax.set_aspect("equal")
    ax.set_title(f"{len(kinds)} qubits, one ancilla per edge, " + ("regular hexagons (drawing only)" if hex_view
                                                                   else "brick wall (the circuits' coordinates)"),
                 fontsize=10)
    handles, _ = ax.get_legend_handles_labels()
    handles += [Patch(facecolor=COLORS[c], alpha=0.25, edgecolor="0.4",
                      label=f"colour {c}: steps {c} and {c + 3} check\nthe edges linking two of these")
                for c in range(3)]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False, fontsize=9)
    return fig


def layout_png(distance: int, hex_view: bool, numbers: bool, png: Path) -> None:
    fig = layout_figure(distance, hex_view, numbers)
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def round_svg(distance: int, r: int, hex_view: bool, numbers: bool = False) -> str:
    """Timeslice view of step r, plaquettes coloured underneath; brick-wall or regular hexagons. Tick 0 is the
    swaps; from tick 1 up to the swaps back, qubits are drawn in their roles after the swaps (see `roles`),
    labelled with `numbers`."""
    circuit = round_circuit(distance, r, hex_view)
    i2pos = {i: complex(*xy) for i, xy in circuit.get_final_qubit_coordinates().items()}
    svg = color_hexes(str(circuit.diagram("timeslice-svg")), hex_centers(distance),
                      HEX_CORNERS if hex_view else QUBIT_CORNERS, i2pos, period(distance))
    ticks = [roles(distance, r, False)] + [roles(distance, r, True)] * (circuit.num_ticks - 1)
    return style_ticks(svg, [k for k, _ in ticks], [lab for _, lab in ticks] if numbers else None)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--distance", type=int, default=2)
    distance = parser.parse_args().distance
    check(distance)
    out = Path(__file__).parent / f"one_ancilla_d{distance}"
    out.mkdir(exist_ok=True)
    for hex_view in (False, True):
        view = "_hex" if hex_view else ""
        layout_png(distance, hex_view, False, out / f"layout{view}.png")
        layout_png(distance, hex_view, True, out / f"layout{view}_numbered.png")
        for r in range(6):
            svg_png(round_svg(distance, r, hex_view), out / f"round{r}{view}.png")
            svg_png(round_svg(distance, r, hex_view, True), out / f"round{r}{view}_numbered.png")
    print(f"wrote {out}/")


if __name__ == "__main__":
    main()
