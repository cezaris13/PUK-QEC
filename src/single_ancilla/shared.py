"""Shared by methods B and C (src/single_ancilla/method_b, src/single_ancilla/method_c): one ancilla per edge, and each check's
syndrome is read out at a corner of a hex instead of the data walking around it (method A, src/single_ancilla/method_a,
whose layout and drawing helpers this builds on).

Step r checks the colour c = r % 3 edges from the hexes of colour c - 1 (on the red hex D0 D12 D13 D14 D2 D1 of
the d = 2 layout: the green edges D0-D1, D12-D13, D2-D14). `checks` gives each one's qubits; a method is a
`scheme(pauli, u, v, a, s, o, corner) -> (layers, syndrome, reference)` that turns one check into gate layers,
and this module runs every check's layer k together, verifies the result (`check`) and draws it (`main`).
"""
import argparse
import re
import sys
from pathlib import Path

import stim

sys.path.insert(0, str(Path(__file__).resolve().parent / "method_a"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "two_ancillas"))
from floquet import HEX_CORNERS, QUBIT_CORNERS, hex_centers, period, qubits, torus
from method_a import edge_list, positions, qubit_kinds, style_ticks
from plots import color_hexes, svg_png

Layers = list[list[tuple[str, list[int]]]]
# Where every single-ancilla run writes: results/single_ancilla/, mirroring src/single_ancilla/.
RESULTS = Path(__file__).resolve().parents[2] / "results" / "single_ancilla"


def checks(distance: int, colour: int) -> list[tuple[int, complex, complex, complex, complex, complex, int]]:
    """(k, u, v, a, s, o, corner) for every colour-`colour` edge, in `edge_list` order (the order memory_circuit
    wants the records in): v-u runs clockwise around a hex of colour `colour` - 1, a is the edge's ancilla (edge
    k), s the ancilla of the edge after it, from u, and o the ancilla of u's third edge, pointing out of the hex.
    `corner` is u's place in QUBIT_CORNERS: 0 right, 2 bottom-left, 4 top-left; u is always some hex's right
    corner. Those hexes hold every colour-`colour` edge once, each s is on a third-colour edge of one of them and
    each o on a `colour` - 1 edge sticking out of one, so no ancilla does two jobs."""
    edges = edge_list(distance)
    at = {frozenset((d0, d1)): (k, anc) for k, (_, d0, d1, anc) in enumerate(edges)}
    out = []
    for h, c in hex_centers(distance).items():
        if c == (colour - 1) % 3:
            ring = [torus(h + corner, distance) for corner in QUBIT_CORNERS]  # clockwise as drawn, y down
            for i in range(6):
                v, u, w = ring[i], ring[(i + 1) % 6], ring[(i + 2) % 6]
                k, a = at[frozenset((v, u))]
                if edges[k][0] == colour:
                    o = next(anc for _, d0, d1, anc in edges if u in (d0, d1) and not {d0, d1} & {v, w})
                    out.append((k, u, v, a, at[frozenset((u, w))][1], o, (i + 1) % 6))
    return sorted(out, key=lambda check: check[0])


def gates(scheme, pauli: str, check: tuple, q2i: dict[complex, int]) -> tuple[Layers, int, int]:
    """`scheme`'s (layers, syndrome, reference) for one of `checks`, qubits as stim indices."""
    _, *qs, corner = check
    return scheme(pauli, *(q2i[q] for q in qs), corner)


def round_circuit(distance: int, r: int, scheme, hex_view: bool = False) -> stim.Circuit:
    """Step r: every colour r % 3 edge's XX (even r) or ZZ (odd r) parity into one record, in `edge_list` order,
    `scheme`'s layer k of every check together, a TICK after each."""
    pos = positions(distance, hex_view)
    q2i = {q: i for i, q in enumerate(qubit_kinds(distance))}
    circuit = stim.Circuit()
    for q, i in q2i.items():
        circuit.append("QUBIT_COORDS", [i], [pos[q].real, pos[q].imag])
    for layer in zip(*[gates(scheme, "XZ"[r % 2], c, q2i)[0] for c in checks(distance, r % 3)]):
        for part in layer:
            for name, targets in part:
                circuit.append(name, targets)
        circuit.append("TICK")
    return circuit


def roles(distance: int, r: int, scheme) -> list[tuple[list[str], list[str]]]:
    """(kinds, labels) by stim index for each tick, as they stand before its gates: D<i> data qubit i, A<k> the
    ancilla of edge k, syndromes drawn as `main` and references as `reff`; each SWAP moves the pair's states,
    kinds and labels together."""
    kinds = qubit_kinds(distance)
    n_data = len(qubits(distance))
    kind = list(kinds.values())
    label = [f"D{i}" if k == "data" else f"A{i - n_data}" for i, k in enumerate(kind)]
    q2i = {q: i for i, q in enumerate(kinds)}
    for c in checks(distance, r % 3):
        kind[gates(scheme, "XZ"[r % 2], c, q2i)[2]] = "reff"
    circuit = round_circuit(distance, r, scheme)
    out = [(kind[:], label[:])]
    for inst in circuit:
        if inst.name == "SWAP":
            t = [x.value for x in inst.targets_copy()]
            for i, j in zip(t[::2], t[1::2]):
                kind[i], kind[j], label[i], label[j] = kind[j], kind[i], label[j], label[i]
        elif inst.name == "TICK":
            out.append((kind[:], label[:]))
    return out[:circuit.num_ticks]


def check(distance: int, scheme) -> None:
    """Each colour's checks cover every data qubit once with that colour's edges and never share an ancilla;
    every two-qubit gate joins a qubit spot to the ancilla spot of an edge next to it; each step reads its j-th
    check's parity into record j and puts every data qubit back; and the memory experiment's detectors are
    deterministic with these steps."""
    import floquet
    edges = edge_list(distance)
    colours = {frozenset((d0, d1)): c for c, d0, d1, _ in edges}
    q2i = {q: i for i, q in enumerate(qubit_kinds(distance))}
    next_to = {frozenset((q2i[d], q2i[anc])) for _, d0, d1, anc in edges for d in (d0, d1)}
    for colour in range(3):
        cs = checks(distance, colour)
        assert sorted(q2i[q] for _, u, v, *_ in cs for q in (u, v)) == list(range(len(qubits(distance)))), colour
        assert {colours[frozenset((u, v))] for _, u, v, *_ in cs} == {colour}, colour
        assert len({anc for _, _, _, a, s, o, _ in cs for anc in (a, s, o)}) == 3 * len(cs), colour
    for r in range(6):
        circuit = round_circuit(distance, r, scheme)
        for inst in circuit:
            if inst.name in ("SWAP", "CX", "CY", "CZ", "MZZ"):
                t = [x.value for x in inst.targets_copy()]
                assert all(frozenset(p) in next_to for p in zip(t[::2], t[1::2])), (r, inst)
        for j, (_, u, v, *_) in enumerate(checks(distance, r % 3)):
            parity = stim.PauliString(circuit.num_qubits)
            parity[q2i[u]] = parity[q2i[v]] = "XZ"[r % 2]
            assert circuit.has_flow(stim.Flow(input=parity, measurements=[j])), (r, j)
        for d in qubits(distance):
            single = stim.PauliString(circuit.num_qubits)
            single[q2i[d]] = "XZ"[r % 2]
            assert circuit.has_flow(stim.Flow(input=single, output=single)), (r, d)
    step = lambda d, r: round_circuit(d, r, scheme)
    floquet.memory_circuit(distance, 2, step).detector_error_model(allow_gauge_detectors=False)  # raises if not


def round_svg(distance: int, r: int, scheme, hex_view: bool, numbers: bool = False, away: bool = False) -> str:
    """Timeslice view of step r, plaquettes coloured underneath, qubits drawn in their roles at each tick, a
    square for each readout's sensor (see `readout_squares`)."""
    circuit = round_circuit(distance, r, scheme, hex_view)
    i2pos = {i: complex(*xy) for i, xy in circuit.get_final_qubit_coordinates().items()}
    svg = color_hexes(str(circuit.diagram("timeslice-svg")), hex_centers(distance),
                      HEX_CORNERS if hex_view else QUBIT_CORNERS, i2pos, period(distance))
    ticks = roles(distance, r, scheme)
    svg = style_ticks(svg, [k for k, _ in ticks], [lab for _, lab in ticks] if numbers else None)
    q2i = {q: i for i, q in enumerate(qubit_kinds(distance))}
    right = {q2i[torus(h + 1, distance)] for h in hex_centers(distance)}
    mzz = [t.value for inst in circuit if inst.name == "MZZ" for t in inst.targets_copy()]
    pairs = [(i, j) if i in right else (j, i) for i, j in zip(mzz[::2], mzz[1::2])]
    return readout_squares(local_pairs(svg, i2pos, period(distance)), pairs, pixel_period(svg, i2pos, period(distance)),
                           away)


def dots(svg: str) -> dict[int, list[complex]]:
    """Qubit -> its dot's pixel position in each tick's panel, in panel order."""
    out = {}
    for q, x, y in re.findall(r'id="qubit_dot:(\d+):[^"]*" cx="([-\d.]+)" cy="([-\d.]+)"', svg):
        out.setdefault(int(q), []).append(complex(float(x), float(y)))
    return out


def pixel_period(svg: str, i2pos: dict[int, complex], torus_period: complex) -> complex:
    """The torus period in pixels, from tick 0's dots, as plots.color_hexes finds its scale."""
    first = {q: ps[0] for q, ps in dots(svg).items()}
    lo, hi = min(first, key=lambda q: i2pos[q].real), max(first, key=lambda q: i2pos[q].real)
    top, bottom = min(first, key=lambda q: i2pos[q].imag), max(first, key=lambda q: i2pos[q].imag)
    return complex(torus_period.real * (first[hi].real - first[lo].real) / (i2pos[hi].real - i2pos[lo].real),
                   torus_period.imag * (first[bottom].imag - first[top].imag) / (i2pos[bottom].imag - i2pos[top].imag))


def segments(a: complex, b: complex, per: complex) -> list[tuple[complex, complex]]:
    """a to b as drawn: straight if they are neighbours inside the picture, else a stub from each end halfway to
    the other end's image across the torus boundary."""
    shift = complex(per.real * round((b - a).real / per.real), per.imag * round((b - a).imag / per.imag))
    return [(a, b)] if not shift else [(a, (a + b - shift) / 2), (b, (a + b + shift) / 2)]


def local_pairs(svg: str, i2pos: dict[int, complex], torus_period: complex) -> str:
    """A pair that wraps around the torus is neighbours, but stim joins its two dots straight across the
    picture. Redraw each such link as `segments`, so every SWAP, CP and MZZ reads as the local gate it is."""
    per = pixel_period(svg, i2pos, torus_period)

    def link(m: re.Match[str]) -> str:
        parts = segments(complex(float(m[1]), float(m[2])), complex(float(m[3]), float(m[4])), per)
        if len(parts) == 1:
            return m[0]
        d = "".join(f"M{p.real},{p.imag} L{q.real},{q.imag} " for p, q in parts)
        return f'<path d="{d}" fill="none" stroke="black" stroke-width="5"/>'
    return re.sub(r'<path d="M([-\d.]+),([-\d.]+) [LC](?:[-\d.]+ [-\d.]+,[-\d.]+ [-\d.]+,)?([-\d.]+)[, ]([-\d.]+) " '
                  r'fill="none" stroke="black" stroke-width="5"/>', link, svg)


def readout_squares(svg: str, pairs: list[tuple[int, int]], per: complex, away: bool = False) -> str:
    """A hollow square for the charge sensor of each spin-blockade readout (corner, partner) in `pairs`, stim
    indices, the corner always some hex's right corner, in every tick: just left of the corner, or with `away`
    on the far side of the corner from its partner."""
    at = dots(svg)
    squares = []
    for u, partner in pairs:
        for p, q in zip(at[u], at[partner]):
            end = segments(p, q, per)[0][1]  # towards the partner, or its image across the boundary
            c = p - 30 * ((end - p) / abs(end - p) if away else 1)
            squares.append(f'<rect x="{c.real - 8}" y="{c.imag - 8}" width="16" height="16" fill="none" '
                           f'stroke="#222222" stroke-width="2"/>')
    k = svg.index('<circle id="qubit_dot')
    return svg[:k] + "\n".join(squares) + "\n" + svg[k:]


def main(scheme, name: str, doc: str, away: bool = False) -> None:
    """`scheme`'s checks, then its 6 steps' timeslices, plain and numbered, brick-wall and hex view, into
    results/single_ancilla/<name>/d<distance>/; `doc` is the --help text."""
    parser = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--distance", type=int, default=2)
    distance = parser.parse_args().distance
    check(distance, scheme)
    out = RESULTS / name / f"d{distance}"
    out.mkdir(parents=True, exist_ok=True)
    for hex_view in (False, True):
        view = "_hex" if hex_view else ""
        for r in range(6):
            svg_png(round_svg(distance, r, scheme, hex_view, False, away), out / f"round{r}{view}.png")
            svg_png(round_svg(distance, r, scheme, hex_view, True, away), out / f"round{r}{view}_numbered.png")
    print(f"wrote {out}/")

