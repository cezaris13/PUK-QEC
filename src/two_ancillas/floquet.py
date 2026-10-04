
from dataclasses import dataclass
from typing import Callable, Iterable

import stim

# Hexagon corners relative to its centre: in qubit (brick-wall) coordinates, and on a regular hexagon
# for drawing. Same order, so zip(QUBIT_CORNERS, HEX_CORNERS) maps one onto the other.
QUBIT_CORNERS = [1, 1 + 1j, -1 + 1j, -1, -1 - 1j, 1 - 1j]
HEX_CORNERS = [4 / 3, 2 / 3 + 1j, -2 / 3 + 1j, -4 / 3, -2 / 3 - 1j, 2 / 3 - 1j]


@dataclass
class EdgeType:
    pauli: str
    hex_to_hex_delta: complex
    hex_to_qubit_delta: complex


EDGE_TYPES = [
    EdgeType(pauli="Z", hex_to_hex_delta=4, hex_to_qubit_delta=1),
    EdgeType(pauli="X", hex_to_hex_delta=2 - 3j, hex_to_qubit_delta=1 - 1j),
    EdgeType(pauli="Y", hex_to_hex_delta=2 + 3j, hex_to_qubit_delta=1 + 1j),
]


def period(distance: int) -> complex:
    """Size of the torus: coordinates wrap at this width (real part) and height (imaginary part)."""
    return complex(4 * distance, 6 * distance)


def torus(c: complex, distance: int) -> complex:
    return complex(c.real % period(distance).real, c.imag % period(distance).imag)


def sorted_complex(xs: Iterable[complex]) -> list[complex]:
    return sorted(xs, key=lambda v: (v.real, v.imag))


def hex_centers(distance: int) -> dict[complex, int]:
    """Plaquette centre -> colour (0, 1 or 2)."""
    centers = {}
    for row in range(3 * distance):
        for col in range(2 * distance):
            centers[torus(row * 2j + 2 * col - 1j * (col % 2), distance)] = (-row - col % 2) % 3
    return centers


def qubits(distance: int) -> list[complex]:
    """Data qubit coordinates, sorted; data qubit i is index i in every circuit, ancillas come after."""
    return sorted_complex({torus(h + e.hex_to_qubit_delta * sign, distance)
                           for h in hex_centers(distance) for e in EDGE_TYPES for sign in (-1, 1)})


@dataclass
class Edge:
    pauli: str
    colour: int  # the sub-round that measures it: the colour of the two plaquettes it joins
    data: tuple[complex, complex]
    # The two spin ancillas between the data qubits, read out together by spin blockade (as in SpinQEC):
    main: complex  # does the calculation (holds the syndrome); sits next to data[0]
    reff: complex  # its readout reference, held in a known state; sits next to data[1]


def edge_ends(h: complex, e: EdgeType, corners: list[complex]) -> tuple[complex, complex]:
    """The data qubits at the two ends of edge `e` leaving plaquette `h`, unwrapped, with each plaquette
    drawn using `corners` (QUBIT_CORNERS for brick-wall coordinates, HEX_CORNERS for regular hexagons)."""
    return (h + corners[QUBIT_CORNERS.index(e.hex_to_qubit_delta)],
            h + e.hex_to_hex_delta + corners[QUBIT_CORNERS.index(-e.hex_to_qubit_delta)])


def thirds(a: complex, b: complex, distance: int) -> tuple[complex, complex]:
    """The two points splitting a -> b into thirds: where an edge's ancillas sit."""
    return torus(a + (b - a) / 3, distance), torus(a + (b - a) * 2 / 3, distance)


def edge_list(distance: int) -> list[Edge]:
    """Every edge once, each with the two spin ancillas that sit between its data qubits."""
    out = []
    for h, colour in hex_centers(distance).items():
        for e in EDGE_TYPES:
            a, b = edge_ends(h, e, QUBIT_CORNERS)
            out.append(Edge(e.pauli, colour, (torus(a, distance), torus(b, distance)), *thirds(a, b, distance)))
    return out


def qubit_kinds(distance: int) -> dict[complex, str]:
    """Every qubit -> "data", "main" or "reff", in index order: a qubit's index in every circuit is its
    position here. Data qubits first (as in `qubits`), then each edge's main and reff."""
    kinds = {q: "data" for q in qubits(distance)}
    for edge in edge_list(distance):
        kinds[edge.main] = "main"
        kinds[edge.reff] = "reff"
    return kinds


def hex_positions(distance: int) -> dict[complex, complex]:
    """Qubit -> where it is drawn when plaquettes are regular hexagons rather than brick-wall rectangles."""
    positions = {}
    for h in hex_centers(distance):
        for qubit_corner, hex_corner in zip(QUBIT_CORNERS, HEX_CORNERS):
            p = torus(h + hex_corner, distance)
            existing = positions.setdefault(torus(h + qubit_corner, distance), p)
            assert abs(existing - p) < 1e-9, "the 3 hexes sharing a qubit must agree on where it goes"
        for e in EDGE_TYPES:  # ancillas stay a third of the way along their edge
            positions.update(zip(thirds(*edge_ends(h, e, QUBIT_CORNERS), distance),
                                 thirds(*edge_ends(h, e, HEX_CORNERS), distance)))
    return positions


def step_colour(r: int) -> int:
    """The colour step r measures. Colours 0, 1, 2 are b, r, g, so the schedule rX gZ bX rZ gX bZ is
    colour 1, 2, 0, 1, 2, 0 with Pauli X, Z alternating."""
    return (r + 1) % 3


def round_circuit(distance: int, r: int, hex_view: bool = False) -> stim.Circuit:
    """Step r of the 6-step schedule rX gZ bX rZ gX bZ: every edge of colour `step_colour(r)` has its XX (even r) or
    ZZ (odd r) parity measured by CNOT, M, CNOT on the pair, no ancillas.

    QUBIT_COORDS are the brick-wall coordinates, or the regular-hexagon ones with `hex_view`; they
    only move where timeslice diagrams draw each qubit.
    """
    edges_color = "rgb"  # colour 0, 1, 2 as plots.py draws them: blue, red (orange), green
    parity_meas = "XZ"  # even steps measure XX, odd steps ZZ


    # Qubit coordinate -> its index in the circuit.
    q2i = {}
    for index, q in enumerate(qubit_kinds(distance)):
        q2i[q] = index

    # Qubit coordinate -> where diagrams draw it.
    if hex_view:
        q2pos = hex_positions(distance)
    else:
        q2pos = {}
        for q in q2i:
            q2pos[q] = q

    # The two data qubits of every edge, grouped by colour. An edge leaves plaquette h at one corner and
    # reaches plaquette `far` at the opposite corner; those two plaquettes share a colour, the edge's colour.
    centers = hex_centers(distance)
    edges: dict[str, list[list[complex]]] = {"r": [], "g": [], "b": []}
    for h, colour in centers.items():
        for e in EDGE_TYPES:
            far = torus(h + e.hex_to_hex_delta, distance)
            assert centers[far] == colour, "an edge must connect two plaquettes of the same colour"
            d1 = torus(h + e.hex_to_qubit_delta, distance)  # where the edge leaves plaquette h
            d2 = torus(far - e.hex_to_qubit_delta, distance)  # where it reaches plaquette far
            edges[edges_color[colour]].append(sorted_complex([d1, d2]))

    circuit = stim.Circuit()
    for q, index in q2i.items():
        circuit.append("QUBIT_COORDS", [index], [q2pos[q].real, q2pos[q].imag])



    checks = [parity_check(parity_meas[r % 2], q2i[e.data[0]], q2i[e.data[1]], q2i[e.main], q2i[e.reff])
              for e in edge_list(distance) if e.colour == step_colour(r)]
    for layer in zip(*checks):  # layer k of every edge at once, one TICK each (same Pauli, same layer count)
        for part in layer:
            for name, targets in part:
                circuit.append(name, targets)
        circuit.append("TICK")
    return circuit


def parity_check(pauli: str, d0: int, d1: int, main: int, reff: int) -> list[list[tuple[str, list[int]]]]:
    """Measure P_d0 P_d1 (`pauli` "X", "Y" or "Z") on one edge through its spin ancilla pair.

    The reference always starts in |0> and the pair is read out as one Z parity (MZZ), like spin blockade:
    Z_syndrome * Z_reference, the reference adding its known +1. So only X or Y errors on the reference flip
    the outcome; it is immune to dephasing.
    X and Y: the syndrome (`main`, next to d0) starts in |+> and controls a P on d0, which kicks P_d0 back
    onto its X. SWAP moves it into the dot next to d1, where it controls a P on d1 and picks up P_d1. A final
    H on that dot (now `reff`) turns its X into the Z that MZZ reads.
    Z: no Hadamards. The syndrome starts in |0> too, and the CXs point the other way: d0, then d1 (after the
    SWAP) controls an X on the syndrome, copying Z_d0 Z_d1 straight into its Z.
    Returns its layers, each a list of (gate, targets), for the caller to TICK between: `round_circuit`
    runs layer k of every edge together. Leaves one measurement record.
    """
    if pauli == "Z":
        reset, first, second, turn = [("R", [main, reff])], ("CX", [d0, main]), ("CX", [d1, reff]), []
    else:
        reset, first, second, turn = ([("RX", [main]), ("R", [reff])], (f"C{pauli}", [main, d0]),
                                      (f"C{pauli}", [reff, d1]), [[("H", [reff])]])
    return [reset, [first], [("SWAP", [main, reff])], [second], *turn, [("MZZ", [main, reff])]]


def memory_circuit(distance: int, rounds: int,
                   step: Callable[[int, int], stim.Circuit] = round_circuit, basis: str = "Z") -> stim.Circuit:
    """Noiseless memory experiment of the CSS honeycomb code in `basis` (Z or X), with DETECTORs and one
    OBSERVABLE.

    Reset every data qubit to |0>, run `rounds` full periods of the 6-step schedule (see `round_circuit`),
    then measure every data qubit in Z. In the X basis the same with X and Z swapped: reset to |+>, start the
    schedule three steps on (step s + 3 checks the same colour as step s, in the other Pauli), read out in X.
    A Z memory never sees Z errors, nor an X memory X errors, so thresholds take the worse of the two.

    `step(distance, r)` builds step r: any circuit that puts the parity of each colour `step_colour(r)` edge into one
    record, in `edge_list` order, keeps data qubit i at stim index i, and sets its own QUBIT_COORDS.
    `round_circuit` by default; the src/single_ancilla/ schemes pass their own.

    The physics in brief, which the comments below refer back to:

    * Edges and plaquettes are 3-coloured. An edge of colour c joins two colour-c plaquettes, so it lies
      on the border of plaquettes of the other two colours. Every plaquette's border has 6 edges,
      alternating between those other two colours, 3 of each.
    * So one step (all edges of one colour, one Pauli) measures, for free, that Pauli's 6-body plaquette
      operator on every plaquette of the other two colours: X_P = product of the 3 XX checks on P's border.
    * A step measuring colour-c edges in Pauli X touches each colour-c plaquette at single corners only
      (its edges stick out of those plaquettes), so it anticommutes with, and randomises, their Z
      plaquette values. The same goes with X and Z swapped. Plaquettes of the other colours survive.
    * Detector = a plaquette value measured twice with nothing randomising it in between: without errors
      the two agree, so their XOR is 0.
    """
    data = qubits(distance)  # data qubit i is stim index i
    colour_of = hex_centers(distance)  # plaquette centre -> colour
    corners = {h: {torus(h + c, distance) for c in QUBIT_CORNERS} for h in colour_of}  # its 6 data qubits

    circuit = stim.Circuit()
    other = "ZX"[basis == "Z"]  # the Pauli whose plaquettes start unknown
    offset = 3 if basis == "X" else 0  # step s of an X memory is the Z memory's step s + 3
    circuit.append("R" if basis == "Z" else "RX", range(len(data)))
    circuit.append("TICK")

    # Measurements are numbered 0, 1, 2, ... in the order they happen; `count` is how many so far.
    # stim refers to earlier measurements relative to the latest one: measurement m is rec[m - count].
    count = 0

    def detector(measurements: list[int]) -> None:
        circuit.append("DETECTOR", [stim.target_rec(m - count) for m in measurements])

    # What we last learned about each plaquette operator: (Pauli, plaquette centre) -> the measurements whose
    # XOR is its value, or None if it is unknown (never measured, or randomised since).
    # After the reset every data qubit is Z = +1, so every Z plaquette is known to be +1: an empty list
    # (XOR of nothing = 0). X plaquettes are random after a Z reset: unknown.
    known = {(basis, h): [] for h in colour_of} | {(other, h): None for h in colour_of}

    # The logical observable: Z along a horizontal loop through data rows y = 0 and y = 1. No single fixed
    # operator survives every step, so its value hops between row 0, row 1 and both rows as the checks run.
    # Following it through the schedule, it works out to: XOR of every ZZ check with both data qubits in
    # rows 0-1, plus the final Z readout of row 1. (Checked numerically: deterministic without noise.)
    observable: list[int] = []
    in_band = lambda pair: pair[0].imag in (0, 1) and pair[1].imag in (0, 1)

    for s in range(6 * rounds):
        # Schedule rX gZ bX rZ gX bZ: Pauli alternates X, Z; colour cycles 1, 2, 0.
        pauli, other_pauli, colour = "XZ"[(s + offset) % 2], "ZX"[(s + offset) % 2], step_colour(s + offset)
        circuit += step(distance, s + offset)

        # One measurement per edge of this colour, in `edge_list` order (as `round_circuit` makes them).
        # edge's two data qubits -> the number of its measurement.
        measured = {e.data: count + i
                    for i, e in enumerate(e for e in edge_list(distance) if e.colour == colour)}
        count += len(measured)

        for h, c in colour_of.items():
            if c == colour:
                # This step's edges stick out of plaquette h: the other Pauli's value on h is now random.
                known[(other_pauli, h)] = None
                continue
            # This step's edges run along h's border: the 3 of them multiply to h's `pauli` plaquette.
            now = [m for pair, m in measured.items() if set(pair) <= corners[h]]
            if known[(pauli, h)] is not None:
                detector(now + known[(pauli, h)])  # new value XOR old value should be 0
            known[(pauli, h)] = now

        if pauli == basis:
            observable += [m for pair, m in measured.items() if in_band(pair)]

    # Final readout: every data qubit in Z. The last step (step 5 of the period) was ZZ on colour-0 edges.
    circuit.append("M" if basis == "Z" else "MX", range(len(data)))
    readout = {q: count + i for i, q in enumerate(data)}
    count += len(data)

    # Each last-step ZZ check vs the product of its two data qubits' readouts.
    for (d0, d1), m in measured.items():
        detector([m, readout[d0], readout[d1]])
    # Z plaquettes of the other two colours are products of those checks, so the line above already
    # covers them. Colour-0 Z plaquettes aren't: compare their last value with their 6 corners' readouts.
    for h, c in colour_of.items():
        if c == colour:
            detector(known[(basis, h)] + [readout[q] for q in corners[h]])

    observable += [readout[q] for q in data if q.imag == 1]
    # ponytail: one of the torus's two logical qubits; add the vertical loop as observable 1 to catch both.
    circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(m - count) for m in observable], 0)
    return circuit
