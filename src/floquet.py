
from dataclasses import dataclass
from typing import Iterable

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


def round_circuit(distance: int, r: int, hex_view: bool = False) -> stim.Circuit:
    """Sub-round r: every colour-r edge's parity measured by CNOT, M, CNOT on the pair, no ancillas.

    QUBIT_COORDS are the brick-wall coordinates, or the regular-hexagon ones with `hex_view`; they
    only move where timeslice diagrams draw each qubit.
    """
    q2i = {q: i for i, q in enumerate(qubit_kinds(distance))}
    q2pos = hex_positions(distance) if hex_view else {q: q for q in q2i}
    edges: dict[str, list[list[complex]]] = {"X": [], "Y": [], "Z": []}
    for h, colour in hex_centers(distance).items():
        if colour == r:
            for e in EDGE_TYPES:
                pair = [torus(h + e.hex_to_qubit_delta, distance),
                        torus(h + e.hex_to_hex_delta - e.hex_to_qubit_delta, distance)]
                edges[e.pauli].append(sorted_complex(pair))

    circuit = stim.Circuit()
    for q, i in q2i.items():
        circuit.append("QUBIT_COORDS", [i], [q2pos[q].real, q2pos[q].imag])
    x_qubits = [q2i[q] for pair in edges["X"] for q in pair]
    y_qubits = [q2i[q] for pair in edges["Y"] for q in pair]
    pairs = [q2i[q] for group in edges.values() for pair in group for q in pair]
    circuit.append("H", x_qubits)  # every parity becomes a ZZ parity...
    circuit.append("H_YZ", y_qubits)
    circuit.append("TICK")
    circuit.append("CNOT", pairs)  # ...and then a single Z on the pair's second qubit
    circuit.append("TICK")
    circuit.append("M", pairs[1::2])
    circuit.append("TICK")
    circuit.append("CNOT", pairs)
    circuit.append("TICK")
    circuit.append("H_YZ", y_qubits)  # restore the bases
    circuit.append("H", x_qubits)
    return circuit


def memory_circuit(distance: int, rounds: int) -> stim.Circuit:
    """Noiseless memory experiment with DETECTORs and OBSERVABLE_INCLUDEs, layers separated by TICKs."""
    return stim.Circuit.generated("surface_code:rotated_memory_z", distance=distance, rounds=rounds)
