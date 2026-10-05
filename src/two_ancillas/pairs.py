"""Two ancillas per edge: a syndrome and a reference spin ancilla between every pair of data qubits, read out
together by spin blockade (as in SpinQEC). The step circuit for floquet.memory_circuit (src/shared/floquet.py).
"""
import sys
from pathlib import Path

import stim

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
from floquet import Edge, edge_list, hex_positions, qubits, step_colour, torus


def thirds(a: complex, b: complex, distance: int) -> tuple[complex, complex]:
    """The two points splitting a -> b into thirds: where an edge's ancillas sit."""
    return torus(a + (b - a) / 3, distance), torus(a + (b - a) * 2 / 3, distance)


def edge_ancillas(distance: int) -> list[tuple[Edge, complex, complex]]:
    """(edge, main, reff) for every edge in `edge_list` order: the two spin ancillas between its data qubits,
    read out together by spin blockade (as in SpinQEC). main does the calculation (holds the syndrome) and sits
    next to data[0]; reff, its readout reference held in a known state, sits next to data[1]."""
    return [(e, *thirds(*e.ends, distance)) for e in edge_list(distance)]


def positions(distance: int) -> dict[complex, complex]:
    """Qubit -> where it is drawn on regular hexagons: the data as `hex_positions`, the ancillas still a third
    of the way along their edge."""
    pos = hex_positions(distance)
    for e in edge_list(distance):
        pos.update(zip(thirds(*e.ends, distance), thirds(*e.hex_ends, distance)))
    return pos


def qubit_kinds(distance: int) -> dict[complex, str]:
    """Every qubit -> "data", "main" or "reff", in index order: a qubit's index in every circuit is its
    position here. Data qubits first (as in `qubits`), then each edge's main and reff."""
    kinds = {q: "data" for q in qubits(distance)}
    for _, main, reff in edge_ancillas(distance):
        kinds[main] = "main"
        kinds[reff] = "reff"
    return kinds


def round_circuit(distance: int, r: int, hex_view: bool = False) -> stim.Circuit:
    """Step r of the 6-step schedule rX gZ bX rZ gX bZ: every edge of colour `step_colour(r)` has its XX (even r) or
    ZZ (odd r) parity measured through its spin ancilla pair (`parity_check`).

    QUBIT_COORDS are the brick-wall coordinates, or the regular-hexagon ones with `hex_view`; they
    only move where timeslice diagrams draw each qubit.
    """
    parity_meas = "XZ"  # even steps measure XX, odd steps ZZ


    # Qubit coordinate -> its index in the circuit.
    q2i = {}
    for index, q in enumerate(qubit_kinds(distance)):
        q2i[q] = index

    # Qubit coordinate -> where diagrams draw it.
    if hex_view:
        q2pos = positions(distance)
    else:
        q2pos = {}
        for q in q2i:
            q2pos[q] = q

    circuit = stim.Circuit()
    for q, index in q2i.items():
        circuit.append("QUBIT_COORDS", [index], [q2pos[q].real, q2pos[q].imag])



    checks = [parity_check(parity_meas[r % 2], q2i[e.data[0]], q2i[e.data[1]], q2i[main], q2i[reff])
              for e, main, reff in edge_ancillas(distance) if e.colour == step_colour(r)]
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
