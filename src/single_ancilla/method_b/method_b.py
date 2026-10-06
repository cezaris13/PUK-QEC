"""Method B (hexes.pdf): each check's own ancilla collects the parity from both data next to it, then the
readout pair gathers just right of the corner u, whose charge sensor sits just left of u. How depends on u's
corner of the checking hex; on the red hex D0 D12 D13 D14 D2 D1 of the d = 2 layout (slides: D1 ... D6
clockwise from the top left):

    top left,    D1-D0 (slides D6-D1):  syndrome A13 (A6) swaps into D0, read against A0 (A1)
    right,       D12-D13 (D2-D3):       syndrome A5 (A2) swaps into D13, read against A15 (A9), sticking out
    bottom left, D14-D2 (D4-D5):        reference A2 (A5) swaps into D2, read against syndrome A12 (A4)

The first two swap the syndrome into u (`_syndrome_to_corner`), the last is method C (src/single_ancilla/method_c). The
slides also leave D4 and A3 swapped at the end; that looks like a slip, so here D14 and A28 stay where they are.

    .venv/bin/python src/single_ancilla/method_b/method_b.py --distance 2
"""

import method_c
import corner_readout
from floquet import collect


def layers(pauli: str, u: int, v: int, a: int, s: int, o: int,
           corner: int) -> tuple[corner_readout.Layers, int, int]:
    """Method B's gates for one check, then (syndrome, reference) as they start: by u's corner (see
    corner_readout.checks), the syndrome into u against s (top left) or o (right), or method C (bottom left)."""
    if corner == 2:
        return method_c.layers(pauli, u, v, a, s)
    return _syndrome_to_corner(pauli, u, v, a, s if corner == 4 else o)


def round_circuit(distance: int, r: int, hex_view: bool = False, code: str = "css"):
    """Step r for floquet.memory_circuit, as corner_readout.round_circuit builds it with these gates."""
    return corner_readout.round_circuit(distance, r, layers, hex_view, code)


def _syndrome_to_corner(pauli: str, u: int, v: int, a: int, ref: int) -> tuple[corner_readout.Layers, int, int]:
    """The syndrome a collects P_u P_v, swaps into u's spot and is read out against `ref` next to it, then u's
    data swaps back; (layers, syndrome, reference). Reset, CPs and H as in method_c.layers."""
    reset, first, second, turn = collect(pauli, u, v, a)
    return [[reset, ("R", [ref])], [first], [second], *turn, [("SWAP", [a, u])], [("MZZ", [u, ref])],
            [("SWAP", [a, u])]], a, ref


if __name__ == "__main__":
    corner_readout.main(layers, "method_b", __doc__)
