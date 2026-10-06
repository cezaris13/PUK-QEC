"""Method C (hexes.pdf): each check's own ancilla collects the parity from both data next to it and stays put;
its reference swaps into the corner u next to it, and the pair is read out by the charge sensor on the far side
of u. Checks as in src/single_ancilla/shared/corner_readout.py: on the red hex D0 D12 D13 D14 D2 D1 of the d = 2
layout (slides: D1 ... D6 clockwise from the top left),

    D1-D0 (slides D6-D1): syndrome A13 (A6), reference A0 (A1) swaps into D0, sensor above D0
    D12-D13 (D2-D3):      syndrome A5 (A2),  reference A28 (A3) swaps into D13, sensor below right of D13
    D14-D2 (D4-D5):       syndrome A12 (A4), reference A2 (A5) swaps into D2, sensor left of D2

    CP   a -> u, a -> v   syndrome picks up P_u P_v
    SWAP s u              reference onto u's spot, u's data onto s's
    MZZ  a u              spin-blockade readout
    SWAP s u              u's data back home

    .venv/bin/python src/single_ancilla/method_c/method_c.py --distance 2
"""

import corner_readout
from floquet import collect


def layers(pauli: str, u: int, v: int, a: int, s: int, *_) -> tuple[corner_readout.Layers, int, int]:
    """Method C's gates for one check, then (syndrome, reference) as they start; `pauli` is two letters, P_u then
    P_v ("XX", "ZZ", or the X3Z3 code's "XZ" / "ZX"); the syndrome's reset, CPs and H as in method_a."""
    reset, first, second, turn = collect(pauli, u, v, a)
    return [[reset, ("R", [s])], [first], [second], *turn, [("SWAP", [s, u])], [("MZZ", [a, u])], [("SWAP", [s, u])]], a, s


def round_circuit(distance: int, r: int, hex_view: bool = False, code: str = "css"):
    """Step r for floquet.memory_circuit, as corner_readout.round_circuit builds it with these gates."""
    return corner_readout.round_circuit(distance, r, layers, hex_view, code)


if __name__ == "__main__":
    corner_readout.main(layers, "method_c", __doc__, away=True)
