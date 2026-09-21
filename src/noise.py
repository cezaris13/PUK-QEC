from dataclasses import dataclass, field
from functools import partial

import stim


@dataclass
class NoiseTerm:
    """Noise attached to one gate type: what happens to the gated qubits, and to the idle ones.

    Leave a name empty ("") to switch that channel off.
    `probs` is whatever the named stim instruction expects: one number for
    DEPOLARIZE1/X_ERROR, 3 for PAULI_CHANNEL_1, 15 for PAULI_CHANNEL_2.
    """

    gate_noise_name: str = ""
    gate_noise_probs: list[float] | None = None
    idle_noise_name: str = ""
    idle_noise_probs: list[float] | None = None


@dataclass
class NoiseModel:
    """All the noise terms that define one error model.

    `round_noise` is applied to the data qubits once per syndrome round, on top of
    the per-layer idle noise. Use it for the dead time while ancillas are being
    read out and reset -- the dominant idling cost for spin qubits, since spin
    readout (~microseconds) is far slower than a gate (~tens of nanoseconds).
    """

    R_gate: NoiseTerm = field(default_factory=NoiseTerm)
    M_gate: NoiseTerm = field(default_factory=NoiseTerm)
    H_gate: NoiseTerm = field(default_factory=NoiseTerm)
    CX_gate: NoiseTerm = field(default_factory=NoiseTerm)
    CZ_gate: NoiseTerm = field(default_factory=NoiseTerm)
    SWAP_gate: NoiseTerm = field(default_factory=NoiseTerm)
    round_noise_name: str = ""
    round_noise_probs: list[float] | None = None


def biased_single_qubit_probabilities(p: float, eta: float) -> list[float]:
    """Split a total single-qubit error rate `p` into [p_x, p_y, p_z] for PAULI_CHANNEL_1.

    Bias is defined as eta = p_z / (p_x + p_y), so eta = 0.5 is depolarizing and
    eta -> infinity is pure dephasing. Returns probabilities summing to `p`.
    """
    dephasing = p * eta / (eta + 1)
    other = p / (2 * (eta + 1))
    return [other, other, dephasing]


def biased_two_qubit_probabilities(p: float, eta: float) -> list[float]:
    """Split a total two-qubit error rate `p` into the 15 PAULI_CHANNEL_2 probabilities.

    Order is stim's: IX IY IZ XI XX XY XZ YI YX YY YZ ZI ZX ZY ZZ. The weight
    landing on the three Z-flavoured terms (IZ, ZI, ZZ) is

        zeta = 3/5 * (eta/(eta+1))**2 + 2/5 * (eta/(eta+1))

    which is what you get by applying the single-qubit bias independently to both
    halves of the pair. Returns probabilities summing to `p`.

    https://doi.org/10.1038/s41534-025-01074-1 page 11 formula 8
    """
    zeta = 3 / 5 * (eta / (eta + 1)) ** 2 + 2 / 5 * (eta / (eta + 1))
    dominant = p * zeta / 3          # IZ, ZI, ZZ
    other = p * (1 - zeta) / 12      # the remaining twelve
    probs = [other] * 15
    for i in (2, 11, 14):            # IZ, ZI, ZZ
        probs[i] = dominant
    return probs


def spin_noise_model(p: float, eta: float = 10.0) -> NoiseModel:
    """The one-parameter circuit-level model used for this repo's threshold sweeps.

    Everything scales off `p`, the two-qubit gate error rate, so a threshold is a
    single crossing point. In descending order:

        readout            5p      reset               2p
        idling per round   2p      SWAP              1.5p
        2-qubit gate        p      1-qubit gate      p/10
        idling per layer  p/10

    "Idling per round" is the dead time while the ancilla pairs are read out and
    reset -- long for spin qubits, hence its size relative to a gate.

    Provenance: these ratios are lifted verbatim from the sweep cell of the
    original ThressholdPlots.ipynb (commit 9e7c994). They are a modelling choice
    with no derivation in this repo, and they are *not* from Hetenyi & Wootton --
    that paper instead computes readout error from the integration time (Eqs. 6-8)
    and idling error from tau, T1 and T2 (Eqs. 3-5) rather than pinning them to
    the gate error rate. Two further departures from the paper worth knowing:

    * it applies depolarizing noise to gates, whereas this biases them by `eta`;
    * the implied single/two-qubit gate ratio here is eta_G = 10, while the
      paper's Table I is quoted at eta_G = 1.

    So numbers from this model are internally consistent but not directly
    comparable to the paper's tables. `eta` biases every channel towards
    dephasing; see the module docstring and `biased_two_qubit_probabilities`.
    """
    idle = dict(idle_noise_name="PAULI_CHANNEL_1",
                idle_noise_probs=biased_single_qubit_probabilities(p / 10, eta))
    two_qubit = NoiseTerm(gate_noise_name="PAULI_CHANNEL_2",
                          gate_noise_probs=biased_two_qubit_probabilities(p, eta), **idle)
    return NoiseModel(
        R_gate=NoiseTerm(gate_noise_name="X_ERROR", gate_noise_probs=[2 * p]),
        M_gate=NoiseTerm(gate_noise_name="X_ERROR", gate_noise_probs=[5 * p]),
        H_gate=NoiseTerm(gate_noise_name="PAULI_CHANNEL_1",
                         gate_noise_probs=biased_single_qubit_probabilities(p / 10, eta),
                         **idle),
        CX_gate=two_qubit,
        CZ_gate=two_qubit,
        SWAP_gate=NoiseTerm(gate_noise_name="PAULI_CHANNEL_2",
                            gate_noise_probs=biased_two_qubit_probabilities(1.5 * p, eta),
                            **idle),
        round_noise_name="PAULI_CHANNEL_1",
        round_noise_probs=biased_single_qubit_probabilities(2 * p, eta),
    )


def hbd_two_qubit_probabilities(p: float, eta: float) -> list[float]:
    """The HBD paper's split of a two-qubit error rate `p` into 15 PAULI_CHANNEL_2 probs.

    Same shape as `biased_two_qubit_probabilities` -- weight concentrated on IZ, ZI,
    ZZ -- but a different amount of it. HBD puts

        zeta = eta / (eta + 1)

    on the three Z-flavoured terms, i.e. the single-qubit bias applied directly to
    the pair, where `biased_two_qubit_probabilities` applies it independently to
    each half and gets 3/5*x^2 + 2/5*x. HBD is the more biased of the two at every
    finite eta (at eta = 10: 0.909 vs 0.859). Returns probabilities summing to `p`.

    https://journals.aps.org/prapplied/abstract/10.1103/q7w6-nljp
    """
    dominant = p * eta / (3 * (1 + eta))    # IZ, ZI, ZZ
    other = p / (12 * (1 + eta))            # the remaining twelve
    probs = [other] * 15
    for i in (2, 11, 14):
        probs[i] = dominant
    return probs


def hbd_noise_model(p: float, eta: float = 10.0, si1000: bool = False) -> NoiseModel:
    """The HBD circuit-level model, as published (`si1000=False`) or SI1000-scaled.

    Everything scales off `p`, the two-qubit gate error rate:

                        HBD     MHBD (si1000=True)
        2-qubit gate      p       p
        1-qubit gate      p       p/10
        reset             p       2p
        readout           p       5p
        idling per round  p       2p

    HBD is flat -- every operation costs `p`. MHBD re-weights it the way SI1000
    (https://quantum-journal.org/papers/q-2021-12-20-605/) re-weights its own
    model: cheap single-qubit gates, expensive readout, idling between the two.
    That makes MHBD, not HBD, the one to compare against `spin_noise_model`, whose
    ratios are the same except for reset (2p in both) and how the gates are biased.

    Two differences from `spin_noise_model` survive in both variants:

    * single-qubit gates are *depolarizing*, not biased -- straight from the paper;
    * there is no per-layer idle noise, only the per-round term. Idle qubits pay
      nothing while a gate layer runs.

    https://journals.aps.org/prapplied/abstract/10.1103/q7w6-nljp
    """
    one_q, reset, readout, rounds = (0.1, 2, 5, 2) if si1000 else (1, 1, 1, 1)
    two_qubit = NoiseTerm(gate_noise_name="PAULI_CHANNEL_2",
                          gate_noise_probs=hbd_two_qubit_probabilities(p, eta))
    return NoiseModel(
        R_gate=NoiseTerm(gate_noise_name="X_ERROR", gate_noise_probs=[reset * p]),
        M_gate=NoiseTerm(gate_noise_name="X_ERROR", gate_noise_probs=[readout * p]),
        H_gate=NoiseTerm(gate_noise_name="DEPOLARIZE1", gate_noise_probs=[one_q * p]),
        CX_gate=two_qubit,
        CZ_gate=two_qubit,
        SWAP_gate=two_qubit,
        round_noise_name="PAULI_CHANNEL_1",
        round_noise_probs=biased_single_qubit_probabilities(rounds * p, eta),
    )


#: name -> f(p, eta) -> NoiseModel. `--noise` on the sweep picks from these.
NOISE_MODELS = {
    "spin": spin_noise_model,
    "HBD": hbd_noise_model,
    "MHBD": partial(hbd_noise_model, si1000=True),
}


def _flip_after_reset(circuit: stim.Circuit, gate_name: str, qubits: list[int], term: NoiseTerm) -> None:
    """Append `term`'s reset error: a flip that actually flips the state `gate_name` prepared."""
    if term.gate_noise_name:
        circuit.append("Z_ERROR" if gate_name.endswith("X") else "X_ERROR", qubits, term.gate_noise_probs[0])


def add_noise(circuit: stim.Circuit, model: NoiseModel, qubits: list[int] | None = None) -> stim.Circuit:
    """A copy of the noiseless `circuit` with `model` applied one layer (TICK to TICK) at a time.

    Two-qubit gates get CZ_gate, SWAP_gate or (any other) CX_gate noise; single-qubit unitaries get
    H_gate. Resets get R_gate right after. Measurements, MPP included, get M_gate as a flipped
    result -- readout misassignment, which leaves the qubit alone. Qubits no gate touches in a layer
    get that layer's idle noise, and in a layer that measures they also get `round_noise`: the wait
    while the ancillas are read out.
    """
    if qubits is None:
        qubits = sorted({t.qubit_value for inst in circuit.flattened() for t in inst.targets_copy()
                         if t.qubit_value is not None})
    out = stim.Circuit()
    touched: set[int] = set()
    terms: list[NoiseTerm] = []
    measured = False

    def end_layer() -> None:
        nonlocal measured
        idle = [q for q in qubits if q not in touched]
        term = next((t for t in terms if t.idle_noise_name), None)
        if idle and term:
            out.append(term.idle_noise_name, idle, term.idle_noise_probs)
        if idle and measured and model.round_noise_name:
            out.append(model.round_noise_name, idle, model.round_noise_probs)
        touched.clear()
        terms.clear()
        measured = False

    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            end_layer()
            out.append(stim.CircuitRepeatBlock(inst.repeat_count, add_noise(inst.body_copy(), model, qubits)))
            continue
        if inst.name == "TICK":
            end_layer()
            out.append(inst)
            continue
        gate = stim.gate_data(inst.name)
        targets = inst.targets_copy()
        gate_qubits = [t.qubit_value for t in targets if t.qubit_value is not None]
        if gate.produces_measurements:
            m = model.M_gate
            out.append(inst.name, targets, m.gate_noise_probs[0]) if m.gate_noise_name else out.append(inst)
            if gate.is_reset:
                _flip_after_reset(out, gate.name, gate_qubits, model.R_gate)
            terms.append(m)
            measured = True
        elif gate.is_reset:
            out.append(inst)
            _flip_after_reset(out, gate.name, gate_qubits, model.R_gate)
            terms.append(model.R_gate)
        elif gate.is_unitary:
            term = ({"CZ": model.CZ_gate, "SWAP": model.SWAP_gate}.get(gate.name, model.CX_gate)
                    if gate.is_two_qubit_gate else model.H_gate)
            out.append(inst)
            if term.gate_noise_name:
                out.append(term.gate_noise_name, targets, term.gate_noise_probs)
            terms.append(term)
        else:  # annotations (DETECTOR, QUBIT_COORDS, ...) pass straight through
            out.append(inst)
            continue
        touched.update(gate_qubits)
    end_layer()
    return out
