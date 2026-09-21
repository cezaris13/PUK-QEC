"""The code under test: builds the noiseless memory-experiment circuit that `memory.py` adds noise to."""
import stim


def memory_circuit(distance: int, rounds: int) -> stim.Circuit:
    """Noiseless memory experiment with DETECTORs and OBSERVABLE_INCLUDEs, layers separated by TICKs."""
    return stim.Circuit.generated("surface_code:rotated_memory_z", distance=distance, rounds=rounds)
