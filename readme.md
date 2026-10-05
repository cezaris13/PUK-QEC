# Floquet honeycomb code with spin-qubit readout

Memory experiments for the Floquet honeycomb code, with readout through spin ancillas. The code compares
several readout schemes:

- **pairs**: a syndrome and a reference ancilla on every edge (`src/two_ancillas/`)
- **method_a / method_b / method_c**: one ancilla per edge, read out at a corner of a hex, following
  hexes.pdf (`src/single_ancilla/method_*/`)

Every scheme goes through the same `floquet.memory_circuit`, noise model and decoder (Stim + PyMatching), so the
only thing that differs between them is the step circuits.

Step 0 of each scheme at distance 1, one panel per tick (`make drawings`, hex view):

| pairs | method_a |
|---|---|
| <img src="docs/readme/pairs.png" width="400"> | <img src="docs/readme/method_a.png" width="400"> |
| **method_b** | **method_c** |
| <img src="docs/readme/method_b.png" width="400"> | <img src="docs/readme/method_c.png" width="400"> |

## Setup

```bash
make venv          # .venv with requirements.txt, plus a Jupyter kernel "floquet_code"
```

## Running

Every result is a `make` target. Settings can be overridden on the command line, and their defaults are
documented at the top of the [Makefile](Makefile).

| Target | What it makes |
|---|---|
| `make decode` | LER vs p for every scheme in `SCHEMES` (default: all four), in one CSV and one plot |
| `make plots` | LER vs p for the pairs scheme only, one CSV per distance |
| `make footprint` | LER vs physical qubit count at one `P`, extrapolated to `TARGET` |
| `make thresholds` | threshold surface and threshold vs bias for `SCHEME` (after IBM's QEC-with-spin-qubits) |
| `make drawings` | layouts and timeslices for every scheme at `DRAW_DISTANCE` |
| `make docs` | the LaTeX write-ups in `docs/` (needs tectonic) |
| `make serve` | interactive viewer of `results/` (needs node) |

Example:

```bash
make decode NOISE=MHBD ETA=100 DISTANCES="3 5 7" SHOTS=100000
```

Output file names include their settings, so runs with different settings don't overwrite each other.

## Layout

```
src/two_ancillas/     floquet.py (circuits), noise.py (spin / HBD / MHBD), memory.py, plots.py
src/single_ancilla/   method_a|b|c/, shared.py, decode.py, footprint.py, thresholds.py, detectors.py
notebooks/            walkthroughs of each scheme and of the memory experiment
results/              all generated output, laid out like src/
docs/                 write-ups: floquet_rounds, decoding, spin_x3z3_readout, one_ancilla_review
viewer/               Vite + React viewer behind `make serve`
```
