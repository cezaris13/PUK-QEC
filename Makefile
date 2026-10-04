PY = .venv/bin/python

# Sweep settings. Override any on the command line:
#   make plots NOISE=MHBD ETA=100 DISTANCES="3 5 7 9" SHOTS=100000
# Run time grows roughly as SHOTS x NUM x (number of distances) x d^3 (~d^2 qubits for d rounds),
# so the largest distance dominates. notebooks/two_ancillas/memory_experiment.ipynb shows each setting at work.

# Code distances. One CSV and one LER line each; rounds = distance. `make timeslices` draws the first.
DISTANCES = 3 5 7
# Noise model from NOISE_MODELS in src/two_ancillas/noise.py: spin, HBD or MHBD.
NOISE = spin
# Bias eta = p_z / (p_x + p_y): 0.5 is depolarizing, large is mostly dephasing. spin biases every
# channel; HBD/MHBD bias two-qubit gates and per-round idling, but keep 1-qubit gates depolarizing.
ETA = 10
# Monte Carlo shots per (distance, p) point. Smallest resolvable LER ~ 1/SHOTS; a point's relative
# error ~ 1/sqrt(logical errors), so 5 errors is +-45%.
SHOTS = 10000
# Physical error rate range. p is the two-qubit gate error; the model scales everything else off it
# (spin: readout 5p, reset 2p, idling per round 2p, 1-qubit gates p/10).
P_MIN = 1e-4
P_MAX = 1e-2
# Number of p values from P_MIN to P_MAX, both included, log-spaced: 9 over two decades = 4 per decade.
NUM = 9
# Honeycomb lattice size for `make honeycomb` drawings: 12*d^2 data qubits plus 2 spin ancillas on each of
# the 18*d^2 edges, so 1 is 12 + 36 qubits and 2 is 48 + 144; bigger gets hard to read.
DRAW_DISTANCE = 1
# Set NUMBERS=1 to label every qubit in the layout and timeslice drawings: Di data, Sk / Rk syndrome / reference ancilla.
NUMBERS = 1

# Readout schemes for `make decode` and `make footprint` (see src/single_ancilla/decode.py): pairs (two ancillas per edge),
# method_a, method_b, method_c (hexes.pdf). One file per set of schemes.
SCHEMES = pairs method_a method_b method_c
# `make footprint`: LER vs physical qubits at this one p, extrapolated to TARGET. Low p needs many shots per point
# before the larger distances see any logical errors to fit.
P = 1e-3
TARGET = 1e-6
FOOTPRINT_SHOTS = 10000000
# Distances for the fit: d = 5 already needs ~10^7 shots at p = 1e-3 to see errors, d = 7 far more.
FOOTPRINT_DISTANCES = 3 5

# All output goes under results/, mirroring src/: the two-ancilla code's runs, and the single-ancilla methods'
# (each method's round pictures in its own folder) with the decode and footprint comparisons.
TWO = results/two_ancillas
ONE = results/single_ancilla

empty :=
space := $(empty) $(empty)
TAG = $(subst $(space),-,$(strip $(SCHEMES)))
FOOTPRINT = $(ONE)/$(NOISE)_eta$(ETA)_p$(P)_d$(subst $(space),-,$(strip $(FOOTPRINT_DISTANCES)))_shots$(FOOTPRINT_SHOTS)_$(TAG)_footprint

RUN_NAME = $(NOISE)_eta$(ETA)_p$(P_MIN)-$(P_MAX)x$(NUM)_shots$(SHOTS)
RUN = $(TWO)/$(RUN_NAME)
DECODE = $(ONE)/$(RUN_NAME)_$(TAG)_decode
CSVS = $(foreach d,$(DISTANCES),$(RUN)_d$(d).csv)
TIMESLICES = $(TWO)/timeslices_d$(firstword $(DISTANCES)).png
HONEYCOMB = $(TWO)/honeycomb$(if $(NUMBERS),_numbered)_d$(DRAW_DISTANCE)/layout.png

venv:
	python3 -m venv .venv
	.venv/bin/pip install -U pip
	.venv/bin/pip install -r requirements.txt
	.venv/bin/python -m ipykernel install --user --name floquet_code

timeslices: $(TIMESLICES)

# Qubit layout, plus the 3 sub-rounds' coloured timeslices in brick-wall and regular-hexagon views.
honeycomb: $(HONEYCOMB)

plots: $(TIMESLICES) $(HONEYCOMB) $(RUN)_ler.png

# Methods A, B and C of hexes.pdf, one ancilla per edge (src/single_ancilla/method_*/): each one's 6 steps' timeslices, plain
# and numbered, at DRAW_DISTANCE, in results/single_ancilla/method_*/d<DRAW_DISTANCE>/; method A also draws the layouts.
experiments: $(ONE)/method_a/d$(DRAW_DISTANCE)/layout.png \
	$(ONE)/method_b/d$(DRAW_DISTANCE)/round0.png $(ONE)/method_c/d$(DRAW_DISTANCE)/round0.png

# One 6-step period only: every TICK is a tile, so d periods would be far too big to draw.
$(TWO)/timeslices_d%.png: src/two_ancillas/floquet.py src/two_ancillas/plots.py
	$(PY) src/two_ancillas/plots.py timeslices --distance $* --rounds 1 --png $@

# layout.png stands in for the whole folder: the round PNGs are written alongside it.
$(TWO)/honeycomb_d%/layout.png: src/two_ancillas/floquet.py src/two_ancillas/plots.py
	$(PY) src/two_ancillas/plots.py honeycomb --distance $* --out $(@D)

# The same drawings with numbered qubits, in their own folder so neither ever passes for the other.
$(TWO)/honeycomb_numbered_d%/layout.png: src/two_ancillas/floquet.py src/two_ancillas/plots.py
	$(PY) src/two_ancillas/plots.py honeycomb --distance $* --numbers --out $(@D)

# layout.png stands in for the folder here too; the script picks the folder from the distance.
$(ONE)/method_a/d%/layout.png: src/single_ancilla/method_a/method_a.py src/two_ancillas/floquet.py src/two_ancillas/plots.py
	$(PY) src/single_ancilla/method_a/method_a.py --distance $*

# Methods B and C, on src/single_ancilla/shared.py; round0.png stands in for each folder's timeslices.
$(ONE)/method_b/d%/round0.png: src/single_ancilla/method_b/method_b.py src/single_ancilla/method_c/method_c.py src/single_ancilla/shared.py \
		src/single_ancilla/method_a/method_a.py src/two_ancillas/floquet.py src/two_ancillas/plots.py
	$(PY) src/single_ancilla/method_b/method_b.py --distance $*

$(ONE)/method_c/d%/round0.png: src/single_ancilla/method_c/method_c.py src/single_ancilla/shared.py src/single_ancilla/method_a/method_a.py \
		src/two_ancillas/floquet.py src/two_ancillas/plots.py
	$(PY) src/single_ancilla/method_c/method_c.py --distance $*

# Every readout scheme (spin-ancilla pairs; one ancilla per edge: methods A, B, C) through
# the same memory experiment and decoder: one CSV and one plot, all DISTANCES, colour = distance, line style = scheme.
decode: $(DECODE).png

$(DECODE).png: src/single_ancilla/decode.py src/single_ancilla/method_a/method_a.py src/single_ancilla/shared.py \
		src/single_ancilla/method_b/method_b.py src/single_ancilla/method_c/method_c.py \
		src/two_ancillas/floquet.py src/two_ancillas/noise.py src/two_ancillas/memory.py
	$(PY) src/single_ancilla/decode.py --schemes $(SCHEMES) --distances $(DISTANCES) --noise $(NOISE) --eta $(ETA) \
		--shots $(SHOTS) --p-min $(P_MIN) --p-max $(P_MAX) --num $(NUM) --out $(DECODE)

# LER vs physical qubits at P for each of SCHEMES, and the qubits each needs for LER TARGET.
footprint: $(FOOTPRINT).png

$(FOOTPRINT).png: src/single_ancilla/footprint.py src/single_ancilla/decode.py src/single_ancilla/method_a/method_a.py src/single_ancilla/shared.py \
		src/single_ancilla/method_b/method_b.py \
		src/single_ancilla/method_c/method_c.py src/two_ancillas/floquet.py src/two_ancillas/noise.py src/two_ancillas/memory.py
	$(PY) src/single_ancilla/footprint.py --schemes $(SCHEMES) --distances $(FOOTPRINT_DISTANCES) --noise $(NOISE) --eta $(ETA) \
		--p $(P) --shots $(FOOTPRINT_SHOTS) --target $(TARGET) --out $(FOOTPRINT)

$(RUN)_d%.csv: src/two_ancillas/floquet.py src/two_ancillas/noise.py src/two_ancillas/memory.py
	$(PY) src/two_ancillas/memory.py --distance $* --rounds $* --noise $(NOISE) --eta $(ETA) --shots $(SHOTS) \
		--p-min $(P_MIN) --p-max $(P_MAX) --num $(NUM) --csv $@

$(RUN)_ler.png: $(CSVS) src/two_ancillas/plots.py
	$(PY) src/two_ancillas/plots.py ler $(CSVS) --png $@

# The write-up on building the spin-ancilla Floquet rounds (needs tectonic).
docs: docs/floquet_rounds.pdf docs/spin_x3z3_readout.pdf docs/decoding.pdf

docs/floquet_rounds.pdf: docs/floquet_rounds.tex docs/spin_subround.png $(TWO)/honeycomb_d1/layout.png
	cd docs && tectonic floquet_rounds.tex

# Write-up of the handwritten notes of 6 July 2026 (the scan sits next to it) and ideas for adding them.
docs/spin_x3z3_readout.pdf: docs/spin_x3z3_readout.tex docs/spin_x3z3_notes_2026-07-06.pdf
	cd docs && tectonic spin_x3z3_readout.tex

# How the memory circuit is decoded: DEM, decomposition, matching, with numbers computed from the repo.
docs/decoding.pdf: docs/decoding.tex
	cd docs && tectonic decoding.tex

clean:
	rm -rf .venv

# A rule that fails leaves no target behind, so the next make retries it instead of calling it done.
.DELETE_ON_ERROR:

.PHONY: venv timeslices honeycomb plots experiments decode footprint docs clean
