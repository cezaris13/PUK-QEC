# Floquet honeycomb code with spin-qubit readout. Every setting below can be overridden on the command line:
#   make plots NOISE=MHBD ETA=100 DISTANCES="3 5 7 9" SHOTS=100000
#   make drawings DRAW_DISTANCE=2 SENSORS=
#
# Targets:
#   venv        Python environment in .venv (run once)
#   drawings    every layout and timeslice picture, two ancillas and methods A/B/C
#   plots       two-ancilla LER vs p, one CSV per distance
#   decode      every scheme in SCHEMES through the same memory experiment and decoder
#   footprint   LER vs physical qubits at one p, per scheme
#   thresholds  threshold surface and threshold vs bias for SCHEME
#   docs        the LaTeX write-ups (needs tectonic)
#   serve       interactive viewer of results/ (needs node)
#   clean       remove .venv

PY = .venv/bin/python

# ============================================================================================================
# Parameters
# ============================================================================================================

# --- Memory experiment (plots, decode) ----------------------------------------------------------------------
# Run time grows roughly as SHOTS x NUM x (number of distances) x d^3 (~d^2 qubits for d rounds), so the
# largest distance dominates. notebooks/two_ancillas/memory_experiment.ipynb shows each setting at work.

# Code distances. One CSV and one LER line each; rounds = distance.
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

# --- Readout schemes (decode, footprint) --------------------------------------------------------------------
# pairs (two ancillas per edge), method_a, method_b, method_c (one ancilla per edge, hexes.pdf); see
# src/single_ancilla/decode.py. One output file per set of schemes.
SCHEMES = pairs method_a method_b method_c

# --- Footprint ----------------------------------------------------------------------------------------------
# LER vs physical qubits at this one p, extrapolated to TARGET. Low p needs many shots per point before the
# larger distances see any logical errors to fit.
P = 1e-3
TARGET = 1e-6
FOOTPRINT_SHOTS = 10000000
# Distances for the fit: d = 5 already needs ~10^7 shots at p = 1e-3 to see errors, d = 7 far more.
FOOTPRINT_DISTANCES = 3 5

# --- Thresholds ---------------------------------------------------------------------------------------------
# One scheme at a time: pairs, method_a, method_b or method_c.
SCHEME = pairs
THRESHOLD_DISTANCES = 3 5

# --- Drawings -----------------------------------------------------------------------------------------------
# Lattice size: 12*d^2 data qubits plus 2 spin ancillas on each of the 18*d^2 edges, so 1 is 12 + 36 qubits
# and 2 is 48 + 144; bigger gets hard to read.
DRAW_DISTANCE = 1
# Set NUMBERS=1 to label every qubit in the two-ancilla drawings: Di data, Sk / Rk syndrome / reference
# ancilla. (Methods A/B/C always draw both plain and numbered.)
NUMBERS = 1
# Set SENSORS= (empty) to leave out the charge-sensor squares in the method A/B/C timeslices.
SENSORS = 1

# ============================================================================================================
# Derived names (not meant to be overridden)
# ============================================================================================================

# All output goes under results/, mirroring src/.
TWO = results/two_ancillas
ONE = results/single_ancilla

# $(call dashed,a b c) -> a-b-c: a list of words as one file-name part.
empty :=
space := $(empty) $(empty)
dashed = $(subst $(space),-,$(strip $(1)))

# Output file names carry their settings, so runs with different settings never overwrite each other.
RUN_NAME = $(NOISE)_eta$(ETA)_p$(P_MIN)-$(P_MAX)x$(NUM)_shots$(SHOTS)
RUN = $(TWO)/$(RUN_NAME)
CSVS = $(foreach d,$(DISTANCES),$(RUN)_d$(d).csv)
DECODE = $(ONE)/$(RUN_NAME)_$(call dashed,$(SCHEMES))_decode
FOOTPRINT = $(ONE)/$(NOISE)_eta$(ETA)_p$(P)_d$(call dashed,$(FOOTPRINT_DISTANCES))_shots$(FOOTPRINT_SHOTS)_$(call dashed,$(SCHEMES))_footprint

# Sources each result depends on: change one and make reruns what it feeds.
TWO_SRC = src/two_ancillas/floquet.py src/two_ancillas/noise.py src/two_ancillas/memory.py
ONE_SRC = $(TWO_SRC) src/single_ancilla/shared.py \
	src/single_ancilla/method_a/method_a.py src/single_ancilla/method_b/method_b.py src/single_ancilla/method_c/method_c.py

# ============================================================================================================
# Setup
# ============================================================================================================

venv:
	python3 -m venv .venv
	.venv/bin/pip install -U pip
	.venv/bin/pip install -r requirements.txt
	.venv/bin/python -m ipykernel install --user --name floquet_code

# ============================================================================================================
# Drawings
# ============================================================================================================

# Every drawing at DRAW_DISTANCE, in one go:
# - two ancillas (results/two_ancillas/): one period's timeslices (every TICK a tile, so one period only), and the
#   qubit layout with the 6 steps' coloured timeslices, brick wall and regular hexagons (numbered with NUMBERS);
# - one ancilla, methods A, B and C of hexes.pdf (results/single_ancilla/method_*/d<DRAW_DISTANCE>/): each one's
#   6 steps, plain and numbered, with charge sensors unless SENSORS is empty; method A also draws the layouts.
# Phony, so it always redraws and picks up NUMBERS / SENSORS.
drawings:
	$(PY) src/two_ancillas/plots.py timeslices --distance $(DRAW_DISTANCE) --rounds 1 --png $(TWO)/timeslices_d$(DRAW_DISTANCE).png
	$(PY) src/two_ancillas/plots.py honeycomb --distance $(DRAW_DISTANCE) $(if $(NUMBERS),--numbers) \
		--out $(TWO)/honeycomb$(if $(NUMBERS),_numbered)_d$(DRAW_DISTANCE)
	for m in a b c; do \
		$(PY) src/single_ancilla/method_$$m/method_$$m.py --distance $(DRAW_DISTANCE) --$(if $(SENSORS),,no-)sensors || exit 1; \
	done

# The plain layout docs/floquet_rounds.tex includes.
$(TWO)/honeycomb_d%/layout.png: src/two_ancillas/floquet.py src/two_ancillas/plots.py
	$(PY) src/two_ancillas/plots.py honeycomb --distance $* --out $(@D)

# ============================================================================================================
# Two ancillas per edge: LER vs p
# ============================================================================================================

plots: $(RUN)_ler.png

# One CSV per distance, d rounds each.
$(RUN)_d%.csv: $(TWO_SRC)
	$(PY) src/two_ancillas/memory.py --distance $* --rounds $* --noise $(NOISE) --eta $(ETA) --shots $(SHOTS) \
		--p-min $(P_MIN) --p-max $(P_MAX) --num $(NUM) --csv $@

$(RUN)_ler.png: $(CSVS) src/two_ancillas/plots.py
	$(PY) src/two_ancillas/plots.py ler $(CSVS) --png $@

# ============================================================================================================
# Comparing readout schemes
# ============================================================================================================

# Every scheme in SCHEMES through the same memory experiment and decoder: one CSV and one plot, all
# DISTANCES, colour = distance, line style = scheme.
decode: $(DECODE).png

$(DECODE).png: src/single_ancilla/decode.py $(ONE_SRC)
	$(PY) src/single_ancilla/decode.py --schemes $(SCHEMES) --distances $(DISTANCES) --noise $(NOISE) --eta $(ETA) \
		--shots $(SHOTS) --p-min $(P_MIN) --p-max $(P_MAX) --num $(NUM) --out $(DECODE)

# LER vs physical qubits at P for each of SCHEMES, and the qubits each needs for LER TARGET.
footprint: $(FOOTPRINT).png

$(FOOTPRINT).png: src/single_ancilla/footprint.py src/single_ancilla/decode.py $(ONE_SRC)
	$(PY) src/single_ancilla/footprint.py --schemes $(SCHEMES) --distances $(FOOTPRINT_DISTANCES) --noise $(NOISE) --eta $(ETA) \
		--p $(P) --shots $(FOOTPRINT_SHOTS) --target $(TARGET) --out $(FOOTPRINT)

# Threshold surface (b) and threshold vs bias (c) after IBM's QEC-with-spin-qubits, for SCHEME, with their
# bias-blind decoder. Writes results/two_ancillas/thresholds/ (pairs) or results/single_ancilla/<SCHEME>/thresholds/.
# Phony: the script names the files from its settings.
thresholds:
	$(PY) src/single_ancilla/thresholds.py --scheme $(SCHEME) --distances $(THRESHOLD_DISTANCES)

# ============================================================================================================
# Write-ups (needs tectonic)
# ============================================================================================================

docs: docs/floquet_rounds.pdf docs/spin_x3z3_readout.pdf docs/decoding.pdf

# Building the spin-ancilla Floquet rounds.
docs/floquet_rounds.pdf: docs/floquet_rounds.tex docs/spin_subround.png $(TWO)/honeycomb_d1/layout.png
	cd docs && tectonic floquet_rounds.tex

# The handwritten notes of 6 July 2026 (the scan sits next to it) and ideas for adding them.
docs/spin_x3z3_readout.pdf: docs/spin_x3z3_readout.tex docs/spin_x3z3_notes_2026-07-06.pdf
	cd docs && tectonic spin_x3z3_readout.tex

# How the memory circuit is decoded: DEM, decomposition, matching, with numbers computed from the repo.
docs/decoding.pdf: docs/decoding.tex
	cd docs && tectonic decoding.tex

# ============================================================================================================
# Viewer and housekeeping
# ============================================================================================================

# Interactive viewer of everything under results/: LER vs p for every scheme, threshold surfaces, footprints.
# A Vite + React app in viewer/ (needs node); its dev server also serves results/.
serve: viewer/node_modules
	cd viewer && npm run dev

viewer/node_modules: viewer/package.json
	cd viewer && npm install
	touch $@

clean:
	rm -rf .venv

# A rule that fails leaves no target behind, so the next make retries it instead of calling it done.
.DELETE_ON_ERROR:

.PHONY: venv drawings plots decode footprint thresholds docs serve clean
