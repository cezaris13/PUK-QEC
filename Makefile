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
#   thresholds-all    thresholds for every scheme in SCHEMES; thresholds-medium / thresholds-paper the same with more shots
#   thresholds-x3z3   thresholds-all for the X3Z3 code (any of them takes CODE=x3z3 too)
#   thresholds-x3z3-medium  thresholds-medium for the X3Z3 code
#   x3z3        X3Z3 vs CSS Floquet code for every scheme: threshold vs bias, LER vs p, LER vs qubits
#   serve       interactive viewer of results/ (needs node)
#   clean       remove .venv

PY = .venv/bin/python

# ============================================================================================================
# Parameters
# ============================================================================================================

# --- Memory experiment (plots, decode) ----------------------------------------------------------------------

# Code distances. One CSV and one LER line each; rounds = distance.
DISTANCES = 3 5 7
# Noise model from NOISE_MODELS in src/shared/noise.py: spin, HBD or MHBD.
NOISE = spin
# Bias eta = p_z / (p_x + p_y): 0.5 is depolarizing, large is mostly dephasing. spin biases every
# channel; HBD/MHBD bias two-qubit gates and per-round idling, but keep 1-qubit gates depolarizing.
ETA = 10
# Monte Carlo shots per (distance, p) point. Smallest resolvable LER ~ 1/SHOTS; a point's relative
# error ~ 1/sqrt(logical errors), so 5 errors is +-45%.
SHOTS = 50000
# Physical error rate range. p is the two-qubit gate error; the model scales everything else off it
# (spin: readout 5p, reset 2p, idling per round 2p, 1-qubit gates p/10).
P_MIN = 1e-4
P_MAX = 1e-2
# Number of p values from P_MIN to P_MAX, both included, log-spaced: 9 over two decades = 4 per decade.
NUM = 9

# --- Readout schemes (decode, footprint) --------------------------------------------------------------------
# pairs (two ancillas per edge), method_a, method_b, method_c (one ancilla per edge, hexes.pdf)
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
# Floquet code: css, or x3z3 (Setiawan & McLauchlan); x3z3 files are named x3z3_<scheme>_...
CODE = css
THRESHOLD_DISTANCES = 3 5
# Grid of the threshold surface: NPHI (NPHI + 1) / 2 directions, each one a full scan, so run time grows ~NPHI^2.
NPHI = 6
# BIAS_ONLY=1 reruns only the bias plot (c) and keeps the surface (b) already saved for these settings.
BIAS_ONLY =
# Points of the bias plot: NBIAS values of eta, log-spaced from 0.01 to 100.
NBIAS = 9
# Error rates per direction of the surface (NUM_P) and of the bias plot (BIAS_NUM_P).
NUM_P = 10
BIAS_NUM_P = 10
# Each point runs batches of THRESHOLD_BATCH shots until THRESHOLD_SHOTS or MAX_FAIL logical errors.
THRESHOLD_SHOTS = 20000
THRESHOLD_BATCH = 2000
MAX_FAIL = 2000
# IBM's paper (arXiv:2306.17786) used NPHI=20 NBIAS~25 NUM_P=30 BIAS_NUM_P=60 THRESHOLD_SHOTS=300000 THRESHOLD_BATCH=3000
# MAX_FAIL=30000: ~100x the defaults' run time.

# --- X3Z3 vs CSS (x3z3) -------------------------------------------------------------------------------------
# Every scheme in SCHEMES, both codes, under NOISE at each of X3Z3_ETAS; p from X3Z3_P_MIN to X3Z3_P_MAX
# (X3Z3_NUM log-spaced) plus P, the p of the LER-vs-qubits figure. SHOTS per point and basis.
X3Z3_DISTANCES = 3 5
X3Z3_ETAS = 0.5 1 3 10 30 100 1000
X3Z3_P_MIN = 2e-4
X3Z3_P_MAX = 1.5e-2
X3Z3_NUM = 12

# --- Drawings -----------------------------------------------------------------------------------------------
# Lattice size: 12*d^2 data qubits plus 2 spin ancillas on each of the 18*d^2 edges, so 1 is 12 + 36 qubits
# and 2 is 48 + 144; bigger gets hard to read.
DRAW_DISTANCE = 1
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

RUN_NAME = $(NOISE)_eta$(ETA)_p$(P_MIN)-$(P_MAX)x$(NUM)_shots$(SHOTS)
RUN = $(TWO)/$(RUN_NAME)
CSVS = $(foreach d,$(DISTANCES),$(RUN)_d$(d).csv)
DECODE = $(ONE)/$(RUN_NAME)_$(call dashed,$(SCHEMES))_decode
FOOTPRINT = $(ONE)/$(NOISE)_eta$(ETA)_p$(P)_d$(call dashed,$(FOOTPRINT_DISTANCES))_shots$(FOOTPRINT_SHOTS)_$(call dashed,$(SCHEMES))_footprint

SHARED_SRC = src/shared/floquet.py src/shared/noise.py src/shared/memory.py
TWO_SRC = $(SHARED_SRC) src/two_ancillas/pairs.py
ONE_SRC = $(TWO_SRC) src/single_ancilla/shared/layout.py src/single_ancilla/shared/corner_readout.py \
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

# Every drawing at DRAW_DISTANCE:
# - two ancillas: one period's timeslices in
#   results/two_ancillas/timeslices_d<DRAW_DISTANCE>.png, and the qubit layout with the 6 steps' coloured
#   timeslices in results/two_ancillas/honeycomb_d<DRAW_DISTANCE>/;
# - one ancilla, methods A, B and C of hexes.pdf: each one's 6 steps, with charge sensors unless SENSORS is
#   empty, in results/single_ancilla/method_*/d<DRAW_DISTANCE>/; method A also draws the layout.
# Each folder splits into rectangle/ (brick wall) and hex/ (regular hexagons), each into plain/ and numbered/
# (every qubit labelled: Di data, Sk / Rk syndrome / reference ancilla), holding layout.png and round0-5.png.
# Phony, so it always redraws and picks up SENSORS.
drawings:
	$(PY) src/two_ancillas/plots.py timeslices --distance $(DRAW_DISTANCE) --rounds 1 --png $(TWO)/timeslices_d$(DRAW_DISTANCE).png
	$(PY) src/two_ancillas/plots.py honeycomb --distance $(DRAW_DISTANCE) --out $(TWO)/honeycomb_d$(DRAW_DISTANCE)
	for m in a b c; do \
		$(PY) src/single_ancilla/method_$$m/method_$$m.py --distance $(DRAW_DISTANCE) --$(if $(SENSORS),,no-)sensors || exit 1; \
	done

# The plain layout docs/floquet_rounds.tex includes.
$(TWO)/honeycomb_d%/rectangle/plain/layout.png: src/shared/floquet.py src/shared/drawing.py src/two_ancillas/pairs.py \
		src/two_ancillas/plots.py
	$(PY) src/two_ancillas/plots.py honeycomb --distance $* --out $(TWO)/honeycomb_d$*

# ============================================================================================================
# Two ancillas per edge: LER vs p
# ============================================================================================================

plots: $(RUN)_ler.png

# One CSV per distance, d rounds each.
# Always run: memory.py appends each p as it finishes and skips those already in the CSV, so a cut-short run picks
# up where it stopped (it warns if TWO_SRC changed since; delete the CSV to start over). Precious: kept on Ctrl-C.
$(RUN)_d%.csv: FORCE
	$(PY) src/shared/memory.py --distance $* --rounds $* --noise $(NOISE) --eta $(ETA) --shots $(SHOTS) \
		--p-min $(P_MIN) --p-max $(P_MAX) --num $(NUM) --csv $@ --sources $(TWO_SRC)
.PRECIOUS: $(RUN)_d%.csv
FORCE:

$(RUN)_ler.png: $(CSVS) src/shared/drawing.py src/two_ancillas/plots.py
	$(PY) src/two_ancillas/plots.py ler $(CSVS) --png $@

# ============================================================================================================
# Comparing readout schemes
# ============================================================================================================

# Every scheme in SCHEMES through the same memory experiment and decoder: one CSV and one plot, all
# DISTANCES, colour = distance, line style = scheme.
decode: $(DECODE).png

$(DECODE).png: src/single_ancilla/decode.py $(ONE_SRC)
	$(PY) src/single_ancilla/decode.py --schemes $(SCHEMES) --distances $(DISTANCES) --noise $(NOISE) --eta $(ETA) \
		--shots $(SHOTS) --p-min $(P_MIN) --p-max $(P_MAX) --num $(NUM) --out $(DECODE) --sources $(ONE_SRC)

# LER vs physical qubits at P for each of SCHEMES, and the qubits each needs for LER TARGET.
footprint: $(FOOTPRINT).png

$(FOOTPRINT).png: src/single_ancilla/footprint.py src/single_ancilla/decode.py $(ONE_SRC)
	$(PY) src/single_ancilla/footprint.py --schemes $(SCHEMES) --distances $(FOOTPRINT_DISTANCES) --noise $(NOISE) --eta $(ETA) \
		--p $(P) --shots $(FOOTPRINT_SHOTS) --target $(TARGET) --out $(FOOTPRINT) --sources $(ONE_SRC)

# Threshold surface (b) and threshold vs bias (c) after IBM's QEC-with-spin-qubits, for SCHEME, with their
# bias-blind decoder. Writes results/two_ancillas/thresholds/ (pairs) or results/single_ancilla/<SCHEME>/thresholds/.
# Phony: the script names the files from its settings.
thresholds:
	$(PY) src/single_ancilla/thresholds.py --scheme $(SCHEME) --code $(CODE) --distances $(THRESHOLD_DISTANCES) --nphi $(NPHI) --nbias $(NBIAS) \
		--num-p $(NUM_P) --bias-num-p $(BIAS_NUM_P) --shots $(THRESHOLD_SHOTS) --batch $(THRESHOLD_BATCH) --max-fail $(MAX_FAIL) $(if $(BIAS_ONLY),--bias-only) \
		--sources $(ONE_SRC)

# thresholds for every scheme in SCHEMES, one after the other; stops at the first that fails.
thresholds-all:
	for s in $(SCHEMES); do $(MAKE) thresholds SCHEME=$$s || exit 1; done

# thresholds-all for the X3Z3 code. For IBM's settings: make thresholds-paper CODE=x3z3.
thresholds-x3z3:
	$(MAKE) thresholds-all CODE=x3z3

# thresholds-medium for the X3Z3 code. Several hours per scheme.
thresholds-x3z3-medium:
	$(MAKE) thresholds-medium CODE=x3z3

# thresholds-all between the defaults and IBM's settings. Several hours per scheme.
thresholds-medium:
	$(MAKE) thresholds-all NPHI=15 NBIAS=17 NUM_P=20 BIAS_NUM_P=40 THRESHOLD_SHOTS=100000 THRESHOLD_BATCH=2000 MAX_FAIL=10000

# thresholds-all with IBM's settings (see THRESHOLD_SHOTS above). Roughly a day per scheme.
thresholds-paper:
	$(MAKE) thresholds-all NPHI=20 NBIAS=25 NUM_P=30 BIAS_NUM_P=60 THRESHOLD_SHOTS=300000 THRESHOLD_BATCH=3000 MAX_FAIL=30000

# X3Z3 (Setiawan & McLauchlan, arXiv:2411.04974) against the CSS Floquet code for every scheme in SCHEMES:
# results/single_ancilla/x3z3/<settings>.csv, and _bias.png (threshold and LER vs eta), _ler_eta<eta>.png and
# _qubits_eta<eta>.png per eta. Phony: the script names its files and resumes from the CSV.
x3z3:
	$(PY) src/single_ancilla/x3z3.py --schemes $(SCHEMES) --noise $(NOISE) --distances $(X3Z3_DISTANCES) \
		--etas $(X3Z3_ETAS) --p-min $(X3Z3_P_MIN) --p-max $(X3Z3_P_MAX) --num $(X3Z3_NUM) --p-fixed $(P) \
		--shots $(SHOTS) --out $(ONE)/x3z3/$(NOISE)_d$(call dashed,$(X3Z3_DISTANCES))_shots$(SHOTS)_$(call dashed,$(SCHEMES))

# ============================================================================================================
# Viewer and housekeeping
# ============================================================================================================

serve: viewer/node_modules
	cd viewer && npm run dev

viewer/node_modules: viewer/package.json
	cd viewer && npm install
	touch $@

clean:
	rm -rf .venv

.DELETE_ON_ERROR:

.PHONY: venv drawings plots decode footprint thresholds thresholds-all thresholds-medium thresholds-paper thresholds-x3z3 thresholds-x3z3-medium x3z3 docs serve clean
