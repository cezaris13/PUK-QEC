PY = .venv/bin/python

# Sweep settings. Override any on the command line:
#   make plots NOISE=MHBD ETA=100 DISTANCES="3 5 7 9" SHOTS=100000
# Run time grows roughly as SHOTS x NUM x (number of distances) x d^3 (~d^2 qubits for d rounds),
# so the largest distance dominates. memory_experiment.ipynb shows each setting at work.

# Code distances. One CSV and one LER line each; rounds = distance. `make timeslices` draws the first.
DISTANCES = 3 5 7
# Noise model from NOISE_MODELS in src/noise.py: spin, HBD or MHBD.
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
# Set NUMBERS=1 to label every qubit in the layout drawings: Di data, Sk / Rk syndrome / reference ancilla.
NUMBERS =

RUN = results/$(NOISE)_eta$(ETA)_p$(P_MIN)-$(P_MAX)x$(NUM)_shots$(SHOTS)
CSVS = $(foreach d,$(DISTANCES),$(RUN)_d$(d).csv)
TIMESLICES = results/timeslices_d$(firstword $(DISTANCES)).png
HONEYCOMB = results/honeycomb$(if $(NUMBERS),_numbered)_d$(DRAW_DISTANCE)/layout.png

venv:
	python3 -m venv .venv
	.venv/bin/pip install -U pip
	.venv/bin/pip install -r requirements.txt
	.venv/bin/python -m ipykernel install --user --name floquet_code

timeslices: $(TIMESLICES)

# Qubit layout, plus the 3 sub-rounds' coloured timeslices in brick-wall and regular-hexagon views.
honeycomb: $(HONEYCOMB)

plots: $(TIMESLICES) $(HONEYCOMB) $(RUN)_ler.png

results/timeslices_d%.png: src/floquet.py src/plots.py
	$(PY) src/plots.py timeslices --distance $* --rounds $* --png $@

# layout.png stands in for the whole folder: the round PNGs are written alongside it.
results/honeycomb_d%/layout.png: src/floquet.py src/plots.py
	$(PY) src/plots.py honeycomb --distance $* --out $(@D)

# The same drawings with numbered qubits, in their own folder so neither ever passes for the other.
results/honeycomb_numbered_d%/layout.png: src/floquet.py src/plots.py
	$(PY) src/plots.py honeycomb --distance $* --numbers --out $(@D)

$(RUN)_d%.csv: src/floquet.py src/noise.py src/memory.py
	$(PY) src/memory.py --distance $* --rounds $* --noise $(NOISE) --eta $(ETA) --shots $(SHOTS) \
		--p-min $(P_MIN) --p-max $(P_MAX) --num $(NUM) --csv $@

$(RUN)_ler.png: $(CSVS) src/plots.py
	$(PY) src/plots.py ler $(CSVS) --png $@

clean:
	rm -rf .venv

# A rule that fails leaves no target behind, so the next make retries it instead of calling it done.
.DELETE_ON_ERROR:

.PHONY: venv timeslices honeycomb plots docs clean
