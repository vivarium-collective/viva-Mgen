"""Fitted constants for viva-mGen, taken verbatim from the Karr et al. 2012
whole-cell model's ``data/parameters.json`` (states.Mass / states.Time /
states.Metabolite) and the *M. genitalium* genome.

These are the *real* fitted values from the published model — the one part of
the original knowledge base that is pure-JSON-loadable (the reaction/species/
kcat tables are locked inside undecodable MATLAB MCOS objects; the metabolic
network is supplied instead by the published Suthers 2009 iPS189
reconstruction, see :mod:`viva_mgen.kb`).
"""

from __future__ import annotations

# --- genome ---------------------------------------------------------------
GENOME_LENGTH_BP = 580070          # M. genitalium G37 chromosome length (bp)
N_GENES = 525                      # protein + RNA genes in the reconstruction

# --- cell mass (states.Mass) ---------------------------------------------
CELL_INITIAL_DRY_WEIGHT_G = 3.93e-15   # g dry weight of a newborn cell
FRACTION_WET_WEIGHT = 0.7              # water is 70% of wet weight
CELL_INITIAL_DRY_WEIGHT_FG = CELL_INITIAL_DRY_WEIGHT_G * 1e15   # = 3.93 fg
CELL_INITIAL_WET_WEIGHT_FG = CELL_INITIAL_DRY_WEIGHT_FG / (1.0 - FRACTION_WET_WEIGHT)

# dry-weight mass fractions (sum ~1.0)
DRY_MASS_FRACTIONS = {
    "DNA": 0.1688,
    "RNA": 0.0929552671058823,
    "protein": 0.619701780705882,
    "lipid": 0.0568059965647059,
    "carbohydrate": 0.00516418150588235,
    "nucleotide": 0.01833748,
    "ion": 0.00980392156862745,
    "polyamine": 0.01323529411764706,
    "vitamin": 0.01519607843137255,
}

# --- cell-cycle timing (states.Time), seconds ----------------------------
CELL_CYCLE_LENGTH_S = 32400.0          # ~9.0 h — the model's mean doubling time
REPLICATION_INITIATION_DURATION_S = 12960.0
REPLICATION_DURATION_S = 15571.0
CYTOKINESIS_DURATION_S = 3869.0

# --- metabolite mean concentrations (states.Metabolite), mM --------------
MEAN_NTP_CONC_MM = 5.0
MEAN_NDP_CONC_MM = 0.5
MEAN_NMP_CONC_MM = 0.2

# --- metabolism (processes.Metabolism) -----------------------------------
GROWTH_ASSOCIATED_MAINTENANCE = 59.81
NON_GROWTH_ASSOCIATED_MAINTENANCE = 8.39
EXCHANGE_UB_CARBON = 12.0
EXCHANGE_UB_NONCARBON = 20.0

# --- simulation ----------------------------------------------------------
STEP_SIZE_SEC = 1.0                    # fixed 1 s timestep (as in the original)
AVOGADRO = 6.022140857e23
