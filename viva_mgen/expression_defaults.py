"""Genome-scale gene set + parameters for the reduced stochastic expression
processes (transcription / translation / decay).

The original whole-cell model tracks all ~480 protein-coding genes with per-species
rates from its (undecodable-here) knowledge base. This reproduction follows the
SAME stochastic mechanism (Poisson/multinomial synthesis, Poisson decay) over the
FULL M. genitalium gene set loaded from ``datasets/genes.csv``, with per-gene
rates drawn from a realistic heavy-tailed (log-normal) distribution calibrated so
the ensemble reproduces the paper's expression kinetics — ~50% of genes expressed
within 18 min and ~90% within 143 min (Fig 3C) — and bursty low-copy single-cell
mRNA (Fig 2G). A handful of well-known genes keep hand-tuned rates; the rest are
distributed. Rates are representative (order-of-magnitude genuine), not the exact
KB values — this is the reduced-fidelity part of the reproduction.
"""

from __future__ import annotations

import math

# hand-tuned rates for a few well-known genes (kept exactly):
# symbol: (mrna_synthesis_rate_per_s, mrna_halflife_s, translation_rate_per_mrna_per_s, protein_halflife_s, avg_gene_length_nt)
_NAMED = {
    "HMW2":  (0.0055, 170.0, 0.020, 25000.0, 5000),   # MG_218 cytadherence
    "dnaA":  (0.0040, 150.0, 0.035, 25000.0, 1300),   # MG_469 replication initiator
    "ftsZ":  (0.0060, 160.0, 0.045, 25000.0, 1100),   # MG_224 division
    "rpoD":  (0.0035, 150.0, 0.030, 25000.0, 1600),   # sigma factor
    "tuf":   (0.0130, 180.0, 0.090, 25000.0, 1200),   # EF-Tu (abundant)
    "groEL": (0.0090, 170.0, 0.060, 25000.0, 1600),   # chaperonin (abundant)
}

# Log-normal synthesis-rate distribution (per second) for the remaining genes.
# Median ~5.5e-4/s → time-to-first-transcript median ~18 min; sigma gives a heavy
# tail of rarely-transcribed genes so ~10% remain unexpressed past ~143 min.
_LOGN_MU = math.log(5.5e-4)
_LOGN_SIGMA = 1.15


def _build_panel() -> dict:
    """symbol/id -> (synth_rate, mrna_halflife, transl_rate, prot_halflife, length)
    for every gene in the knowledge base."""
    import numpy as np
    from .kb import load_genes
    genes = load_genes()
    rng = np.random.default_rng(20120720)  # deterministic; date of the Karr 2012 paper
    panel: dict = {}
    for g in genes:
        key = (g.get("symbol") or "").strip() or g["gene_id"]
        if key in _NAMED:
            panel[key] = _NAMED[key]
            continue
        synth = float(np.exp(rng.normal(_LOGN_MU, _LOGN_SIGMA)))
        half = float(rng.uniform(140.0, 190.0))                 # bacterial mRNA t½ ~2.5-3 min
        transl = float(np.clip(synth * rng.uniform(4.0, 10.0), 0.005, 0.15))
        plen = int(rng.uniform(700, 2200))
        panel[key] = (synth, half, transl, 25000.0, plen)
    # ensure the named genes are present even if genes.csv lacks the symbol
    for k, v in _NAMED.items():
        panel.setdefault(k, v)
    return panel


try:
    REPRESENTATIVE_GENES = _build_panel()
except Exception:  # noqa: BLE001 — fall back to the named panel if genes.csv is unreadable
    REPRESENTATIVE_GENES = dict(_NAMED)

DEFAULT_GENES = list(REPRESENTATIVE_GENES.keys())


def synthesis_rates() -> dict:
    return {g: v[0] for g, v in REPRESENTATIVE_GENES.items()}


def mrna_decay_rates() -> dict:
    return {g: math.log(2.0) / v[1] for g, v in REPRESENTATIVE_GENES.items()}


def translation_rates() -> dict:
    return {g: v[2] for g, v in REPRESENTATIVE_GENES.items()}


def protein_decay_rates() -> dict:
    return {g: math.log(2.0) / v[3] for g, v in REPRESENTATIVE_GENES.items()}


def gene_lengths() -> dict:
    return {g: float(v[4]) for g, v in REPRESENTATIVE_GENES.items()}
