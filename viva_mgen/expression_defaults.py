"""Representative gene set + parameters for the reduced stochastic expression
processes (transcription / translation / decay).

The original whole-cell model tracks all ~480 protein-coding genes with
per-species rates drawn from the (undecodable-here) knowledge base. This
reduced reproduction follows the SAME stochastic mechanism (Poisson/multinomial
synthesis and Poisson decay) on a small, representative panel of genes with
order-of-magnitude-genuine parameters, sufficient to reproduce the *qualitative*
single-cell expression phenomena of Fig 2G/2H (bursty low-copy mRNA, stable
higher-copy protein). Parameters are representative, not the exact KB values —
this is the reduced-fidelity part of the reproduction.
"""

from __future__ import annotations

# symbol: (mrna_synthesis_rate_per_s, mrna_halflife_s, translation_rate_per_mrna_per_s, protein_halflife_s, avg_gene_length_nt)
REPRESENTATIVE_GENES = {
    "HMW2":  (0.0055, 170.0, 0.020, 25000.0, 5000),   # MG_218 cytadherence
    "dnaA":  (0.0040, 150.0, 0.035, 25000.0, 1300),   # MG_469 replication initiator
    "ftsZ":  (0.0060, 160.0, 0.045, 25000.0, 1100),   # MG_224 division
    "rpoD":  (0.0035, 150.0, 0.030, 25000.0, 1600),   # sigma factor
    "tuf":   (0.0130, 180.0, 0.090, 25000.0, 1200),   # EF-Tu (abundant)
    "groEL": (0.0090, 170.0, 0.060, 25000.0, 1600),   # chaperonin (abundant)
}

DEFAULT_GENES = list(REPRESENTATIVE_GENES.keys())


def synthesis_rates() -> dict:
    return {g: v[0] for g, v in REPRESENTATIVE_GENES.items()}


def mrna_decay_rates() -> dict:
    # decay rate = ln2 / half-life
    import math
    return {g: math.log(2.0) / v[1] for g, v in REPRESENTATIVE_GENES.items()}


def translation_rates() -> dict:
    return {g: v[2] for g, v in REPRESENTATIVE_GENES.items()}


def protein_decay_rates() -> dict:
    import math
    return {g: math.log(2.0) / v[3] for g, v in REPRESENTATIVE_GENES.items()}


def gene_lengths() -> dict:
    return {g: float(v[4]) for g, v in REPRESENTATIVE_GENES.items()}
