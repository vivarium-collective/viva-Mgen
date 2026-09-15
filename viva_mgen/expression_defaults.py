"""Genome-scale gene set + parameters for the stochastic expression processes
(transcription / translation / decay).

The whole-cell model tracks all ~480 protein-coding genes with per-gene rates.
This reproduction follows the SAME stochastic mechanism (Poisson/multinomial
synthesis, Poisson decay) over the FULL M. genitalium gene set loaded from
``datasets/genes.csv`` (~522 genes). Per-gene rates are computed by the native
ParCa (:mod:`viva_mgen.parca`, ``_parca_panel``) from the REAL observed
expression profile and REAL KB per-gene half-lives — the knowledge base is
decoded natively by :mod:`viva_mgen.kb_decode`, no MATLAB required. If the
decoded expression is unavailable the module falls back to a heavy-tailed
(log-normal) draw calibrated to the paper's expression kinetics (~50% of genes
expressed within 18 min, ~90% within 143 min, Fig 3C). A handful of well-known
genes keep hand-tuned rates. Protein half-lives are a single representative
value because the KB carries no per-monomer protein half-lives.
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


def _parca_panel() -> dict:
    """Per-gene panel computed by the native ParCa from the REAL observed
    expression (viva_mgen.parca), preferred over the log-normal draw. The named
    genes keep their hand-tuned rates."""
    from .parca import calculate_parameters
    panel = calculate_parameters()
    panel.update(_NAMED)   # keep the well-characterized genes exact
    return panel


try:
    # prefer the ParCa (real observed expression); fall back to the log-normal
    # panel, then to the named genes, so the module always yields a usable panel.
    try:
        REPRESENTATIVE_GENES = _parca_panel()
    except Exception:  # noqa: BLE001 — no decoded expression available
        REPRESENTATIVE_GENES = _build_panel()
except Exception:  # noqa: BLE001 — fall back to the named panel if genes.csv is unreadable
    REPRESENTATIVE_GENES = dict(_NAMED)

DEFAULT_GENES = list(REPRESENTATIVE_GENES.keys())


# Whole-cell stable-RNA synthesis calibration. In the emergent inventory the RNA
# mass is dominated by the STABLE RNAs — rRNA (~73%) and tRNA (~19%); real mRNA is
# only ~9%. Those stable species have no packaging/charging sink here (rRNA is not
# consumed into ribosomes, tRNA not sequestered while charged), and their
# half-lives (rRNA ~20 h, tRNA ~45 min) exceed the ~9 h cycle, so they accumulate
# roughly linearly instead of reaching a bounded steady state — overshooting the
# real cell's stable-RNA content ~2.7×. This factor rescales rRNA/tRNA synthesis
# (a proxy for that missing sink) so the emergent protein:DNA:RNA dry-mass
# fractions land on Karr's ~0.70:0.19:0.11 (fig2 report-card bands). mRNA synthesis
# is left at full rate, so gene-expression timing (fig3 t50/t90) and the mRNA pool
# feeding translation are unchanged. Fitted jointly with metabolism.gtp_base_supply.
# See docs/FIDELITY_GAPS.md gap #2.
STABLE_RNA_SYNTHESIS_SCALE = 0.37


def _stable_rna_keys() -> frozenset:
    """Panel keys whose gene product is a non-coding (stable) RNA — rRNA, tRNA, or
    the SRP RNA — classified by :func:`viva_mgen.parca.rna_type` (which matches the
    RNA products themselves, not the proteins that process them). Cached; empty if
    the gene table is unreadable."""
    cached = getattr(_stable_rna_keys, "_cache", None)
    if cached is not None:
        return cached
    keys = set()
    try:
        from .kb import load_genes
        from .parca import rna_type
        for g in load_genes():
            if rna_type(g) != "mRNA":
                keys.add((g.get("symbol") or "").strip() or g["gene_id"])
    except Exception:  # noqa: BLE001 — no gene table: classify nothing (mRNA-only panel)
        keys = set()
    _stable_rna_keys._cache = frozenset(keys)
    return _stable_rna_keys._cache


def synthesis_rates() -> dict:
    stable = _stable_rna_keys()
    return {g: v[0] * (STABLE_RNA_SYNTHESIS_SCALE if g in stable else 1.0)
            for g, v in REPRESENTATIVE_GENES.items()}


def unscaled_synthesis_rates() -> dict:
    """Per-gene mRNA-synthesis rates WITHOUT STABLE_RNA_SYNTHESIS_SCALE. The scale
    is a proxy for the missing stable-RNA degradation/packaging sink — it lowers
    the emergent rRNA/tRNA POOL, not the true transcription RATE. For quantities
    that depend on how fast the cell actually transcribes (e.g. the Fig 5 energy
    budget, where the cell really does synthesize rRNA/tRNA at full rate and then
    degrades them), use these full rates rather than the pool-calibrated ones."""
    return {g: v[0] for g, v in REPRESENTATIVE_GENES.items()}


def mrna_decay_rates() -> dict:
    return {g: math.log(2.0) / v[1] for g, v in REPRESENTATIVE_GENES.items()}


def translation_rates() -> dict:
    # Only mRNA (protein-coding) genes are translated. The panel carries a nominal
    # translation rate for every gene, but ribosomes translate mRNA — not the tRNA/
    # rRNA/SRP-RNA gene products — so non-coding genes are excluded here (otherwise
    # e.g. a 76-nt tRNA gene would be "translated" into a spurious 25-aa protein,
    # which dominated the emergent protein count). Their transcripts still count as
    # RNA via rna_counts; they simply have no protein product.
    stable = _stable_rna_keys()
    return {g: v[2] for g, v in REPRESENTATIVE_GENES.items() if g not in stable}


def protein_decay_rates() -> dict:
    return {g: math.log(2.0) / v[3] for g, v in REPRESENTATIVE_GENES.items()}


def gene_lengths() -> dict:
    return {g: float(v[4]) for g, v in REPRESENTATIVE_GENES.items()}


def reference_protein_counts() -> dict:
    """Expected steady-state protein count per gene, from the panel:
    mRNA_ss = synthesis_rate / mrna_decay_rate; protein_ss = translation_rate * mRNA_ss / protein_decay_rate.
    Keyed like protein_counts (gene symbol/id). Used as the reference the dynamic
    metabolism coupling scales the live proteome against."""
    synth = synthesis_rates()
    mdec = mrna_decay_rates()
    transl = translation_rates()
    pdec = protein_decay_rates()
    out = {}
    for g in synth:
        md = mdec.get(g, 0.0)
        pd = pdec.get(g, 0.0)
        if md <= 0 or pd <= 0:
            out[g] = 0.0
            continue
        mrna_ss = synth[g] / md
        out[g] = transl.get(g, 0.0) * mrna_ss / pd
    return out


# Emergent whole-cell proteome scale (gap #6). The analytic reference_protein_counts
# above (~3.0M total = steady-state translation·mRNA/decay) far exceeds the
# GTP-throttled EMERGENT proteome the integrated cell actually reaches (~135k at
# division). Two consumers need the emergent scale, not the analytic one:
#   * the BIRTH proteome a daughter is seeded with, so the cell cycle is
#     birth→double rather than accumulate-from-zero; and
#   * enzyme_coupling's live/reference ratio, which must be ~1 at the normal
#     proteome (else every reaction clips low and metabolism starves).
# So scale the analytic distribution down to the emergent totals. Fitted jointly
# with metabolism.gtp_base_supply. See docs/FIDELITY_GAPS.md gap #6.
BIRTH_PROTEOME_SCALE = 0.0225  # analytic reference -> ~half the emergent division proteome (~67.5k)


def birth_proteome() -> dict:
    """The proteome a daughter cell is born with: ~half the emergent division
    proteome (~67.5k), keyed like protein_counts, in the analytic distribution
    scaled to the emergent total. Seeding it makes the cycle birth→double."""
    return {g: v * BIRTH_PROTEOME_SCALE for g, v in reference_protein_counts().items()}


def coupling_reference() -> dict:
    """Reference proteome for enzyme_coupling's live/reference ratio: the emergent
    DIVISION proteome (~135k = 2×birth), so live/ref ≈ 0.5 at birth and ≈ 1.0 near
    division and reactions scale gracefully with enzyme abundance."""
    return {g: 2.0 * v for g, v in birth_proteome().items()}
