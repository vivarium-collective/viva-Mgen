"""Parameter calculator (ParCa) for viva-Mgen — a native reproduction of the
Karr 2012 ``FitConstants`` fitting, computing per-gene expression/decay/synthesis
parameters from the experimentally observed inputs instead of hand-tuned guesses.

This mirrors the initialization stage of
``+edu/+stanford/+covert/+cell/+sim/+util/FitConstants.m``:

* **observed gene expression** — the Weiner et al. 2003 transcription profiles
  decoded from the knowledge base (:func:`viva_mgen.kb.load_gene_expression`),
  imputed for missing genes and normalized to a distribution (``geneExp``);
* **decay rates** — from the knowledge base's per-RNA-type half-life scheme
  (mRNA ≈ 2.5 min, tRNA ≈ 45 min, rRNA ≈ 1200 min — the values recovered from the
  KB, grouped by RNA class), ``decay = ln 2 / halfLife``;
* **production/synthesis** — the FitConstants "match expression, decay rates"
  balance: to hold a species at its expression against dilution + decay, the
  production rate is ``expression × (ln2 / cellCycleLength + decay)``.

The result is the per-gene parameter panel the stochastic transcription /
translation / decay processes consume — the same shape as the legacy
``expression_defaults`` panel, but *computed from real data*. Full FitConstants
fidelity (the analytic/QP fit with the RNA-mass / NMP / AA constraints and the
FtsZ/DnaA/topoisomerase held constraints) is layered on top of this initialization
in follow-up work; this module is the faithful initialization stage.
"""
from __future__ import annotations

import functools
import math

from .kb import load_gene_expression, load_genes

# Karr knowledge-base mRNA/tRNA/rRNA half-life scheme (minutes), recovered from the
# KB grouped by RNA class. mRNAs are short-lived; stable RNAs are long-lived.
_HALFLIFE_MIN = {"mRNA": 2.5, "tRNA": 45.0, "rRNA": 1200.0, "sRNA": 45.0}
_CELL_CYCLE_S = 9.0 * 3600.0            # M. genitalium cycle ~9 h (dilution term)
# translation: protein copies scale with mRNA; keep the well-known abundant genes'
# translation efficiency and a default for the rest (protein half-life ~ stable).
_PROTEIN_HALFLIFE_S = 25000.0
# Overall synthesis-rate scale: the ensemble median mRNA synthesis rate that the
# reduced stochastic model was calibrated against (RNA t50 ≈ 18 min, Fig 3C). The
# ParCa preserves this median while taking the *relative* per-gene rates from the
# real observed expression, so kinetics stay in range but the distribution is real.
_MEDIAN_MRNA_SYNTH_PER_S = 5.5e-4


def rna_type(gene: dict) -> str:
    """Classify a gene's RNA product from its name/id (rRNA / tRNA / mRNA)."""
    text = f"{gene.get('name', '')} {gene.get('gene_id', '')}".lower()
    if "ribosomal rna" in text or "rrna" in text or "rrn" in text.replace("_", ""):
        return "rRNA"
    if "trna" in text:
        return "tRNA"
    return "mRNA"


@functools.lru_cache(maxsize=1)
def calculate_parameters() -> dict:
    """Run the ParCa initialization → per-gene parameter panel.

    Returns ``{key: (mrna_synthesis_rate_per_s, mrna_halflife_s,
    translation_rate_per_mrna_per_s, protein_halflife_s, avg_gene_length_nt)}``
    keyed by gene symbol (falling back to gene id) — the same shape the
    transcription / translation / decay processes read.
    """
    genes = load_genes()
    expr = load_gene_expression()
    if not expr:
        raise RuntimeError(
            "no observed gene expression available (datasets/karr_gene_expression.csv "
            "missing); run scripts/extract_kb_genes.py")

    def _rt(g):
        """Real KB RNA type when decoded, else name-based classification."""
        return (expr.get(g["gene_id"], {}).get("rna_type") or "").strip() or rna_type(g)

    def _half_s(g):
        """Real KB per-gene half-life (min → s), else the by-type scheme."""
        hl = expr.get(g["gene_id"], {}).get("half_life_min")
        if hl and hl > 0:
            return hl * 60.0
        return _HALFLIFE_MIN.get(_rt(g), _HALFLIFE_MIN["mRNA"]) * 60.0

    def _len(g):
        return int(expr.get(g["gene_id"], {}).get("length_nt") or 1000)

    # observed expression (37 C column, mean fallback); impute missing/zero genes
    # with half the smallest observed value (FitConstants imputeMissingData).
    obs = {}
    for g in genes:
        e = expr.get(g["gene_id"], {})
        v = e.get("expression_37C")
        if v is None or v <= 0:
            v = e.get("expression_mean")
        obs[g["gene_id"]] = v if (v and v > 0) else None
    positive = [v for v in obs.values() if v]
    floor = (min(positive) * 0.5) if positive else 1.0

    # relative synthesis-rate weight per gene ∝ expression × (dilution + decay)
    # (FitConstants "match expression, decay rates" balance), scaled so the median
    # mRNA gene lands on the calibrated ensemble median synthesis rate.
    weight = {}
    for g in genes:
        gid = g["gene_id"]
        e = obs[gid] if obs[gid] else floor
        decay = math.log(2.0) / _half_s(g)
        weight[gid] = e * (math.log(2.0) / _CELL_CYCLE_S + decay)

    mrna_weights = sorted(weight[g["gene_id"]] for g in genes if _rt(g) == "mRNA")
    med_w = mrna_weights[len(mrna_weights) // 2] if mrna_weights else 1.0
    scale = _MEDIAN_MRNA_SYNTH_PER_S / med_w if med_w else 1.0

    panel = {}
    for g in genes:
        gid = g["gene_id"]
        key = (g.get("symbol") or "").strip() or gid
        half_s = _half_s(g)
        synth = weight[gid] * scale
        # translation rate: proportional to synthesis (more mRNA → more protein),
        # in the same band the reduced model used (0.005–0.15 /mrna/s).
        transl = min(0.15, max(0.005, synth * 120.0))
        panel[key] = (synth, half_s, transl, _PROTEIN_HALFLIFE_S, _len(g))
    return panel


def synthesis_rates() -> dict:
    return {g: v[0] for g, v in calculate_parameters().items()}


def mrna_decay_rates() -> dict:
    return {g: math.log(2.0) / v[1] for g, v in calculate_parameters().items()}


def translation_rates() -> dict:
    return {g: v[2] for g, v in calculate_parameters().items()}
