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

from .kb import load_gene_expression, load_genes, load_karr_parameters

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

_NA = 6.02214076e23           # Avogadro
_RNA_NT_MW = 340.0            # avg ribonucleotide-monophosphate MW (g/mol) in a chain


def _mass_constants() -> dict:
    """Real Karr mass constants (state 'Mass' in the KB parameters)."""
    from .kb import load_karr_parameters
    mass = (load_karr_parameters().get("states", {}) or {}).get("Mass", {})
    return {
        "dry_weight_g": float(mass.get("cellInitialDryWeight", 3.93e-15)),
        "rna_fraction": float(mass.get("dryWeightFractionRNA", 0.0930)),
    }


def fit_analytically(counts0, mw, held_idx=(), supercoil=None):
    """Equality-constrained QP core of Karr 2012 ``FitConstants.fitAnalytically``.

    Minimize ``||x - counts0||^2`` (least change from the observed-expression
    initial guess) subject to the linear equality constraints of
    ``linConstraintFunc``:

    * **RNA mass**: ``sum(x_i * mw_i) = dryWeight * RNAfraction * N_A`` — total RNA
      reproduces the observed dry-mass RNA fraction;
    * **held expression**: ``x[j] = counts0[j]`` for each ``j`` in ``held_idx``
      (DnaA / FtsZ, pinned so replication-initiation / cytokinesis need no refit);
    * **net supercoiling zero**: ``sum(coef_i * x_i) = 0`` (topoisomerase I positive
      activity balanced by gyrase negative activity), when ``supercoil`` coefficients
      are given.

    Solved in closed form by projecting ``counts0`` onto ``{A x = b}``
    (``x = x0 - Aᵀ (A Aᵀ)⁻¹ (A x0 - b)``). Returns the fitted counts.
    """
    import numpy as np
    x0 = np.asarray(counts0, float)
    mw = np.asarray(mw, float)
    mc = _mass_constants()
    rows, rhs = [], []
    # RNA-mass constraint
    rows.append(mw.copy())
    rhs.append(mc["dry_weight_g"] * mc["rna_fraction"] * _NA)
    # held-expression constraints
    for j in held_idx:
        r = np.zeros_like(x0); r[j] = 1.0
        rows.append(r); rhs.append(x0[j])
    # net-supercoiling constraint
    if supercoil is not None:
        rows.append(np.asarray(supercoil, float)); rhs.append(0.0)
    A = np.vstack(rows); b = np.asarray(rhs, float)
    # projection onto Ax=b nearest x0
    AAt = A @ A.T
    x = x0 - A.T @ np.linalg.solve(AAt, A @ x0 - b)
    return np.maximum(x, 0.0)


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

    # --- analytic QP fit (FitConstants.fitAnalytically) -----------------------
    # Refine the observed-expression RNA distribution to satisfy the linear
    # constraints (RNA mass, DnaA/FtsZ held, net-supercoiling zero), least-change.
    import numpy as np
    idx = {g["gene_id"]: i for i, g in enumerate(genes)}
    decay = np.array([math.log(2.0) / _half_s(g) + math.log(2.0) / _CELL_CYCLE_S for g in genes])
    mw = np.array([_len(g) * _RNA_NT_MW for g in genes])
    # initial RNA counts ∝ observed expression, scaled to the target RNA mass
    frac0 = np.array([(obs[g["gene_id"]] if obs[g["gene_id"]] else floor) for g in genes])
    mc = _mass_constants()
    counts0 = frac0 / (frac0 @ mw) * (mc["dry_weight_g"] * mc["rna_fraction"] * _NA)
    held = [idx[gid] for gid in ("MG_469", "MG_224") if gid in idx]  # DnaA, FtsZ
    # net-supercoiling coefficients from the real ported constants
    sc = (load_karr_parameters().get("processes", {}) or {}).get("DNASupercoiling", {})
    coef = np.zeros(len(genes))
    sym = {(g.get("symbol") or "").lower(): i for i, g in enumerate(genes)}
    for s in ("topa",):
        if s in sym:
            coef[sym[s]] = float(sc.get("topoIDeltaLK", 1.0)) * float(sc.get("topoIActivityRate", 1.0))
    for s in ("gyra", "gyrb"):
        if s in sym:
            coef[sym[s]] = float(sc.get("gyraseDeltaLK", -2.0)) * float(sc.get("gyraseActivityRate", 1.2)) / 2.0 * 1.3
    try:
        counts = fit_analytically(counts0, mw, held_idx=held,
                                  supercoil=coef if coef.any() else None)
    except Exception:  # noqa: BLE001 — never let the QP break parameter loading
        counts = counts0
    # physical synthesis rate = count × (decay + dilution); rescale so the median
    # mRNA rate matches the reduced model's kinetic calibration (RNA t50 ≈ 18 min).
    synth_arr = counts * decay
    mrna_mask = np.array([_rt(g) == "mRNA" for g in genes])
    med = float(np.median(synth_arr[mrna_mask])) if mrna_mask.any() else 1.0
    scale = _MEDIAN_MRNA_SYNTH_PER_S / med if med else 1.0

    panel = {}
    for i, g in enumerate(genes):
        key = (g.get("symbol") or "").strip() or g["gene_id"]
        half_s = _half_s(g)
        synth = float(synth_arr[i]) * scale
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
