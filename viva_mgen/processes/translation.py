"""Translation submodel — clean-room reproduction of Karr 2012.

Reproduces the stochastic protein-synthesis mechanism of the Karr 2012
``Translation`` submodel: ribosomes translate mRNAs into protein monomers,
consuming GTP (and, in full, aa-tRNAs). The original is a full ribosome state
machine with tmRNA stalling; this version keeps the essential mRNA-proportional
Poisson synthesis over the FULL M. genitalium gene set (~522 genes), with
per-gene rates fitted by the native ParCa, consuming GTP at ~2 per peptide bond
— sufficient for the burst dynamics of Fig 2G and the translation share of the
energy budget (Fig 5). The reduction is the ribosome state machine itself, not
the gene coverage.
"""

from __future__ import annotations

import numpy as np
from process_bigraph import Process

from ..expression_defaults import translation_rates
from .allocation import select_budget, demand_entry


class TranslationReproductionProcess(Process):
    """Stochastic protein synthesis proportional to mRNA copy number, full gene set.

    Inputs
    ------
    rna_counts : map[float]
        Current per-gene mRNA counts (translation propensity per gene).
    gtp : float
        Available GTP pool (molecules); ~2 GTP consumed per peptide elongation
        cap on synthesis.

    Outputs
    -------
    protein_counts : map[float]
        Per-gene protein count deltas (additive).
    gtp : float
        GTP consumed this interval (negative delta).
    """

    description = (
        "Stochastic translation — reproduction of Karr 2012 Translation.\n"
        "For every gene g in the full M. genitalium gene set (~522 genes), new protein per step is\n"
        "Poisson in the current mRNA copy number\n"
        "    n_g ~ Poisson(r_g · mRNA_g · Δt)\n"
        "capped by GTP supply (~2 GTP per peptide bond, ≈ gtp_per_protein per chain). Per-gene\n"
        "rates r_g are ParCa-fitted.\n"
        "Contract — in: rna_counts (per-gene mRNA copies), gtp (pool, molecules). "
        "out: protein_counts (per-gene protein Δ, additive map), gtp (Δ consumed, negative).\n"
        "Fidelity: FAITHFUL gene coverage (all ~522 genes) and ParCa-fitted per-gene rates. The\n"
        "reduction is the ribosome state machine (elongation/tmRNA stalling collapsed into one\n"
        "mRNA-proportional, GTP-limited propensity); drives the Fig 2G/2H mRNA↔protein decoupling\n"
        "and the translation share of the Fig 5 energy budget. Consumption is arbitrated by the\n"
        "whole-cell resource allocator (Karr hybrid partitioning): capped each tick at its allocated\n"
        "GTP budget from the finite metabolism-replenished pool."
    )

    config_schema = {
        "translation_rates": {"_type": "map[float]", "_default": {}},
        "gtp_per_protein": {"_type": "float", "_default": 600.0},  # ~2*avg aa length
        "seed": {"_type": "integer", "_default": 1},
        "consumer_id": {"_type": "string", "_default": "translation"},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rates = dict(self.config["translation_rates"]) or translation_rates()
        self._rng = np.random.default_rng(int(self.config["seed"]))
        self._cid = self.config["consumer_id"]

    def inputs(self):
        return {"rna_counts": "map[float]", "gtp": "float", "alloc__gtp": "map[float]"}

    def outputs(self):
        return {"protein_counts": "map[float]", "gtp": "float", "demand__gtp": "map[float]"}

    def initial_state(self):
        return {"rna_counts": {}, "gtp": 1e7}

    def update(self, state, interval):
        mrna = dict(state.get("rna_counts", {}) or {})
        gtp_avail = float(state.get("gtp", 0.0))
        cost = float(self.config["gtp_per_protein"])
        budget = select_budget(state.get("alloc__gtp", {}), self._cid)
        gtp_cap = min(gtp_avail, budget)  # gtp_avail still bounds as a floor safety
        # 1) draw each gene's would-be synthesis (Poisson in its mRNA copies) and
        #    the whole-process GTP demand, independent of the budget.
        drawn = {}
        want_gtp = 0.0
        for gene, rate in self._rates.items():
            copies = float(mrna.get(gene, 0.0))
            if copies <= 0:
                continue
            expected = rate * copies * interval
            want_gtp += expected * cost  # pre-clamp Poisson-expectation GTP demand
            n = int(self._rng.poisson(max(expected, 0.0)))
            if n > 0:
                drawn[gene] = n
        # 2) share the GTP budget PROPORTIONALLY across genes rather than
        #    first-come per dict order — otherwise whichever high-demand gene is
        #    visited first monopolizes the whole budget and starves the rest
        #    (Karr's translation spreads elongation across all mRNAs; no single
        #    gene consumes all ribosomes/GTP). When the draw fits in budget every
        #    gene translates in full; when it does not, each gene keeps the same
        #    fraction of its draw.
        total_need = sum(n * cost for n in drawn.values())
        scale = 1.0 if total_need <= gtp_cap else (gtp_cap / total_need if total_need > 0 else 0.0)
        new_counts = {}
        gtp_used = 0.0
        for gene, n in drawn.items():
            m = n if scale >= 1.0 else int(n * scale)
            if m <= 0:
                continue
            new_counts[gene] = float(m)
            gtp_used += m * cost
        return {"protein_counts": new_counts, "gtp": -gtp_used,
                "demand__gtp": demand_entry(self._cid, want_gtp)}
