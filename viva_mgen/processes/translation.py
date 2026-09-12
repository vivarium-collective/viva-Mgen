"""Translation submodel — clean-room reproduction (reduced).

Reproduces the stochastic protein-synthesis mechanism of the Karr 2012
``Translation`` submodel: ribosomes translate mRNAs into protein monomers,
consuming GTP (and, in full, aa-tRNAs). The original is a full ribosome state
machine with tmRNA stalling; this reduced version keeps the essential
mRNA-proportional Poisson synthesis on the representative gene panel, consuming
GTP at ~2 per peptide bond, sufficient for the burst dynamics of Fig 2G and the
translation share of the energy budget (Fig 5).
"""

from __future__ import annotations

import numpy as np
from process_bigraph import Process

from ..expression_defaults import translation_rates


class TranslationReproductionProcess(Process):
    """Stochastic protein synthesis proportional to mRNA copy number.

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
        "Stochastic translation — reduced reproduction of Karr 2012 Translation.\n"
        "For each gene g, new protein per step is Poisson in the current mRNA copy number\n"
        "    n_g ~ Poisson(r_g · mRNA_g · Δt)\n"
        "capped by GTP supply (~2 GTP per peptide bond, ≈ gtp_per_protein per chain). Keeps the "
        "mRNA-proportional, GTP-limited synthesis of the full ribosome state machine.\n"
        "Contract — in: rna_counts (per-gene mRNA copies), gtp (pool, molecules). "
        "out: protein_counts (per-gene protein Δ, additive map), gtp (Δ consumed, negative).\n"
        "Fidelity: REDUCED — representative panel; drives the Fig 2G/2H mRNA↔protein decoupling "
        "and the translation share of the Fig 5 energy budget."
    )

    config_schema = {
        "translation_rates": {"_type": "map[float]", "_default": {}},
        "gtp_per_protein": {"_type": "float", "_default": 600.0},  # ~2*avg aa length
        "seed": {"_type": "integer", "_default": 1},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rates = dict(self.config["translation_rates"]) or translation_rates()
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"rna_counts": "map[float]", "gtp": "float"}

    def outputs(self):
        return {"protein_counts": "map[float]", "gtp": "float"}

    def initial_state(self):
        return {"rna_counts": {}, "gtp": 1e7}

    def update(self, state, interval):
        mrna = dict(state.get("rna_counts", {}) or {})
        gtp_avail = float(state.get("gtp", 0.0))
        cost = float(self.config["gtp_per_protein"])
        new_counts = {}
        gtp_used = 0.0
        for gene, rate in self._rates.items():
            copies = float(mrna.get(gene, 0.0))
            if copies <= 0:
                continue
            expected = rate * copies * interval
            n = int(self._rng.poisson(max(expected, 0.0)))
            if n <= 0:
                continue
            need = n * cost
            if gtp_used + need > gtp_avail:
                n = int(max(0, (gtp_avail - gtp_used) // cost))
                need = n * cost
            if n <= 0:
                continue
            new_counts[gene] = float(n)
            gtp_used += need
        return {"protein_counts": new_counts, "gtp": -gtp_used}
