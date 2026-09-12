"""RNA and protein decay submodels — clean-room reproduction.

Reproduces the Poisson-decay mechanism of the Karr 2012 ``RNADecay`` and
``ProteinDecay`` submodels: each RNA / protein species is degraded as a Poisson
process at its own rate (ln2 / half-life). The original also salvages
NMPs/amino acids and consumes protease/water; this reduced version keeps the
species-level Poisson decay that sets steady-state copy numbers (Fig 2G/2H) and
the balance against synthesis.
"""

from __future__ import annotations

import numpy as np
from process_bigraph import Process

from ..expression_defaults import mrna_decay_rates, protein_decay_rates


class _PoissonDecay(Process):
    """Shared Poisson species-decay Process. Subclasses set the port name and
    default rates."""

    _PORT = "counts"
    config_schema = {
        "decay_rates": {"_type": "map[float]", "_default": {}},
        "seed": {"_type": "integer", "_default": 2},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rates = dict(self.config["decay_rates"]) or self._default_rates()
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def _default_rates(self):
        return {}

    def inputs(self):
        return {self._PORT: "map[float]"}

    def outputs(self):
        return {self._PORT: "map[float]"}

    def initial_state(self):
        return {self._PORT: {}}

    def update(self, state, interval):
        counts = dict(state.get(self._PORT, {}) or {})
        deltas = {}
        for sp, n in counts.items():
            n = float(n)
            if n <= 0:
                continue
            rate = self._rates.get(sp)
            if rate is None:
                continue
            expected = n * rate * interval
            decayed = min(n, float(self._rng.poisson(max(expected, 0.0))))
            if decayed > 0:
                deltas[sp] = -decayed
        return {self._PORT: deltas}


class RnaDecayReproductionProcess(_PoissonDecay):
    """Poisson decay of mRNA species (reproduction of Karr 2012 RNADecay)."""

    description = (
        "Stochastic mRNA decay — reproduction of Karr 2012 RNADecay.\n"
        "Each species s degrades as a Poisson process at rate ln2 / half-life:\n"
        "    decayed_s ~ min(n_s, Poisson(n_s · λ_s · Δt))\n"
        "setting steady-state copy numbers together with transcription.\n"
        "Contract — in/out: rna_counts (per-gene mRNA; a negative-Δ map).\n"
        "Fidelity: REDUCED — species-level Poisson decay on the representative panel "
        "(NMP salvage/water accounting of the original omitted)."
    )
    _PORT = "rna_counts"

    def _default_rates(self):
        return mrna_decay_rates()

    def inputs(self):
        return {"rna_counts": "map[float]"}

    def outputs(self):
        return {"rna_counts": "map[float]"}

    def initial_state(self):
        return {"rna_counts": {}}


class ProteinDecayReproductionProcess(_PoissonDecay):
    """Poisson decay of protein species (reproduction of Karr 2012 ProteinDecay)."""

    description = (
        "Stochastic protein decay — reproduction of Karr 2012 ProteinDecay.\n"
        "Each species s degrades as a Poisson process at rate ln2 / half-life:\n"
        "    decayed_s ~ min(n_s, Poisson(n_s · λ_s · Δt))\n"
        "balancing translation to set steady-state protein copy numbers.\n"
        "Contract — in/out: protein_counts (per-gene protein; a negative-Δ map).\n"
        "Fidelity: REDUCED — species-level Poisson decay on the representative panel "
        "(protease/peptidase + ATP/amino-acid accounting of the original omitted)."
    )
    _PORT = "protein_counts"

    def _default_rates(self):
        return protein_decay_rates()

    def inputs(self):
        return {"protein_counts": "map[float]"}

    def outputs(self):
        return {"protein_counts": "map[float]"}

    def initial_state(self):
        return {"protein_counts": {}}
