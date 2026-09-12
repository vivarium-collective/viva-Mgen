"""Transcription submodel — clean-room reproduction (reduced).

Reproduces the stochastic RNA-synthesis mechanism of the Karr 2012
``Transcription`` submodel: RNA polymerase binds transcription units and
synthesizes mRNA, consuming NTPs and releasing pyrophosphate. The original is a
full RNA-polymerase state machine over all transcription units; this reduced
version keeps the essential stochastic (Poisson) synthesis on a representative
gene panel (see :mod:`viva_mgen.expression_defaults`), sufficient for the
single-cell burst dynamics of Fig 2G/2H.
"""

from __future__ import annotations

import numpy as np
from process_bigraph import Process

from ..expression_defaults import synthesis_rates, gene_lengths


class TranscriptionReproductionProcess(Process):
    """Stochastic mRNA synthesis on a representative gene panel.

    Inputs
    ------
    ntp : float
        Available NTP pool (molecules). Synthesis is capped by NTP supply.
    rna_pol : float
        Available RNA-polymerase count (scales synthesis propensity).

    Outputs
    -------
    rna_counts : map[float]
        Per-gene mRNA count deltas (additive) — newly synthesized transcripts.
    ntp : float
        NTP consumed this interval (negative delta), composes with metabolism.
    """

    config_schema = {
        "synthesis_rates": {"_type": "map[float]", "_default": {}},
        "gene_lengths": {"_type": "map[float]", "_default": {}},
        "rna_pol_reference": {"_type": "float", "_default": 100.0},
        "seed": {"_type": "integer", "_default": 0},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rates = dict(self.config["synthesis_rates"]) or synthesis_rates()
        self._lengths = dict(self.config["gene_lengths"]) or gene_lengths()
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"ntp": "float", "rna_pol": "float"}

    def outputs(self):
        return {"rna_counts": "map[float]", "ntp": "float"}

    def initial_state(self):
        return {"ntp": 1e7, "rna_pol": self.config["rna_pol_reference"]}

    def update(self, state, interval):
        pol = float(state.get("rna_pol", self.config["rna_pol_reference"]))
        pol_factor = pol / float(self.config["rna_pol_reference"]) if self.config["rna_pol_reference"] else 1.0
        ntp_avail = float(state.get("ntp", 0.0))
        new_counts = {}
        ntp_used = 0.0
        for gene, rate in self._rates.items():
            expected = rate * pol_factor * interval
            n = int(self._rng.poisson(max(expected, 0.0)))
            if n <= 0:
                continue
            length = self._lengths.get(gene, 1200.0)
            need = n * length
            if ntp_used + need > ntp_avail:  # cap by NTP supply
                n = int(max(0, (ntp_avail - ntp_used) // length))
                need = n * length
            if n <= 0:
                continue
            new_counts[gene] = float(n)
            ntp_used += need
        return {"rna_counts": new_counts, "ntp": -ntp_used}
