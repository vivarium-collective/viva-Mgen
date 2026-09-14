"""Transcription submodel — clean-room reproduction of Karr 2012.

Reproduces the stochastic RNA-synthesis mechanism of the Karr 2012
``Transcription`` submodel: RNA polymerase binds transcription units and
synthesizes mRNA, consuming NTPs and releasing pyrophosphate. The original is a
full RNA-polymerase state machine over all transcription units; this version
keeps the essential stochastic (Poisson) single-molecule synthesis over the FULL
M. genitalium gene set (~522 genes, see :mod:`viva_mgen.expression_defaults`),
with per-gene rates fitted by the native ParCa from the real observed expression
profile — sufficient for the single-cell burst dynamics of Fig 2G/2H. The
reduction is the RNA-polymerase state machine itself (initiation/elongation/
termination are collapsed into one propensity), not the gene coverage.
"""

from __future__ import annotations

import numpy as np
from process_bigraph import Process

from ..expression_defaults import synthesis_rates, gene_lengths
from .allocation import select_budget, demand_entry


class TranscriptionReproductionProcess(Process):
    """Stochastic mRNA synthesis over the full M. genitalium gene set.

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

    description = (
        "Stochastic transcription — reproduction of Karr 2012 Transcription.\n"
        "For every gene g in the full M. genitalium gene set (~522 genes), new mRNA per step is a\n"
        "Poisson draw\n"
        "    n_g ~ Poisson(k_g · (RNApol / RNApol_ref) · Δt)\n"
        "capped by NTP supply (each transcript consumes length_g NTP; PPi released). Per-gene\n"
        "rates k_g are ParCa-fitted from the real observed expression profile.\n"
        "Contract — in: ntp (pool, molecules), rna_pol (available polymerase count). "
        "out: rna_counts (per-gene mRNA Δ, additive map), ntp (Δ consumed, negative).\n"
        "Fidelity: FAITHFUL gene coverage (all ~522 genes) and ParCa-fitted per-gene rates. The\n"
        "reduction is the RNA-polymerase state machine — initiation/elongation/termination are\n"
        "collapsed into a single Poisson propensity scaled by available polymerase; reproduces the\n"
        "bursty-mRNA behaviour of Fig 2G. Consumption is arbitrated by the whole-cell resource\n"
        "allocator (Karr hybrid partitioning): capped each tick at its allocated NTP budget from the\n"
        "finite metabolism-replenished pool."
    )

    config_schema = {
        "synthesis_rates": {"_type": "map[float]", "_default": {}},
        "gene_lengths": {"_type": "map[float]", "_default": {}},
        "rna_pol_reference": {"_type": "float", "_default": 100.0},
        "seed": {"_type": "integer", "_default": 0},
        "consumer_id": {"_type": "string", "_default": "transcription"},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rates = dict(self.config["synthesis_rates"]) or synthesis_rates()
        self._lengths = dict(self.config["gene_lengths"]) or gene_lengths()
        self._rng = np.random.default_rng(int(self.config["seed"]))
        self._cid = self.config["consumer_id"]

    def inputs(self):
        return {"ntp": "float", "rna_pol": "float", "alloc__ntp": "map[float]"}

    def outputs(self):
        return {"rna_counts": "map[float]", "ntp": "float", "demand__ntp": "map[float]"}

    def initial_state(self):
        return {"ntp": 1e7, "rna_pol": self.config["rna_pol_reference"]}

    def update(self, state, interval):
        pol = float(state.get("rna_pol", self.config["rna_pol_reference"]))
        pol_factor = pol / float(self.config["rna_pol_reference"]) if self.config["rna_pol_reference"] else 1.0
        ntp_avail = float(state.get("ntp", 0.0))
        budget = select_budget(state.get("alloc__ntp", {}), self._cid)
        ntp_cap = min(ntp_avail, budget)  # ntp_avail still bounds as a floor safety
        new_counts = {}
        ntp_used = 0.0
        want_ntp = 0.0
        for gene, rate in self._rates.items():
            expected = rate * pol_factor * interval
            length = self._lengths.get(gene, 1200.0)
            want_ntp += expected * length  # pre-clamp Poisson-expectation NTP demand
            n = int(self._rng.poisson(max(expected, 0.0)))
            if n <= 0:
                continue
            need = n * length
            if ntp_used + need > ntp_cap:  # cap by NTP supply/budget
                n = int(max(0, (ntp_cap - ntp_used) // length))
                need = n * length
            if n <= 0:
                continue
            new_counts[gene] = float(n)
            ntp_used += need
        return {"rna_counts": new_counts, "ntp": -ntp_used,
                "demand__ntp": demand_entry(self._cid, want_ntp)}
