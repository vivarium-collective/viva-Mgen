"""DNA replication submodel — clean-room reproduction (reduced mechanism).

Reproduces the cell-cycle-regulation mechanism the Karr 2012 model uncovered
(Fig 4): the *M. genitalium* cell cycle has three phases — replication
initiation, replication, and cytokinesis — and the model showed an *emergent*,
genetically-unregulated control of cell-cycle duration:

* Replication initiation proceeds by stochastic accumulation of DnaA into the
  oriC complex; once a threshold is reached, replication begins. More initial
  DnaA ⇒ shorter initiation (Fig 4C).
* During the (deterministic, dNTP-limited) replication phase, DNA polymerase
  consumes dNTPs. Crucially, dNTPs accumulate *during* the initiation phase
  (synthesized but not yet consumed), so a longer initiation builds a larger
  dNTP surplus, which makes the subsequent replication *faster* — an inverse
  relationship between initiation and replication durations (Fig 4D/4E) that
  buffers total cell-cycle length.

This is a reduced but mechanistically genuine reproduction of that emergent
regulation, not the full replisome/Okazaki-fragment machinery of the original
``Replication`` submodel.
"""

from __future__ import annotations

import numpy as np
from process_bigraph import Process

from .. import constants as C

_INIT, _REPL, _DONE = 0.0, 1.0, 2.0


class ReplicationReproductionProcess(Process):
    """dNTP-buffered three-phase replication with emergent duration control.

    Inputs
    ------
    dntp_synthesis_scale : float
        Multiplier (default 1.0) on dNTP synthesis rate — a sibling could
        couple metabolism to replication through it.

    Outputs (all current-value snapshots)
    -------
    replicated_fraction : overwrite[float]   fraction of the chromosome copied
    dntp_pool : overwrite[float]             current dNTP molecules
    dnaA_complex : overwrite[float]          DnaA monomers in the oriC complex
    phase_code : overwrite[float]            0=initiation, 1=replication, 2=done
    chromosome_copy : overwrite[float]       1 → 2 upon completion
    initiation_duration : overwrite[float]   seconds (0 until initiation ends)
    replication_duration : overwrite[float]  seconds (0 until replication ends)
    dntp_at_replication_start : overwrite[float]
    """

    description = (
        "dNTP-buffered three-phase replication — reduced reproduction of the emergent "
        "cell-cycle regulation Karr 2012 uncovered (Fig 4).\n"
        "INITIATION: DnaA accumulates stochastically (Poisson) into the oriC complex while "
        "dNTPs build up (at a reduced pre-S-phase rate) unconsumed; replication begins when "
        "DnaA ≥ threshold. REPLICATION: dNTP synthesis up-regulates to its full rate and DNA "
        "polymerase advances min(rate·Δt, dNTP) per step, consuming 1 dNTP/nt; done at terC "
        "after ~4 h (dNTP-synthesis-limited).\n"
        "Because dNTPs accumulate during initiation, a longer initiation builds a larger surplus "
        "that makes replication FASTER — an emergent inverse initiation↔replication relationship "
        "that buffers total cell-cycle length (Fig 4C/D/E).\n"
        "Contract — in: dntp_synthesis_scale (a sibling can couple metabolism→replication). "
        "out (snapshots): replicated_fraction, dntp_pool, dnaA_complex, phase_code "
        "(0=init,1=repl,2=done), chromosome_copy, initiation_duration, replication_duration, "
        "dntp_at_replication_start.\n"
        "Fidelity: REDUCED mechanism (genuine emergent-regulation dynamics; not the full "
        "replisome/Okazaki machinery)."
    )

    config_schema = {
        "genome_length_bp": {"_type": "float", "_default": float(C.GENOME_LENGTH_BP)},
        "dnaA_threshold": {"_type": "float", "_default": 30.0},
        "dnaA_synthesis_rate": {"_type": "float", "_default": 30.0 / 14400.0},  # → ~3.6 h mean initiation
        "initial_dnaA": {"_type": "float", "_default": 0.0},
        "dntp_synthesis_rate": {"_type": "float", "_default": 580070.0 / 15571.0},  # nt-equiv/s
        # dNTP synthesis runs SLOWER during initiation than during replication:
        # dNTP synthesis is up-regulated for S-phase, so only a modest surplus
        # accumulates before replication (a large surplus would let the pol race
        # through the genome in ~1 h instead of the ~4.33 h the paper reports).
        # This keeps replication dNTP-limited near 4.33 h while the surplus still
        # MODULATES it — preserving the emergent inverse initiation↔replication
        # relationship (Fig 4E) as cell-to-cell variation rather than dominating it.
        "init_synth_fraction": {"_type": "float", "_default": 0.08},
        "dna_pol_rate": {"_type": "float", "_default": 250.0},  # nt/s (both replisomes), pol cap
        "initial_dntp": {"_type": "float", "_default": 0.0},
        "seed": {"_type": "integer", "_default": 0},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rng = np.random.default_rng(int(self.config["seed"]))
        self._phase = _INIT
        self._pos = 0.0
        self._dnaA = float(self.config["initial_dnaA"])
        self._dntp = float(self.config["initial_dntp"])
        self._t = 0.0
        self._init_dur = 0.0
        self._repl_dur = 0.0
        self._dntp_at_start = 0.0

    def inputs(self):
        return {"dntp_synthesis_scale": "float"}

    def outputs(self):
        return {
            "replicated_fraction": "overwrite[float]",
            "dntp_pool": "overwrite[float]",
            "dnaA_complex": "overwrite[float]",
            "phase_code": "overwrite[float]",
            "chromosome_copy": "overwrite[float]",
            "initiation_duration": "overwrite[float]",
            "replication_duration": "overwrite[float]",
            "dntp_at_replication_start": "overwrite[float]",
        }

    def initial_state(self):
        return {"dntp_synthesis_scale": 1.0}

    def update(self, state, interval):
        scale = float(state.get("dntp_synthesis_scale", 1.0) or 1.0)
        synth = self.config["dntp_synthesis_rate"] * scale * interval
        self._t += interval
        genome = self.config["genome_length_bp"]

        if self._phase == _INIT:
            # DnaA accumulates stochastically; dNTP accumulates (no consumption yet)
            # but at the reduced pre-S-phase rate, so only a modest surplus builds.
            self._dnaA += float(self._rng.poisson(max(self.config["dnaA_synthesis_rate"] * interval, 0.0)))
            self._dntp += synth * self.config["init_synth_fraction"]
            if self._dnaA >= self.config["dnaA_threshold"]:
                self._phase = _REPL
                self._init_dur = self._t
                self._dntp_at_start = self._dntp
        elif self._phase == _REPL:
            # polymerase advances, capped by both pol rate and dNTP availability
            nt = min(self.config["dna_pol_rate"] * interval, self._dntp)
            self._pos += nt
            self._dntp = self._dntp - nt + synth
            if self._pos >= genome:
                self._pos = genome
                self._phase = _DONE
                self._repl_dur = self._t - self._init_dur

        return {
            "replicated_fraction": self._pos / genome,
            "dntp_pool": self._dntp,
            "dnaA_complex": self._dnaA,
            "phase_code": self._phase,
            "chromosome_copy": 2.0 if self._phase == _DONE else 1.0,
            "initiation_duration": self._init_dur,
            "replication_duration": self._repl_dur,
            "dntp_at_replication_start": self._dntp_at_start,
        }
