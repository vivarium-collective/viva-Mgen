"""DNA / chromosome submodels — clean-room reproduction (reduced mechanisms).

Reduced but mechanistically genuine viva-native reproductions of the DNA-related
submodels of the Karr 2012 *M. genitalium* whole-cell model. Each class captures
the *essential* algorithm of the corresponding original ``evolveState()`` and
reduces the full chromosome/molecule bookkeeping (the CircularSparseMat
representation of bound proteins, damaged sites, linking numbers, per-base
footprints, water/H+/Pi balancing) to the state variable that carries the
biology.

These are written from the *described mechanism*, not by porting the MATLAB
line-for-line — a clean-room reproduction, not the original code.

Originals reproduced (``+edu/+stanford/+covert/+cell/+sim/+process/``):
  ReplicationInitiation.m, DNASupercoiling.m, ChromosomeCondensation.m,
  ChromosomeSegregation.m, DNADamage.m, DNARepair.m
"""

from __future__ import annotations

import numpy as np
from process_bigraph import Process

from .. import constants as C


# ---------------------------------------------------------------------------
# 1. Replication initiation — DnaA-ATP cooperative assembly of the oriC complex
# ---------------------------------------------------------------------------
class ReplicationInitiationReproductionProcess(Process):
    """DnaA-ATP monomers bind the oriC boxes and cooperatively build the oriC
    complex; when it reaches threshold, replication initiation fires.

    Inputs (current-value sensors)
    ------
    atp : float          ATP pool (gates activation of free DnaA → DnaA-ATP)
    dnaA_free : float    free DnaA monomers available to be recruited

    Outputs (current-value snapshots)
    -------
    complex_size : overwrite[float]      DnaA-ATP monomers in the oriC complex
    initiation_ready : overwrite[float]  1.0 once complex_size ≥ threshold, else 0.0
    """

    description = (
        "Replication-initiation oriC assembly — reduced reproduction of Karr 2012 "
        "ReplicationInitiation.\n"
        "DnaA is activated to DnaA-ATP (ATP-gated) and recruited to the oriC R1–R5 boxes, "
        "polymerizing COOPERATIVELY: the per-step recruitment scales with (1 + complex/threshold), "
        "so the already-bound DnaA promotes further binding. Reduced dynamics:\n"
        "    recruited ~ min(dnaA_free, Poisson(k · (1 + complex/thr) · Δt))   [ATP>0 required]\n"
        "    complex += recruited;  initiation_ready = 1 when complex ≥ threshold.\n"
        "Contract — in: atp (activation gate), dnaA_free (recruitable monomers). "
        "out (snapshots): complex_size, initiation_ready (0/1).\n"
        "Fidelity: REDUCED — cooperative single-pool assembly with a fire threshold; omits the "
        "per-box CircularSparseMat binding/unbinding, DnaA-ADP reactivation, and supercoiling gate."
    )

    config_schema = {
        "complex_threshold": {"_type": "float", "_default": 30.0},  # DnaA in complete oriC complex
        "recruit_rate": {"_type": "float", "_default": 30.0 / 12960.0},  # → ~init duration (Fig 4C)
        "initial_complex": {"_type": "float", "_default": 0.0},
        "seed": {"_type": "integer", "_default": 0},
        # The reduced cooperative recruitment uses (1 + complex/threshold) as its
        # cooperativity factor. The real Hill-type cooperativity and DnaA-ATP/ADP
        # binding constants — Karr 2012 parameters.json ReplicationInitiation
        # siteCooperativity 170, stateCooperativity 2, kb1ATP 25 / kb1ADP 2.5,
        # kb2ATP 0.61 / kb2ADP 0.61, kd1ATP 20 / kd1ADP 20, k_Regen 125,
        # K_Regen_P4 0.018 — parameterize the per-box binding/unbinding this single
        # cooperative pool omits, so they have no default slot here; available via
        # viva_mgen.kb.karr_process_params("ReplicationInitiation").
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rng = np.random.default_rng(int(self.config["seed"]))
        self._complex = float(self.config["initial_complex"])

    def inputs(self):
        return {"atp": "float", "dnaA_free": "float"}

    def outputs(self):
        return {
            "complex_size": "overwrite[float]",
            "initiation_ready": "overwrite[float]",
        }

    def initial_state(self):
        return {"atp": 1.0e6, "dnaA_free": 50.0}

    def update(self, state, interval):
        atp = float(state.get("atp", 0.0) or 0.0)
        dnaA_free = float(state.get("dnaA_free", 0.0) or 0.0)
        thr = self.config["complex_threshold"]
        # MATLAB: activateFreeDnaA + cooperative bindAndPolymerizeDnaAATP at oriC.
        # Cooperativity: bound DnaA promotes further recruitment (1 + complex/thr).
        if atp > 0.0 and self._complex < thr:
            coop = 1.0 + self._complex / thr
            expected = max(self.config["recruit_rate"] * coop * interval, 0.0)
            recruited = min(dnaA_free, float(self._rng.poisson(expected)))
            self._complex += recruited
        ready = 1.0 if self._complex >= thr else 0.0
        return {"complex_size": self._complex, "initiation_ready": ready}


# ---------------------------------------------------------------------------
# 2. DNA supercoiling — gyrase drives sigma toward the negative setpoint (ATP)
# ---------------------------------------------------------------------------
class DNASupercoilingReproductionProcess(Process):
    """Gyrase (and topoisomerases) set the superhelical density: gyrase introduces
    negative supercoils (~2 ATP → 2 supercoils per act), relaxing the density
    toward the maintained negative setpoint.

    Inputs (current-value sensors)
    ------
    gyrase : float   active gyrase count (rate of supercoiling activity)
    atp : float      ATP pool (limits gyrase activity)

    Outputs
    -------
    superhelical_density : overwrite[float]   sigma (current snapshot)
    atp : float                               negative delta (ATP consumed)
    """

    description = (
        "DNA supercoiling homeostasis — reduced reproduction of Karr 2012 DNASupercoiling.\n"
        "Gyrase introduces negative supercoils (~2 ATP → 2 supercoils per catalytic act) while "
        "topoisomerases relax; the net effect drives the superhelical density σ toward the "
        "maintained negative setpoint. Reduced dynamics (first-order relaxation, gyrase- and "
        "ATP-limited):\n"
        "    acts = min(gyrase·k·Δt, atp/2);  σ += acts/genome_turns · sign(setpoint − σ) …\n"
        "    σ → setpoint;  Δatp = −2·|Δsupercoils|.\n"
        "Contract — in: gyrase (activity), atp (energy cap). "
        "out: superhelical_density (snapshot σ), atp (negative Δ).\n"
        "Fidelity: REDUCED — lumped σ relaxation toward one setpoint; omits per-region "
        "linking-number tracking, topoI/topoIV/gyrase sigma-limit gating, and dwell-time binding."
    )

    config_schema = {
        "setpoint": {"_type": "float", "_default": -0.06},  # maintained negative superhelicity
        # supercoiling acts / gyrase / s = real gyrase activity rate
        # (Karr 2012 parameters.json DNASupercoiling.gyraseActivityRate = 1.2).
        # The reduced model lumps all topoisomerase activity into this one gyrase
        # relaxation; the real topoI/topoIV activity rates (topoIActivityRate 1.0,
        # topoIVActivityRate 2.5), per-act linking-number changes (gyraseDeltaLK -2,
        # topoIDeltaLK 1, topoIVDeltaLK -2), and gyraseMeanDwellTime 45 s have no
        # per-enzyme slot here — available via
        # viva_mgen.kb.karr_process_params("DNASupercoiling").
        "gyrase_rate": {"_type": "float", "_default": 1.2},
        # ATP per gyrase catalytic act (Karr 2012 parameters.json DNASupercoiling.gyraseATPCost = 2.0).
        "atp_per_act": {"_type": "float", "_default": 2.0},
        "relaxed_bp_per_turn": {"_type": "float", "_default": 10.5},
        "genome_length_bp": {"_type": "float", "_default": float(C.GENOME_LENGTH_BP)},
        "initial_sigma": {"_type": "float", "_default": 0.0},  # start relaxed, gyrase supercoils it
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._sigma = float(self.config["initial_sigma"])

    def inputs(self):
        return {"gyrase": "float", "atp": "float"}

    def outputs(self):
        return {"superhelical_density": "overwrite[float]", "atp": "float"}

    def initial_state(self):
        return {"gyrase": 100.0, "atp": 1.0e6}

    def update(self, state, interval):
        gyrase = max(float(state.get("gyrase", 0.0) or 0.0), 0.0)
        atp = max(float(state.get("atp", 0.0) or 0.0), 0.0)
        setpoint = self.config["setpoint"]
        # total supercoils the relaxed chromosome can hold, to normalize σ ↔ supercoil count
        turns = self.config["genome_length_bp"] / self.config["relaxed_bp_per_turn"]
        gap = setpoint - self._sigma  # how far below setpoint we still need to go (σ<0)
        # gyrase catalytic acts this step, capped by ATP (2 ATP/act) — MATLAB gyrase binding
        acts_wanted = abs(gap) * turns  # supercoils still needed
        acts = min(gyrase * self.config["gyrase_rate"] * interval,
                   atp / self.config["atp_per_act"],
                   acts_wanted)
        dsigma = np.sign(gap) * acts / turns
        self._sigma += dsigma
        atp_used = self.config["atp_per_act"] * acts
        return {"superhelical_density": self._sigma, "atp": -atp_used}


# ---------------------------------------------------------------------------
# 3. Chromosome condensation — SMC binding compacts DNA (loop formation)
# ---------------------------------------------------------------------------
class ChromosomeCondensationReproductionProcess(Process):
    """SMC complexes bind DNA and compact it by loop formation; the condensed
    fraction rises toward full compaction as SMCs bind (saturating).

    Inputs (current-value sensor)
    ------
    smc : float   available SMC complexes

    Outputs
    -------
    condensed_fraction : overwrite[float]   0 → 1 (snapshot)
    """

    description = (
        "Chromosome condensation — reduced reproduction of Karr 2012 ChromosomeCondensation.\n"
        "SMC complexes bind DNA at ~smcSepNt spacing and compact it by loop formation; binding is "
        "limited by available SMC (and, in the original, ATP/water). Reduced dynamics — saturating "
        "approach to full compaction driven by SMC binding:\n"
        "    condensed_fraction += k · smc · (1 − condensed_fraction) · Δt   (clamped to [0,1]).\n"
        "Contract — in: smc (available complexes). out: condensed_fraction (snapshot 0→1).\n"
        "Fidelity: REDUCED — lumped compaction fraction; omits per-site stochastic binding at "
        "smcSepNt spacing, the SMC↔SMC-ADP ATPase cycle, and ATP/water/H+/Pi balancing."
    )

    config_schema = {
        # per-SMC compaction rate (lumped). The real SMC binding-site geometry —
        # Karr 2012 parameters.json ChromosomeCondensation.smcSepNt 7130 (nt between
        # bound SMC complexes) and smcSepProbCenter 2800 — are spacings, not a rate,
        # so they have no slot in this saturating-fraction reduction; available via
        # viva_mgen.kb.karr_process_params("ChromosomeCondensation").
        "bind_rate": {"_type": "float", "_default": 1.0e-3},
        "initial_fraction": {"_type": "float", "_default": 0.0},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._frac = float(self.config["initial_fraction"])

    def inputs(self):
        return {"smc": "float"}

    def outputs(self):
        return {"condensed_fraction": "overwrite[float]"}

    def initial_state(self):
        return {"smc": 100.0}

    def update(self, state, interval):
        smc = max(float(state.get("smc", 0.0) or 0.0), 0.0)
        # MATLAB: bindProteinToChromosomeStochastically limited by available SMC → more bound = more
        # compacted. Saturating first-order approach: remaining uncondensed DNA compacts w/ SMC.
        self._frac += self.config["bind_rate"] * smc * (1.0 - self._frac) * interval
        self._frac = float(np.clip(self._frac, 0.0, 1.0))
        return {"condensed_fraction": self._frac}


# ---------------------------------------------------------------------------
# 4. Chromosome segregation — decatenate/segregate once replication completes
# ---------------------------------------------------------------------------
class ChromosomeSegregationReproductionProcess(Process):
    """Decatenation/segregation of the daughter chromosomes: proceeds only once
    replication is complete, after which the segregated fraction rises to 1.

    Inputs (current-value sensor)
    ------
    replicated_fraction : float   fraction of the chromosome copied (0 → 1)

    Outputs
    -------
    segregated_fraction : overwrite[float]   0 → 1 (snapshot)
    """

    description = (
        "Chromosome segregation — reduced reproduction of Karr 2012 ChromosomeSegregation.\n"
        "The original fires a single all-or-none segregation event once the chromosome is fully "
        "polymerized AND supercoiled (consuming GTP). Reduced dynamics — gate on replication "
        "completion, then progress segregation to completion:\n"
        "    if replicated_fraction ≥ complete_threshold: segregated_fraction += rate · Δt  (→1).\n"
        "Contract — in: replicated_fraction (replication progress). "
        "out: segregated_fraction (snapshot 0→1).\n"
        "Fidelity: REDUCED — replication-gated ramp; omits the supercoiled-state precondition, "
        "the GTP/water→GDP/Pi energetic cost, and the discrete decatenation event."
    )

    config_schema = {
        "complete_threshold": {"_type": "float", "_default": 0.999},  # replication ~complete
        "segregation_rate": {"_type": "float", "_default": 1.0 / 3869.0},  # ~cytokinesis duration
        "initial_fraction": {"_type": "float", "_default": 0.0},
        # The reduced ramp omits the segregation event's energetic cost; the real
        # per-event GTP cost (Karr 2012 parameters.json ChromosomeSegregation.gtpCost 1.0)
        # has no slot here — available via
        # viva_mgen.kb.karr_process_params("ChromosomeSegregation").
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._frac = float(self.config["initial_fraction"])

    def inputs(self):
        return {"replicated_fraction": "float"}

    def outputs(self):
        return {"segregated_fraction": "overwrite[float]"}

    def initial_state(self):
        return {"replicated_fraction": 0.0}

    def update(self, state, interval):
        replicated = float(state.get("replicated_fraction", 0.0) or 0.0)
        # MATLAB: segregate only when polymerizedRegions == full chromosome (replication done).
        if replicated >= self.config["complete_threshold"]:
            self._frac += self.config["segregation_rate"] * interval
            self._frac = float(np.clip(self._frac, 0.0, 1.0))
        return {"segregated_fraction": self._frac}


# ---------------------------------------------------------------------------
# 5. DNA damage — stochastic introduction of lesions (Poisson)
# ---------------------------------------------------------------------------
class DNADamageReproductionProcess(Process):
    """Stochastic introduction of DNA lesions (strand breaks, abasic sites). The
    per-step lesion count is Poisson, with rate raised by any damaging agent.

    Inputs (current-value sensor)
    ------
    damaging_agent : float   e.g. radiation / reactive-species level (default 0)

    Outputs
    -------
    lesions : float   additive delta — new lesions introduced this step
    """

    description = (
        "DNA damage — reduced reproduction of Karr 2012 DNADamage.\n"
        "The original selects damage sites with probability stepSize · rate · [radiation/agent] per "
        "reaction. Reduced dynamics — one lumped Poisson lesion source (spontaneous baseline plus an "
        "agent-driven term):\n"
        "    new_lesions ~ Poisson((base_rate + agent_rate · damaging_agent) · Δt).\n"
        "Contract — in: damaging_agent (radiation/reactive-species level, default 0). "
        "out: lesions (additive Δ, new lesions this step).\n"
        "Fidelity: REDUCED — single Poisson lesion pool; omits per-reaction damage types, "
        "vulnerable-motif site selection, and small-molecule reactant/product stoichiometry."
    )

    config_schema = {
        "base_rate": {"_type": "float", "_default": 1.0e-3},  # spontaneous lesions / s
        "agent_rate": {"_type": "float", "_default": 1.0e-2},  # lesions / s / unit agent
        "seed": {"_type": "integer", "_default": 3},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"damaging_agent": "float"}

    def outputs(self):
        return {"lesions": "float"}

    def initial_state(self):
        return {"damaging_agent": 0.0}

    def update(self, state, interval):
        agent = max(float(state.get("damaging_agent", 0.0) or 0.0), 0.0)
        # MATLAB: selectionProbability = stepSize · reactionBound · [radiation]; here a lumped Poisson.
        rate = self.config["base_rate"] + self.config["agent_rate"] * agent
        new_lesions = float(self._rng.poisson(max(rate * interval, 0.0)))
        return {"lesions": new_lesions}


# ---------------------------------------------------------------------------
# 6. DNA repair — excision/recombination repair of lesions, consuming ATP
# ---------------------------------------------------------------------------
class DNARepairReproductionProcess(Process):
    """Base/nucleotide-excision and recombination repair of lesions: enzymes
    repair up to rate·enzymes·Δt lesions per step, consuming ATP (polymerize +
    ligate).

    Inputs (current-value sensors)
    ------
    lesions : float          current lesion count
    repair_enzyme : float    repair enzyme count (repair capacity)
    atp : float              ATP pool (limits repair)

    Outputs
    -------
    lesions : float   negative delta — lesions repaired this step
    atp : float       negative delta — ATP consumed
    """

    description = (
        "DNA repair — reduced reproduction of Karr 2012 DNARepair.\n"
        "The original dispatches BER/NER/HR + polymerize + ligate subroutines over damaged sites, "
        "consuming ATP/dNTP. Reduced dynamics — enzyme- and ATP-limited lesion clearance:\n"
        "    repaired = min(lesions, rate · repair_enzyme · Δt, atp / atp_per_repair);\n"
        "    Δlesions = −repaired;  Δatp = −repaired · atp_per_repair.\n"
        "Contract — in: lesions (current count), repair_enzyme (capacity), atp (energy cap). "
        "out: lesions (negative Δ), atp (negative Δ).\n"
        "Fidelity: REDUCED — single repair flux with an ATP cost; omits the distinct BER/NER/HR "
        "pathways, DisA scanning, per-base polymerize/ligate steps, and dNTP accounting."
    )

    config_schema = {
        "repair_rate": {"_type": "float", "_default": 1.0e-2},  # lesions / enzyme / s
        "atp_per_repair": {"_type": "float", "_default": 4.0},  # ATP per lesion (excise+polymerize+ligate)
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)

    def inputs(self):
        return {"lesions": "float", "repair_enzyme": "float", "atp": "float"}

    def outputs(self):
        return {"lesions": "float", "atp": "float"}

    def initial_state(self):
        return {"lesions": 0.0, "repair_enzyme": 50.0, "atp": 1.0e6}

    def update(self, state, interval):
        lesions = max(float(state.get("lesions", 0.0) or 0.0), 0.0)
        enzyme = max(float(state.get("repair_enzyme", 0.0) or 0.0), 0.0)
        atp = max(float(state.get("atp", 0.0) or 0.0), 0.0)
        atp_per = self.config["atp_per_repair"]
        # MATLAB: repair subroutines clear damagedSites limited by enzymes + ATP/dNTP availability.
        capacity = self.config["repair_rate"] * enzyme * interval
        repaired = min(lesions, capacity, atp / atp_per if atp_per > 0 else lesions)
        return {"lesions": -repaired, "atp": -repaired * atp_per}
