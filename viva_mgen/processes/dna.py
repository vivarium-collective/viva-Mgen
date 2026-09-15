"""DNA / chromosome submodels — clean-room reproduction of Karr 2012.

Viva-native reproductions of the DNA-related submodels of the Karr 2012
*M. genitalium* whole-cell model. Each class reproduces the corresponding
original ``evolveState()`` mechanism with its real KB kinetic constants
(gyrase activity/ATP cost, DnaA cooperativity, repair energetics), and carries
the biology on a lumped state variable in place of the full per-molecule
chromosome bookkeeping.

Shared modeling choice, stated precisely per class: the Karr chromosome is a
``CircularSparseMat`` tracking every bound protein footprint, damaged site, and
per-region linking number. That per-site representation is a structural feature
of the original; these processes instead evolve the aggregate state it produces
(superhelical density σ, condensed/segregated fraction, lesion count, oriC
complex size). Each ``description`` names exactly what is faithful (mechanism +
real constants) and what the lumped state stands in for.

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
from ..chromosome_state import (N_CHROMOSOME_BINS, add_lesions, repair_sites, n_lesions,
                                N_SUPERCOIL_REGIONS, empty_linking_map, mean_sigma, relax_regions)
from .allocation import select_budget, demand_entry


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
        "Replication-initiation oriC assembly — reproduction of Karr 2012 "
        "ReplicationInitiation.\n"
        "DnaA is activated to DnaA-ATP (ATP-gated) and recruited to the oriC R1–R5 boxes, "
        "polymerizing COOPERATIVELY: the per-step recruitment scales with (1 + complex/threshold), "
        "so the already-bound DnaA promotes further binding. Dynamics:\n"
        "    recruited ~ min(dnaA_free, Poisson(k · (1 + complex/thr) · Δt))   [ATP>0 required]\n"
        "    complex += recruited;  initiation_ready = 1 when complex ≥ threshold.\n"
        "Contract — in: atp (activation gate), dnaA_free (recruitable monomers). "
        "out (snapshots): complex_size, initiation_ready (0/1).\n"
        "Fidelity: the cooperative, ATP-gated, threshold-fired assembly is faithful. The oriC "
        "complex is tracked as one cooperative DnaA-ATP pool rather than the KB's five R1–R5 boxes "
        "with per-box binding/unbinding + DnaA-ADP reactivation (parameters available via "
        "kb.karr_process_params('ReplicationInitiation'); a single cooperativity factor stands in "
        "for the KB siteCooperativity 170 / stateCooperativity 2 constants)."
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
        "DNA supercoiling homeostasis — reproduction of Karr 2012 DNASupercoiling.\n"
        "Gyrase introduces negative supercoils at its real KB activity rate (1.2 acts/gyrase/s, "
        "2 ATP → 2 supercoils per act) while topoisomerases relax; the net effect drives the "
        "superhelical density σ toward the maintained negative setpoint (KB −0.06). Dynamics "
        "(first-order relaxation, gyrase- and ATP-limited):\n"
        "    acts = min(gyrase·1.2·Δt, atp/2, supercoils_needed);  σ += sign(setpoint−σ)·acts/genome_turns;\n"
        "    σ → setpoint;  Δatp = −2·acts.\n"
        "Contract — in: gyrase (activity), atp (energy cap). "
        "out: superhelical_density (snapshot σ), atp (negative Δ).\n"
        "Fidelity: FAITHFUL constants (real KB gyraseActivityRate 1.2, gyraseATPCost 2.0, "
        "setpoint −0.06, 10.5 bp/turn). σ is now tracked PER-REGION on the shared chromosome "
        "structure (linking_number: {region -> σ} over N_SUPERCOIL_REGIONS topological regions, "
        "the viva-native stand-in for CircularSparseMat's per-region linking numbers), and the "
        "genome-wide superhelical_density observable is their mean; gyrase acts are distributed "
        "across regions, each relaxed toward the setpoint. Regions are homogeneous until "
        "replication/transcription wire per-region perturbation in (staged, FIDELITY_GAPS gap #3), "
        "so the mean equals the former single-pool σ (no behavior change). topoI/topoIV activity "
        "is lumped into the net gyrase relaxation (per-enzyme rates available via "
        "kb.karr_process_params('DNASupercoiling')). Consumption is arbitrated by the whole-cell "
        "resource allocator (Karr hybrid partitioning): capped each tick at its allocated ATP budget."
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
        "n_regions": {"_type": "integer", "_default": N_SUPERCOIL_REGIONS},
        "consumer_id": {"_type": "string", "_default": "supercoiling"},
        "seed": {"_type": "integer", "_default": 0},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._cid = self.config["consumer_id"]
        self._n_regions = int(self.config["n_regions"])
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"gyrase": "float", "atp": "float", "alloc__atp": "map[float]",
                "linking_number": "map[float]"}

    def outputs(self):
        return {"superhelical_density": "overwrite[float]", "atp": "float",
                "demand__atp": "map[float]", "linking_number": "overwrite[map[float]]"}

    def initial_state(self):
        return {"gyrase": 100.0, "atp": 1.0e6,
                "linking_number": empty_linking_map(self._n_regions, self.config["initial_sigma"])}

    def update(self, state, interval):
        gyrase = max(float(state.get("gyrase", 0.0) or 0.0), 0.0)
        atp = max(float(state.get("atp", 0.0) or 0.0), 0.0)
        setpoint = self.config["setpoint"]
        regions = dict(state.get("linking_number", {}) or {})
        if not regions:  # unwired/first tick fallback
            regions = empty_linking_map(self._n_regions, self.config["initial_sigma"])
        sigma = mean_sigma(regions)  # genome-wide σ = mean of per-region linking numbers
        # total supercoils the relaxed chromosome can hold, to normalize σ ↔ supercoil count
        turns = self.config["genome_length_bp"] / self.config["relaxed_bp_per_turn"]
        gap = setpoint - sigma  # how far below setpoint we still need to go (σ<0)
        # gyrase catalytic acts this step, capped by ATP (2 ATP/act) — MATLAB gyrase binding
        acts_wanted = abs(gap) * turns  # supercoils still needed
        want_acts = gyrase * self.config["gyrase_rate"] * interval
        want_atp = min(want_acts, acts_wanted) * self.config["atp_per_act"]
        budget = select_budget(state.get("alloc__atp", {}), self._cid)
        atp_cap = min(want_atp, budget, atp)  # atp still bounds as a floor safety
        acts = min(want_acts, atp_cap / self.config["atp_per_act"], acts_wanted)
        # distribute the acts across regions (one supercoil per turns_per_region per act);
        # the mean over regions moves by acts/turns exactly as the former single pool did.
        turns_per_region = turns / self._n_regions if self._n_regions else turns
        changed = relax_regions(regions, setpoint, acts, turns_per_region, self._rng)
        updated = dict(regions)
        updated.update(changed)
        atp_used = self.config["atp_per_act"] * acts
        return {"superhelical_density": mean_sigma(updated), "atp": -atp_used,
                "demand__atp": demand_entry(self._cid, want_atp),
                "linking_number": updated}


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
        "Chromosome condensation — reproduction of Karr 2012 ChromosomeCondensation.\n"
        "SMC complexes bind DNA at ~smcSepNt spacing (KB 7130 nt) and compact it by loop formation; "
        "binding is limited by available SMC (and, in the original, ATP/water). Dynamics — saturating "
        "approach to full compaction driven by SMC binding:\n"
        "    condensed_fraction += k · smc · (1 − condensed_fraction) · Δt   (clamped to [0,1]).\n"
        "Contract — in: smc (available complexes). out: condensed_fraction (snapshot 0→1).\n"
        "Fidelity: the SMC-limited saturating compaction is faithful. Compaction is carried as one "
        "genome-wide condensed fraction rather than per-site occupancy at smcSepNt spacing (7130 nt, "
        "via kb.karr_process_params('ChromosomeCondensation')); the SMC↔SMC-ADP ATPase cycle is "
        "delegated to the metabolite pools."
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
        "Chromosome segregation — reproduction of Karr 2012 ChromosomeSegregation.\n"
        "The original fires a single all-or-none segregation event once the chromosome is fully "
        "polymerized AND supercoiled (consuming GTP). Dynamics — gate on replication completion, "
        "then progress segregation to completion:\n"
        "    if replicated_fraction ≥ complete_threshold: segregated_fraction += rate · Δt  (→1).\n"
        "Contract — in: replicated_fraction (replication progress). "
        "out: segregated_fraction (snapshot 0→1).\n"
        "Fidelity: the replication-completion gate driving segregation to completion is faithful. "
        "Segregation is a continuous ramp rather than the KB's discrete decatenation event; the "
        "per-event GTP cost (KB gtpCost 1.0, via kb.karr_process_params('ChromosomeSegregation')) "
        "is delegated to the metabolite pools."
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
    lesion_map : map[float]  declared for spec/store-wiring symmetry with the
                             additive delta this process emits; not read here —
                             DNADamage only adds lesions, never consults existing ones.

    Outputs
    -------
    lesions : float          additive delta — new lesions introduced this step
    lesion_map : map[float]  additive per-site delta placing the new lesions at random bins
    """

    description = (
        "DNA damage — reproduction of Karr 2012 DNADamage.\n"
        "The original selects damage sites with probability stepSize · rate · [radiation/agent] per "
        "reaction. Dynamics — one lumped Poisson lesion source (spontaneous baseline plus an "
        "agent-driven term):\n"
        "    new_lesions ~ Poisson((base_rate + agent_rate · damaging_agent) · Δt).\n"
        "Contract — in: damaging_agent (radiation/reactive-species level, default 0). "
        "out: lesions (additive Δ, new lesions this step).\n"
        "Fidelity: the agent-driven Poisson damage law is faithful. Lesions are counted as one pool "
        "rather than typed and placed at vulnerable-motif sites; the per-reaction small-molecule "
        "reactant/product stoichiometry is delegated to the metabolite pools. Lesions are now tracked "
        "per-site on the shared chromosome structure (lesion_map) that damage adds to and repair "
        "clears from the same sites — the first consumer of the unified per-site chromosome "
        "(remaining DNA submodels staged). Constant provenance recorded in the gap-#5 audit "
        "(docs/CONSTANT_PROVENANCE.md)."
    )

    config_schema = {
        "base_rate": {"_type": "float", "_default": 1.0e-3},  # spontaneous lesions / s
        "agent_rate": {"_type": "float", "_default": 1.0e-2},  # lesions / s / unit agent
        "seed": {"_type": "integer", "_default": 3},
        "n_bins": {"_type": "integer", "_default": N_CHROMOSOME_BINS},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"damaging_agent": "float", "lesion_map": "map[float]"}

    def outputs(self):
        return {"lesions": "float", "lesion_map": "map[float]"}

    def initial_state(self):
        return {"damaging_agent": 0.0, "lesion_map": {}}

    def update(self, state, interval):
        agent = max(float(state.get("damaging_agent", 0.0) or 0.0), 0.0)
        # MATLAB: selectionProbability = stepSize · reactionBound · [radiation]; here a lumped Poisson.
        rate = self.config["base_rate"] + self.config["agent_rate"] * agent
        new_lesions = float(self._rng.poisson(max(rate * interval, 0.0)))
        return {"lesions": new_lesions,
                "lesion_map": add_lesions(self._rng, int(self.config["n_bins"]), int(new_lesions))}


# ---------------------------------------------------------------------------
# 6. DNA repair — excision/recombination repair of lesions, consuming ATP
# ---------------------------------------------------------------------------
class DNARepairReproductionProcess(Process):
    """Base/nucleotide-excision and recombination repair of lesions: enzymes
    repair up to rate·enzymes·Δt lesions per step, consuming ATP (polymerize +
    ligate).

    Inputs (current-value sensors)
    ------
    lesion_map : map[float]  per-site lesion counts — the source of truth for repair demand
    repair_enzyme : float    repair enzyme count (repair capacity)
    atp : float              ATP pool (limits repair)

    Outputs
    -------
    lesions : float   negative delta — lesions repaired this step
    atp : float       negative delta — ATP consumed
    """

    description = (
        "DNA repair — reproduction of Karr 2012 DNARepair.\n"
        "The original dispatches BER/NER/HR + polymerize + ligate subroutines over damaged sites, "
        "consuming ATP/dNTP. Dynamics — enzyme- and ATP-limited lesion clearance:\n"
        "    present = n_lesions(lesion_map);\n"
        "    repaired = min(present, rate · repair_enzyme · Δt, atp / atp_per_repair);\n"
        "    lesion_map -= repaired sites (weighted-random selection);  Δlesions = −repaired;  "
        "Δatp = −repaired · atp_per_repair.\n"
        "Contract — in: lesion_map (per-site counts, source of truth for repair demand), "
        "repair_enzyme (capacity), atp (energy cap). "
        "out: lesions (negative Δ), atp (negative Δ), lesion_map (negative per-site Δ).\n"
        "Fidelity: the enzyme- and ATP-limited repair flux is faithful. The distinct BER/NER/HR "
        "pathways, DisA scanning, and per-base polymerize/ligate steps are lumped into one clearance "
        "flux; dNTP accounting is delegated to the metabolite pools. Consumption is arbitrated by the "
        "whole-cell resource allocator (Karr hybrid partitioning): capped each tick at its allocated "
        "ATP budget from the finite metabolism-replenished pool. Lesions are now tracked per-site on "
        "the shared chromosome structure (lesion_map) that damage adds to and repair clears from the "
        "same sites — the first consumer of the unified per-site chromosome (remaining DNA submodels "
        "staged)."
    )

    config_schema = {
        "repair_rate": {"_type": "float", "_default": 1.0e-2},  # lesions / enzyme / s
        "atp_per_repair": {"_type": "float", "_default": 4.0},  # ATP per lesion (excise+polymerize+ligate)
        "consumer_id": {"_type": "string", "_default": "dna_repair"},
        "seed": {"_type": "integer", "_default": 5},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._cid = self.config["consumer_id"]
        self._rng = np.random.default_rng(int(self.config.get("seed", 5)))

    def inputs(self):
        return {"repair_enzyme": "float", "atp": "float",
                "alloc__atp": "map[float]", "lesion_map": "map[float]"}

    def outputs(self):
        return {"lesions": "float", "atp": "float", "demand__atp": "map[float]",
                "lesion_map": "map[float]"}

    def initial_state(self):
        return {"repair_enzyme": 50.0, "atp": 1.0e6, "lesion_map": {}}

    def update(self, state, interval):
        enzyme = max(float(state.get("repair_enzyme", 0.0) or 0.0), 0.0)
        atp = max(float(state.get("atp", 0.0) or 0.0), 0.0)
        atp_per = self.config["atp_per_repair"]
        lesion_map = state.get("lesion_map", {}) or {}
        # lesion_map is the sole source of truth for repair demand — the "lesions"
        # float store is an output-only observable (aggregate of damage/repair
        # deltas), never read back here.
        present = n_lesions(lesion_map)
        # MATLAB: repair subroutines clear damagedSites limited by enzymes + ATP/dNTP availability.
        capacity = self.config["repair_rate"] * enzyme * interval
        want_repaired = min(present, capacity)  # unconstrained-by-ATP desired repair
        want_atp = want_repaired * atp_per
        budget = select_budget(state.get("alloc__atp", {}), self._cid)
        atp_cap = min(want_atp, budget, atp)  # atp still bounds as a floor safety
        atp_cap_lesions = atp_cap / atp_per if atp_per > 0 else present
        capacity = min(present, capacity, atp_cap_lesions)
        delta, repaired = repair_sites(lesion_map, capacity, self._rng)
        return {"lesions": -repaired, "atp": -repaired * atp_per,
                "demand__atp": demand_entry(self._cid, want_atp),
                "lesion_map": delta}
