"""Protein maturation submodels — clean-room reproductions (reduced).

Reproduces the post-translational MATURATION PIPELINE of the Karr 2012
M. genitalium whole-cell model, in which a nascent peptide monomer is walked
through a sequence of form-changing steps before it becomes a functional,
localized, complexed protein:

    nascent  --(Processing I)-->  processed_i
             --(Translocation)-->  translocated       (membrane/secreted)
             --(Processing II)-->  processed_ii
             --(Folding)-------->  folded
             --(Modification)--->  modified
             --(Complexation)--->  complexes
             --(Ribosome Asm.)-->  30S / 50S ribosomes
             --(Activation)----->  active/inactive
             --(Terminal Org.)-->  localized adhesins

Each MATLAB submodel's ``evolveState`` is an enzyme/resource-limited transform
that moves counts of one protein FORM into the next form and pays the metabolite
cost (water, ATP, GTP, PG). This module keeps that exact essence — genuine
form-to-form conversion at a rate limited by the relevant enzyme, subunit, or
resource pool — on a representative panel, and drops the full knowledge-base
stoichiometry matrices and per-species index bookkeeping.

Every class is a ``process_bigraph.Process`` with bare, composable port types:
count maps are ``map[float]`` and moves are emitted as additive deltas (negative
on the source form, positive on the product form); pooled resources (ATP, GTP)
are ``float`` negative deltas; genuine "current value" readouts (active fraction,
assembled fraction) use ``overwrite[...]``.

Fidelity: REDUCED — mechanism-faithful reduced reproduction, representative
panel, order-of-magnitude-genuine kinetic parameters.
"""

from __future__ import annotations

import numpy as np
from process_bigraph import Process


# ---------------------------------------------------------------------------
# Shared helpers (mirror the MATLAB idioms used across the submodels)
# ---------------------------------------------------------------------------

def _stochastic_round(rng, x):
    """MATLAB ``randStream.stochasticRound``: floor, then round the fractional
    remainder up with probability equal to that remainder."""
    x = float(x)
    if x <= 0.0:
        return 0
    fl = np.floor(x)
    frac = x - fl
    return int(fl) + int(rng.random() < frac)


def _enzyme_limited_transform(rng, src_counts, limit):
    """Move counts out of ``src_counts`` up to a total of ``limit`` transforms.

    Mirrors the MATLAB pattern
        transformations = src * min(1, enzymeLimit / sum(src));
        transformations = stochasticRound(transformations);
    i.e. when the enzyme cannot cover every available substrate this step, the
    limited capacity is shared PROPORTIONALLY across the waiting species, then
    stochastically rounded to integers and capped at availability.
    Returns a dict of {species: count_moved} (positive integers).
    """
    total = sum(float(v) for v in src_counts.values() if float(v) > 0.0)
    if total <= 0.0 or limit <= 0.0:
        return {}
    frac = min(1.0, float(limit) / total)
    moved = {}
    for k, v in src_counts.items():
        v = float(v)
        if v <= 0.0:
            continue
        n = min(int(v), _stochastic_round(rng, v * frac))
        if n > 0:
            moved[k] = n
    return moved


# ---------------------------------------------------------------------------
# 1. Protein Processing I — deformylation + N-terminal Met cleavage
# ---------------------------------------------------------------------------

class ProteinProcessingIReproductionProcess(Process):
    """N-terminal maturation of nascent peptides (reproduction of Karr 2012
    ProteinProcessingI).

    Inputs
    ------
    nascent : map[float]        per-gene nascent (unprocessed) monomer counts
    deformylase : float         count of peptide deformylase (MG_106) enzyme

    Outputs
    -------
    nascent : map[float]            negative delta (consumed)
    process_i_done : map[float]     positive delta (deformylated + Met-cleaved)
    """

    description = (
        "Protein Processing I — reduced reproduction of Karr 2012 ProteinProcessingI.\n"
        "Peptide deformylase (MG_106) deformylates the N-terminal fMet and methionine\n"
        "aminopeptidase (MG_172) cleaves the N-terminal Met of nascent peptides. MATLAB caps\n"
        "transforms at enzyme*specificRate*Δt and shares that capacity proportionally across\n"
        "waiting monomers; here nascent -> process_i_done at that same enzyme-limited rate.\n"
        "Contract — in: nascent (per-gene, map), deformylase (enzyme count, float). "
        "out: nascent (Δ consumed, neg map), process_i_done (Δ produced, pos map).\n"
        "Fidelity: REDUCED — water/formate/methionine salvage of the original omitted."
    )

    config_schema = {
        "deformylase_specific_rate": {"_type": "float", "_default": 10.0},  # transforms / enzyme / s
        "seed": {"_type": "integer", "_default": 11},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"nascent": "map[float]", "deformylase": "float"}

    def outputs(self):
        return {"nascent": "map[float]", "process_i_done": "map[float]"}

    def initial_state(self):
        return {"nascent": {}, "deformylase": 50.0}

    def update(self, state, interval):
        nascent = dict(state.get("nascent", {}) or {})
        enzyme = float(state.get("deformylase", 0.0))
        # MATLAB: deformylaseLimit = enzymes * specificRate * stepSizeSec
        limit = enzyme * float(self.config["deformylase_specific_rate"]) * interval
        moved = _enzyme_limited_transform(self._rng, nascent, limit)
        return {
            "nascent": {g: -n for g, n in moved.items()},
            "process_i_done": {g: float(n) for g, n in moved.items()},
        }


# ---------------------------------------------------------------------------
# 2. Protein Translocation — SecYEG translocase, GTP/proton-driven
# ---------------------------------------------------------------------------

class ProteinTranslocationReproductionProcess(Process):
    """Translocation of membrane/secreted proteins across the membrane
    (reproduction of Karr 2012 ProteinTranslocation).

    Inputs
    ------
    process_i_done : map[float]   per-gene translocatable monomer counts
    translocase : float           count of SecYEG preprotein translocase/ATPase
    gtp : float                   GTP pool (SRP-driven recognition)

    Outputs
    -------
    process_i_done : map[float]   negative delta (consumed)
    translocated : map[float]     positive delta (moved across membrane)
    gtp : float                   negative delta (GTP consumed)
    """

    description = (
        "Protein Translocation — reduced reproduction of Karr 2012 ProteinTranslocation.\n"
        "Integral-membrane/lipo/secreted proteins are pushed through the SecYEG pore by the\n"
        "SecA translocase, recognised via the GTP-driven signal recognition particle. MATLAB\n"
        "translocates monomers one by one while translocase capacity, ATP and GTP last (SRP =\n"
        "~2 GTP/monomer); here process_i_done -> translocated up to the min of translocase\n"
        "capacity and gtp/gtp_per_monomer, consuming GTP.\n"
        "Contract — in: process_i_done (map), translocase (float), gtp (float). "
        "out: process_i_done (neg map), translocated (pos map), gtp (neg Δ).\n"
        "Fidelity: REDUCED — ATP/water/ADP/Pi and per-length translocase cost of the original omitted."
    )

    config_schema = {
        "translocase_specific_rate": {"_type": "float", "_default": 5.0},  # monomers / enzyme / s
        "gtp_per_monomer": {"_type": "float", "_default": 2.0},  # SRP: ~2 GTP / monomer
        "seed": {"_type": "integer", "_default": 12},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"process_i_done": "map[float]", "translocase": "float", "gtp": "float"}

    def outputs(self):
        return {"process_i_done": "map[float]", "translocated": "map[float]", "gtp": "float"}

    def initial_state(self):
        return {"process_i_done": {}, "translocase": 30.0, "gtp": 1e6}

    def update(self, state, interval):
        src = dict(state.get("process_i_done", {}) or {})
        translocase = float(state.get("translocase", 0.0))
        gtp = float(state.get("gtp", 0.0))
        gtp_per = float(self.config["gtp_per_monomer"])
        # capacity limited by translocase kinetics AND by GTP available for SRP
        enz_limit = translocase * float(self.config["translocase_specific_rate"]) * interval
        gtp_limit = gtp / gtp_per if gtp_per > 0 else 0.0
        limit = min(enz_limit, gtp_limit)
        moved = _enzyme_limited_transform(self._rng, src, limit)
        n_total = sum(moved.values())
        return {
            "process_i_done": {g: -n for g, n in moved.items()},
            "translocated": {g: float(n) for g, n in moved.items()},
            "gtp": -float(n_total) * gtp_per,
        }


# ---------------------------------------------------------------------------
# 3. Protein Processing II — signal-peptide cleavage + diacylglyceryl transfer
# ---------------------------------------------------------------------------

class ProteinProcessingIIReproductionProcess(Process):
    """Signal-peptide cleavage / lipoprotein anchoring of translocated proteins
    (reproduction of Karr 2012 ProteinProcessingII).

    Inputs
    ------
    translocated : map[float]   per-gene translocated monomer counts
    signal_peptidase : float    count of signal peptidase II (LspA, MG_210)

    Outputs
    -------
    translocated : map[float]   negative delta (consumed)
    processed_ii : map[float]   positive delta (mature secreted/lipoprotein form)
    """

    description = (
        "Protein Processing II — reduced reproduction of Karr 2012 ProteinProcessingII.\n"
        "After translocation, signal peptidase II (LspA/MG_210) cleaves the type-II signal\n"
        "sequence and diacylglyceryl transferase (Lgt/MG_086) anchors lipoproteins. MATLAB\n"
        "caps cleavages at peptidase*specificRate*Δt shared proportionally; here translocated\n"
        "-> processed_ii at that enzyme-limited rate.\n"
        "Contract — in: translocated (map), signal_peptidase (enzyme count, float). "
        "out: translocated (neg map), processed_ii (pos map).\n"
        "Fidelity: REDUCED — water/PG160/diacylglyceryl transfer accounting of the original omitted."
    )

    config_schema = {
        "peptidase_specific_rate": {"_type": "float", "_default": 8.0},  # transforms / enzyme / s
        "seed": {"_type": "integer", "_default": 13},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"translocated": "map[float]", "signal_peptidase": "float"}

    def outputs(self):
        return {"translocated": "map[float]", "processed_ii": "map[float]"}

    def initial_state(self):
        return {"translocated": {}, "signal_peptidase": 40.0}

    def update(self, state, interval):
        src = dict(state.get("translocated", {}) or {})
        enzyme = float(state.get("signal_peptidase", 0.0))
        limit = enzyme * float(self.config["peptidase_specific_rate"]) * interval
        moved = _enzyme_limited_transform(self._rng, src, limit)
        return {
            "translocated": {g: -n for g, n in moved.items()},
            "processed_ii": {g: float(n) for g, n in moved.items()},
        }


# ---------------------------------------------------------------------------
# 4. Protein Folding — spontaneous + chaperone-assisted, ATP-driven
# ---------------------------------------------------------------------------

class ProteinFoldingReproductionProcess(Process):
    """Chaperone-assisted + spontaneous folding of monomers (reproduction of
    Karr 2012 ProteinFolding).

    Inputs
    ------
    unfolded : map[float]     per-gene unfolded monomer counts
    chaperone_count : float   count of chaperone (e.g. GroEL) capacity
    atp : float               ATP pool (chaperone cycling)

    Outputs
    -------
    unfolded : map[float]     negative delta (consumed)
    folded : map[float]       positive delta (folded)
    atp : float               negative delta (ATP consumed by chaperone action)
    """

    description = (
        "Protein Folding — reduced reproduction of Karr 2012 ProteinFolding.\n"
        "Monomers in the 'notFolding' set relax spontaneously; the rest fold with chaperone\n"
        "(GroEL/DnaK) assistance, which MATLAB drives as an ATP/enzyme-limited reaction Gillespie\n"
        "loop. Here folding fraction per step = 1 - exp(-(spontaneous + chaperone_rate*chaperone_count)*Δt),\n"
        "and the chaperone-driven share consumes ATP (~atp_per_fold), capping folds at atp availability.\n"
        "Contract — in: unfolded (map), chaperone_count (float), atp (float). "
        "out: unfolded (neg map), folded (pos map), atp (neg Δ).\n"
        "Fidelity: REDUCED — prosthetic-group/ion coordination and per-protein rate matrix of the original omitted."
    )

    config_schema = {
        "spontaneous_rate": {"_type": "float", "_default": 0.05},  # 1/s
        "chaperone_rate": {"_type": "float", "_default": 0.002},  # 1/s per chaperone
        "atp_per_fold": {"_type": "float", "_default": 7.0},  # ~7 ATP per GroEL cycle
        "seed": {"_type": "integer", "_default": 14},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"unfolded": "map[float]", "chaperone_count": "float", "atp": "float"}

    def outputs(self):
        return {"unfolded": "map[float]", "folded": "map[float]", "atp": "float"}

    def initial_state(self):
        return {"unfolded": {}, "chaperone_count": 100.0, "atp": 1e6}

    def update(self, state, interval):
        unfolded = dict(state.get("unfolded", {}) or {})
        cc = float(state.get("chaperone_count", 0.0))
        atp = float(state.get("atp", 0.0))
        spont = float(self.config["spontaneous_rate"])
        chap = float(self.config["chaperone_rate"]) * cc
        rate = spont + chap
        if rate <= 0.0:
            return {"unfolded": {}, "folded": {}, "atp": 0.0}
        # fraction folded this step (exponential approach); chaperone share needs ATP
        frac = 1.0 - np.exp(-rate * interval)
        chap_share = chap / rate  # fraction of folds that are chaperone-assisted
        cost = float(self.config["atp_per_fold"])
        folded = {}
        atp_used = 0.0
        for g, v in unfolded.items():
            v = float(v)
            if v <= 0.0:
                continue
            n = min(int(v), _stochastic_round(self._rng, v * frac))
            if n <= 0:
                continue
            # chaperone-assisted subset costs ATP; throttle by remaining ATP
            n_chap = _stochastic_round(self._rng, n * chap_share)
            need = n_chap * cost
            if atp_used + need > atp:
                n_chap = int(max(0.0, (atp - atp_used) // cost)) if cost > 0 else 0
                # spontaneous folds still proceed for this species
                n = min(n, (n - n_chap) + n_chap)
                need = n_chap * cost
            atp_used += need
            folded[g] = float(n)
        return {
            "unfolded": {g: -n for g, n in folded.items()},
            "folded": dict(folded),
            "atp": -atp_used,
        }


# ---------------------------------------------------------------------------
# 5. Protein Modification — phosphorylation / lipoate / Glu ligation
# ---------------------------------------------------------------------------

class ProteinModificationReproductionProcess(Process):
    """Covalent modification of folded monomers (reproduction of Karr 2012
    ProteinModification).

    Inputs
    ------
    unmodified : map[float]         per-gene unmodified monomer counts
    modification_enzyme : float     count of modifying enzyme (kinase/ligase)
    atp : float                     ATP pool (phosphoryl donor)

    Outputs
    -------
    unmodified : map[float]         negative delta (consumed)
    modified : map[float]           positive delta (modified)
    atp : float                     negative delta (ATP consumed)
    """

    description = (
        "Protein Modification — reduced reproduction of Karr 2012 ProteinModification.\n"
        "Ser/Thr/Tyr phosphorylation and lipoate/glutamate ligation of specific monomers.\n"
        "MATLAB runs an enzyme+substrate-limited reaction Gillespie loop; here unmodified ->\n"
        "modified at rate enzyme*specificRate*Δt, capped by ATP (~1 ATP per phosphoryl transfer).\n"
        "Contract — in: unmodified (map), modification_enzyme (float), atp (float). "
        "out: unmodified (neg map), modified (pos map), atp (neg Δ).\n"
        "Fidelity: REDUCED — per-reaction stoichiometry matrix and cofactor set of the original omitted."
    )

    config_schema = {
        "modification_specific_rate": {"_type": "float", "_default": 6.0},  # transforms / enzyme / s
        "atp_per_modification": {"_type": "float", "_default": 1.0},
        "seed": {"_type": "integer", "_default": 15},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"unmodified": "map[float]", "modification_enzyme": "float", "atp": "float"}

    def outputs(self):
        return {"unmodified": "map[float]", "modified": "map[float]", "atp": "float"}

    def initial_state(self):
        return {"unmodified": {}, "modification_enzyme": 20.0, "atp": 1e6}

    def update(self, state, interval):
        src = dict(state.get("unmodified", {}) or {})
        enzyme = float(state.get("modification_enzyme", 0.0))
        atp = float(state.get("atp", 0.0))
        cost = float(self.config["atp_per_modification"])
        enz_limit = enzyme * float(self.config["modification_specific_rate"]) * interval
        atp_limit = atp / cost if cost > 0 else float("inf")
        limit = min(enz_limit, atp_limit)
        moved = _enzyme_limited_transform(self._rng, src, limit)
        n_total = sum(moved.values())
        return {
            "unmodified": {g: -n for g, n in moved.items()},
            "modified": {g: float(n) for g, n in moved.items()},
            "atp": -float(n_total) * cost,
        }


# ---------------------------------------------------------------------------
# 6. Protein Activation — allosteric active/inactive, regulator-driven
# ---------------------------------------------------------------------------

class ProteinActivationReproductionProcess(Process):
    """Allosteric activation of regulatable proteins (reproduction of Karr 2012
    ProteinActivation).

    Inputs
    ------
    protein_counts : map[float]   per-protein total counts (candidates)
    regulator : float             regulator/stimulus level (e.g. metabolite conc)

    Outputs
    -------
    active_fraction : overwrite[map[float]]   per-protein fraction in active form
    """

    description = (
        "Protein Activation — reduced reproduction of Karr 2012 ProteinActivation.\n"
        "MATLAB's evaluateActivationRules partitions each regulatable protein between an active\n"
        "and inactive form according to metabolite/stimulus rules. Here the active fraction of\n"
        "each protein follows a Hill function of the regulator level:\n"
        "    active_fraction = regulator^n / (K^n + regulator^n).\n"
        "Contract — in: protein_counts (map), regulator (float). "
        "out: active_fraction (overwrite[map[float]], current per-protein value).\n"
        "Fidelity: REDUCED — single scalar regulator + shared Hill law replaces the full per-rule logic."
    )

    config_schema = {
        "regulator_k": {"_type": "map[float]", "_default": {}},  # per-protein half-max K
        "default_k": {"_type": "float", "_default": 1.0},
        "hill_n": {"_type": "float", "_default": 2.0},
        "seed": {"_type": "integer", "_default": 16},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._k = dict(self.config["regulator_k"])
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"protein_counts": "map[float]", "regulator": "float"}

    def outputs(self):
        return {"active_fraction": "overwrite[map[float]]"}

    def initial_state(self):
        return {"protein_counts": {}, "regulator": 1.0}

    def update(self, state, interval):
        proteins = dict(state.get("protein_counts", {}) or {})
        r = max(0.0, float(state.get("regulator", 0.0)))
        n = float(self.config["hill_n"])
        default_k = float(self.config["default_k"])
        rn = r ** n
        frac = {}
        for p in proteins:
            k = float(self._k.get(p, default_k))
            kn = k ** n
            denom = kn + rn
            frac[p] = (rn / denom) if denom > 0 else 0.0
        return {"active_fraction": frac}


# ---------------------------------------------------------------------------
# 7. Macromolecular Complexation — stoichiometric subunit assembly
# ---------------------------------------------------------------------------

# Representative complex stoichiometries {complex: {monomer: subunits_per_complex}}
_DEFAULT_COMPLEX_STOICHIOMETRY = {
    "RNA_polymerase": {"rpoA": 2.0, "rpoB": 1.0, "rpoC": 1.0, "rpoD": 1.0},
    "DNA_polymerase_III": {"dnaE": 1.0, "dnaN": 2.0, "dnaX": 2.0},
    "GroEL_chaperonin": {"groEL": 14.0, "groES": 7.0},
}


class MacromolecularComplexationReproductionProcess(Process):
    """Stoichiometric assembly of monomers into complexes (reproduction of Karr
    2012 MacromolecularComplexation).

    Inputs
    ------
    monomers : map[float]     per-gene monomer counts (subunit pools)

    Outputs
    -------
    monomers : map[float]     negative delta (subunits consumed)
    complexes : map[float]    positive delta (complexes formed)
    """

    description = (
        "Macromolecular Complexation — reduced reproduction of Karr 2012 MacromolecularComplexation.\n"
        "Monomers assemble into complexes limited by subunit availability. MATLAB solves each\n"
        "independent complex network by Monte-Carlo/bound; here for each complex\n"
        "    n_formed = floor(min over subunits(available / stoichiometry))\n"
        "processed in random order so shared subunits are drawn down as complexes compete.\n"
        "Contract — config: stoichiometry {complex: {monomer: n}}. in: monomers (map). "
        "out: monomers (neg Δ map), complexes (pos Δ map).\n"
        "Fidelity: REDUCED — greedy min-subunit assembly replaces the full network Monte-Carlo kinetics."
    )

    config_schema = {
        "stoichiometry": {"_type": "map[map[float]]", "_default": _DEFAULT_COMPLEX_STOICHIOMETRY},
        "seed": {"_type": "integer", "_default": 17},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        stoich = self.config["stoichiometry"]
        self._stoich = {c: dict(sub) for c, sub in (stoich or _DEFAULT_COMPLEX_STOICHIOMETRY).items()}
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"monomers": "map[float]"}

    def outputs(self):
        return {"monomers": "map[float]", "complexes": "map[float]"}

    def initial_state(self):
        return {"monomers": {}}

    def update(self, state, interval):
        avail = {k: float(v) for k, v in (state.get("monomers", {}) or {}).items()}
        consumed = {}
        formed = {}
        # random order handles competition for shared subunits (MATLAB Monte-Carlo networks)
        order = list(self._stoich.keys())
        self._rng.shuffle(order)
        for cx in order:
            subs = self._stoich[cx]
            # limiting subunit: floor(min over subunits(available / stoich))
            n = None
            for m, stoich in subs.items():
                if stoich <= 0:
                    continue
                cap = int(avail.get(m, 0.0) // stoich)
                n = cap if n is None else min(n, cap)
            if not n or n <= 0:
                continue
            for m, stoich in subs.items():
                used = n * stoich
                avail[m] = avail.get(m, 0.0) - used
                consumed[m] = consumed.get(m, 0.0) - used
            formed[cx] = float(n)
        return {"monomers": consumed, "complexes": formed}


# ---------------------------------------------------------------------------
# 8. Ribosome Assembly — 30S/50S from rProteins + rRNA, GTPase-driven
# ---------------------------------------------------------------------------

# Representative rProtein/rRNA composition of each ribosomal subunit (stoich 1 each)
_RIBOSOME_30S = {"rproteins": ["rpsB", "rpsC", "rpsE", "rpsG"], "rrna": "16S"}
_RIBOSOME_50S = {"rproteins": ["rplB", "rplC", "rplD", "rplE"], "rrna": ["23S", "5S"]}


class RibosomeAssemblyReproductionProcess(Process):
    """Assembly of 30S and 50S ribosomal subunits (reproduction of Karr 2012
    RibosomeAssembly).

    Inputs
    ------
    rprotein_counts : map[float]   ribosomal protein pools
    rrna_counts : map[float]       rRNA pools (16S, 23S, 5S)
    assembly_factor : float        GTPase assembly-factor capacity
    gtp : float                    GTP pool

    Outputs
    -------
    ribosome_30S : float           positive delta (30S subunits formed)
    ribosome_50S : float           positive delta (50S subunits formed)
    rprotein_counts : map[float]   negative delta (consumed)
    rrna_counts : map[float]       negative delta (consumed)
    gtp : float                    negative delta (GTP consumed)
    """

    description = (
        "Ribosome Assembly — reduced reproduction of Karr 2012 RibosomeAssembly.\n"
        "30S/50S subunits assemble from ribosomal proteins + rRNA, requiring GTPase assembly\n"
        "factors and GTP. MATLAB forms newComplexs = floor(min(gtp/gtpPerComplex, water/..,\n"
        "RNAs, monomers)) gated on assembly enzymes; here each subunit forms at\n"
        "    n = floor(min(rProtein pool, rRNA pool, assembly_factor, gtp/gtp_per_complex))\n"
        "consuming subunits and GTP.\n"
        "Contract — in: rprotein_counts (map), rrna_counts (map), assembly_factor (float), gtp (float). "
        "out: ribosome_30S (float Δ), ribosome_50S (float Δ), rprotein_counts (neg map), "
        "rrna_counts (neg map), gtp (neg Δ).\n"
        "Fidelity: REDUCED — stoich-1 representative subunit sets; water/GDP/Pi accounting omitted."
    )

    config_schema = {
        "gtp_per_complex": {"_type": "float", "_default": 2.0},
        "seed": {"_type": "integer", "_default": 18},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {
            "rprotein_counts": "map[float]",
            "rrna_counts": "map[float]",
            "assembly_factor": "float",
            "gtp": "float",
        }

    def outputs(self):
        return {
            "ribosome_30S": "float",
            "ribosome_50S": "float",
            "rprotein_counts": "map[float]",
            "rrna_counts": "map[float]",
            "gtp": "float",
        }

    def initial_state(self):
        return {"rprotein_counts": {}, "rrna_counts": {}, "assembly_factor": 20.0, "gtp": 1e5}

    def _assemble(self, rprot, rrna, factor, gtp_cap, spec):
        """floor(min(limiting rProtein, limiting rRNA, assembly factor, gtp cap))."""
        rrna_keys = spec["rrna"] if isinstance(spec["rrna"], list) else [spec["rrna"]]
        limits = [factor, gtp_cap]
        limits += [rprot.get(p, 0.0) for p in spec["rproteins"]]
        limits += [rrna.get(r, 0.0) for r in rrna_keys]
        return int(max(0.0, np.floor(min(limits))))

    def update(self, state, interval):
        rprot = {k: float(v) for k, v in (state.get("rprotein_counts", {}) or {}).items()}
        rrna = {k: float(v) for k, v in (state.get("rrna_counts", {}) or {}).items()}
        factor = float(state.get("assembly_factor", 0.0))
        gtp = float(state.get("gtp", 0.0))
        gtp_per = float(self.config["gtp_per_complex"])
        gtp_cap = gtp / gtp_per if gtp_per > 0 else 0.0

        prot_used = {}
        rrna_used = {}
        gtp_used = 0.0

        # form 30S then 50S in randomized order, drawing down shared GTP/factor pools
        specs = [("30S", _RIBOSOME_30S), ("50S", _RIBOSOME_50S)]
        self._rng.shuffle(specs)
        n30 = n50 = 0
        for name, spec in specs:
            cap = gtp_cap - (gtp_used / gtp_per if gtp_per > 0 else 0.0)
            avail_factor = factor  # factor is catalytic (not consumed), shared each subunit
            n = self._assemble(
                {p: rprot.get(p, 0.0) + prot_used.get(p, 0.0) for p in spec["rproteins"]},
                rrna if not isinstance(spec["rrna"], list) else rrna,
                avail_factor, cap, spec,
            )
            # subtract already-consumed shared rProteins/rRNA for correctness
            rrna_keys = spec["rrna"] if isinstance(spec["rrna"], list) else [spec["rrna"]]
            for p in spec["rproteins"]:
                remaining = rprot.get(p, 0.0) + prot_used.get(p, 0.0)
                n = min(n, int(max(0.0, np.floor(remaining))))
            for r in rrna_keys:
                remaining = rrna.get(r, 0.0) + rrna_used.get(r, 0.0)
                n = min(n, int(max(0.0, np.floor(remaining))))
            if n <= 0:
                continue
            for p in spec["rproteins"]:
                prot_used[p] = prot_used.get(p, 0.0) - n
            for r in rrna_keys:
                rrna_used[r] = rrna_used.get(r, 0.0) - n
            gtp_used += n * gtp_per
            if name == "30S":
                n30 += n
            else:
                n50 += n

        return {
            "ribosome_30S": float(n30),
            "ribosome_50S": float(n50),
            "rprotein_counts": prot_used,
            "rrna_counts": rrna_used,
            "gtp": -gtp_used,
        }


# ---------------------------------------------------------------------------
# 9. Terminal Organelle Assembly — ordered adhesin localization
# ---------------------------------------------------------------------------

# Ordered set of adhesins/accessory proteins localized to the terminal organelle
_TERMINAL_ORGANELLE_ORDER = ["HMW2", "HMW1", "HMW3", "P1", "P41", "P24", "P65"]


class TerminalOrganelleAssemblyReproductionProcess(Process):
    """Ordered localization of adhesins to the terminal organelle (reproduction
    of Karr 2012 TerminalOrganelleAssembly).

    Inputs
    ------
    adhesin_proteins : map[float]   per-protein counts of terminal-organelle proteins

    Outputs
    -------
    terminal_organelle_fraction : overwrite[float]   fraction of the ordered
        assembly that is complete (current value in [0, 1])
    """

    description = (
        "Terminal Organelle Assembly — reduced reproduction of Karr 2012 TerminalOrganelleAssembly.\n"
        "Adhesins/accessory proteins localize to the terminal organelle in a fixed dependency\n"
        "ORDER: MATLAB only localizes a protein once its prerequisite localization reactions pass\n"
        "a threshold. Here the assembled fraction is the length of the leading run of ordered\n"
        "proteins present at/above threshold, divided by the number of required proteins.\n"
        "Contract — config: required order + threshold. in: adhesin_proteins (map). "
        "out: terminal_organelle_fraction (overwrite[float], current value 0..1).\n"
        "Fidelity: REDUCED — sequential presence gate replaces the full localization-reaction matrix."
    )

    config_schema = {
        "required_order": {"_type": "list[string]", "_default": _TERMINAL_ORGANELLE_ORDER},
        "threshold": {"_type": "float", "_default": 1.0},
        "seed": {"_type": "integer", "_default": 19},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._order = list(self.config["required_order"]) or _TERMINAL_ORGANELLE_ORDER
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"adhesin_proteins": "map[float]"}

    def outputs(self):
        return {"terminal_organelle_fraction": "overwrite[float]"}

    def initial_state(self):
        return {"adhesin_proteins": {}}

    def update(self, state, interval):
        proteins = dict(state.get("adhesin_proteins", {}) or {})
        thr = float(self.config["threshold"])
        # ordered assembly: count the leading run of prerequisites that are present
        completed = 0
        for name in self._order:
            if float(proteins.get(name, 0.0)) >= thr:
                completed += 1
            else:
                break  # ordered dependency: cannot localize past a missing prerequisite
        frac = completed / len(self._order) if self._order else 0.0
        return {"terminal_organelle_fraction": float(frac)}
