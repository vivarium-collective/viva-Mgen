"""RNA maturation & regulation submodels — clean-room reproduction (reduced).

Clean-room viva-native reproductions of four RNA-handling submodels of the Karr
2012 *M. genitalium* whole-cell model, each written as a
:class:`process_bigraph.Process`. They mirror the essential mechanism of the
corresponding MATLAB submodel (``evolveState``) at reduced fidelity — the
qualitative flux is genuine, but the exhaustive per-species reaction/stoichiometry
matrices and knowledge-base kinetic constants of the originals are collapsed to
representative rates.

Submodels reproduced (original → class):
- ``TranscriptionalRegulation.m`` → :class:`TranscriptionalRegulationReproductionProcess`
- ``RNAProcessing.m``             → :class:`RNAProcessingReproductionProcess`
- ``RNAModification.m``           → :class:`RNAModificationReproductionProcess`
- ``tRNAAminoacylation.m``        → :class:`TRNAAminoacylationReproductionProcess`

Convention (mirroring transcription.py / decay.py): BARE composable port types —
``map[float]`` / ``float`` for additive deltas, ``overwrite[...]`` only for a
current-value sensor readout. numpy ``default_rng`` seeded from config ``seed``
where the mechanism is stochastic.
"""

from __future__ import annotations

import numpy as np
from process_bigraph import Process

from ..expression_defaults import DEFAULT_GENES


# ---------------------------------------------------------------------------
# 1. Transcriptional regulation
# ---------------------------------------------------------------------------
class TranscriptionalRegulationReproductionProcess(Process):
    """TF binding → per-gene RNA-pol binding-probability fold change.

    Inputs
    ------
    tf_activity : map[float]
        Fractional activity (~occupancy, 0..1) of each transcription factor.

    Outputs
    -------
    fold_change : overwrite[map[float]]
        Per-gene multiplicative fold change on RNA-polymerase promoter binding
        probability, ready to be consumed by TranscriptionReproductionProcess.
    """

    description = (
        "Transcriptional regulation — reduced reproduction of Karr 2012 TranscriptionalRegulation.\n"
        "TFs bind accessible promoter sites and multiply each transcription unit's RNA-pol binding\n"
        "probability by a fold change. Reduced mechanism, per gene g:\n"
        "    fold_change_g = prod over regulating TF f of  (fc_{f,g} ** activity_f)\n"
        "mirroring calcBindingProbabilityFoldChange (foldChange = prod of tfActivities;\n"
        "otherFoldChanges = prod(otherActivities .^ (enzymes>0))). An activity of 0 leaves the gene\n"
        "unregulated (fold change 1); an activated TF pushes fc_{f,g}>1, a repressor fc_{f,g}<1.\n"
        "Contract — in: tf_activity (per-TF fractional activity map). "
        "out: fold_change (per-gene multiplicative Δ-free fold change, an overwrite sensor map).\n"
        "Fidelity: REDUCED — a handful of representative TF→gene edges with representative fold\n"
        "changes, not the full promoter-occupancy / chromosome-accessibility state machine."
    )

    config_schema = {
        # gene -> {tf_name: fold_change_at_full_activity}. >1 activates, <1 represses.
        "gene_regulation": {"_type": "map[map[float]]", "_default": {}},
        "seed": {"_type": "integer", "_default": 0},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._reg = {g: dict(m) for g, m in (self.config["gene_regulation"] or {}).items()}
        if not self._reg:
            self._reg = self._default_regulation()
        self._rng = np.random.default_rng(int(self.config["seed"]))

    @staticmethod
    def _default_regulation() -> dict:
        # Representative TF→gene edges on the expression panel: one activator
        # (rpoD, the sigma factor, boosts broadly) and one repressor (hrcA-like
        # heat-shock repressor damping the chaperonin groEL).
        genes = DEFAULT_GENES
        reg = {g: {"rpoD_sigma": 1.4} for g in genes}
        if "groEL" in reg:
            reg["groEL"]["hrcA_repressor"] = 0.3
        if "ftsZ" in reg:
            reg["ftsZ"]["ftsZ_activator"] = 1.8
        return reg

    def inputs(self):
        return {"tf_activity": "map[float]"}

    def outputs(self):
        return {"fold_change": "overwrite[map[float]]"}

    def initial_state(self):
        return {"fold_change": {g: 1.0 for g in self._reg}}

    def update(self, state, interval):
        activity = dict(state.get("tf_activity", {}) or {})
        fold = {}
        for gene, tf_map in self._reg.items():
            # foldChange_g = prod_f fc_{f,g} ** activity_f   (MATLAB: product of
            # bound tfActivities; here activity is continuous occupancy in [0,1]).
            fc = 1.0
            for tf, base_fc in tf_map.items():
                a = max(0.0, float(activity.get(tf, 0.0)))
                if a > 0.0 and base_fc > 0.0:
                    fc *= base_fc ** a
            fold[gene] = fc
        return {"fold_change": fold}


# ---------------------------------------------------------------------------
# 2. RNA processing (cleavage of polycistronic transcripts)
# ---------------------------------------------------------------------------
class RNAProcessingReproductionProcess(Process):
    """RNase cleavage of nascent r/t/tmRNA transcripts into mature species.

    Inputs
    ------
    nascent_rna : map[float]
        Counts of nascent (unprocessed) transcripts per species.

    Outputs
    -------
    nascent_rna : map[float]
        Negative-Δ map — nascent transcripts consumed by cleavage.
    mature_rna : map[float]
        Positive-Δ map — mature species produced (same keys as nascent).
    """

    description = (
        "RNA processing — reduced reproduction of Karr 2012 RNAProcessing.\n"
        "Ribonucleases (RNase III/P/J, RsgA, DeaD) cleave nascent polycistronic r/t/tmRNA\n"
        "transcripts into mature species. Reduced mechanism, per nascent species s:\n"
        "    matured_s ~ min(n_s, Poisson(n_s · k_proc · Δt))\n"
        "counts moving nascent_rna → mature_rna. Mirrors evolveState_Helper, which stochastically\n"
        "matures unprocessedRNAs → processedRNAs up to enzyme/substrate limits; here a single\n"
        "first-order processing rate stands in for the pooled RNase kcats.\n"
        "Contract — in: nascent_rna (per-species count map). "
        "out: nascent_rna (negative Δ), mature_rna (positive Δ).\n"
        "Fidelity: REDUCED — one lumped first-order cleavage rate, no per-enzyme kcat / cofactor\n"
        "(Mg2+, Zn2+) / ATP-GTP accounting and no intergenic-fragment bookkeeping."
    )

    config_schema = {
        "processing_rate": {"_type": "float", "_default": 0.5},  # per second, pooled RNase kcat
        "seed": {"_type": "integer", "_default": 1},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rate = float(self.config["processing_rate"])
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"nascent_rna": "map[float]"}

    def outputs(self):
        return {"nascent_rna": "map[float]", "mature_rna": "map[float]"}

    def initial_state(self):
        return {"nascent_rna": {}, "mature_rna": {}}

    def update(self, state, interval):
        nascent = dict(state.get("nascent_rna", {}) or {})
        consumed = {}
        produced = {}
        for sp, n in nascent.items():
            n = float(n)
            if n <= 0:
                continue
            # Poisson first-order maturation, capped by available nascent count —
            # the reduced form of the greedy enzyme-limited cleavage loop.
            expected = n * self._rate * interval
            matured = min(n, float(self._rng.poisson(max(expected, 0.0))))
            if matured > 0:
                consumed[sp] = -matured
                produced[sp] = matured
        return {"nascent_rna": consumed, "mature_rna": produced}


# ---------------------------------------------------------------------------
# 3. RNA modification (methylation / pseudouridylation / etc.)
# ---------------------------------------------------------------------------
class RNAModificationReproductionProcess(Process):
    """Enzymatic modification of unmodified rRNA/tRNA → modified species.

    Inputs
    ------
    unmodified_rna : map[float]
        Counts of unmodified rRNA/tRNA species awaiting modification.
    modification_enzyme : float
        Available modification-enzyme count (methyltransferases, pseudouridine
        synthases, ...). Caps total modifications this step.

    Outputs
    -------
    unmodified_rna : map[float]
        Negative-Δ map — species consumed by modification.
    modified_rna : map[float]
        Positive-Δ map — modified species produced (same keys).
    """

    description = (
        "RNA modification — reduced reproduction of Karr 2012 RNAModification.\n"
        "Modification enzymes formylate / lysidinate / methylate / pseudouridylate / thiolate\n"
        "specific r/tRNA bases. Reduced mechanism, per species s:\n"
        "    modified_s ~ min(n_s, Poisson(n_s · k_mod · Δt)),   sum_s modified_s ≤ enzyme · kcat · Δt\n"
        "moving counts unmodified_rna → modified_rna, with the total capped by enzyme availability.\n"
        "Mirrors evolveState's greedy loop that modifies unmodifiedRNAs → modifiedRNAs while enzyme\n"
        "and substrate limits allow (reactionLimits gated on min over enzymes/substrates).\n"
        "Contract — in: unmodified_rna (count map), modification_enzyme (available enzyme count). "
        "out: unmodified_rna (negative Δ), modified_rna (positive Δ).\n"
        "Fidelity: REDUCED — a single enzyme pool + one per-species rate stand in for 13 enzymes\n"
        "over 86 specific base modifications; metabolite substrate/byproduct accounting omitted."
    )

    config_schema = {
        "modification_rate": {"_type": "float", "_default": 0.4},  # per second per RNA
        "enzyme_kcat": {"_type": "float", "_default": 2.0},        # modifications / enzyme / second
        "seed": {"_type": "integer", "_default": 2},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rate = float(self.config["modification_rate"])
        self._kcat = float(self.config["enzyme_kcat"])
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {"unmodified_rna": "map[float]", "modification_enzyme": "float"}

    def outputs(self):
        return {"unmodified_rna": "map[float]", "modified_rna": "map[float]"}

    def initial_state(self):
        return {"unmodified_rna": {}, "modification_enzyme": 50.0}

    def update(self, state, interval):
        unmod = dict(state.get("unmodified_rna", {}) or {})
        enzyme = float(state.get("modification_enzyme", 0.0))
        # Global enzyme-limited budget for this step (reactionLimits gated on enzyme count).
        budget = max(0.0, enzyme * self._kcat * interval)
        consumed = {}
        produced = {}
        for sp, n in unmod.items():
            n = float(n)
            if n <= 0 or budget <= 0:
                continue
            # Poisson per-species propensity, then clip to remaining enzyme budget.
            expected = n * self._rate * interval
            modified = min(n, float(self._rng.poisson(max(expected, 0.0))))
            modified = min(modified, budget)
            if modified > 0:
                consumed[sp] = -modified
                produced[sp] = modified
                budget -= modified
        return {"unmodified_rna": consumed, "modified_rna": produced}


# ---------------------------------------------------------------------------
# 4. tRNA aminoacylation (charging) — runs BEFORE Translation
# ---------------------------------------------------------------------------
class TRNAAminoacylationReproductionProcess(Process):
    """Synthetases charge free tRNAs with amino acids (consumes ATP).

    Inputs
    ------
    free_trna : map[float]
        Counts of free (uncharged) tRNA species.
    amino_acid : float
        Available amino-acid pool (charging consumes one AA per tRNA).
    atp : float
        Available ATP pool (charging costs 1 ATP per aminoacylation).
    synthetase : float
        Available aminoacyl-tRNA synthetase count (caps the charging rate).

    Outputs
    -------
    free_trna : map[float]
        Negative-Δ map — free tRNAs consumed by charging.
    aminoacylated_trna : map[float]
        Positive-Δ map — charged tRNAs produced (same keys).
    atp : float
        Negative delta — ATP consumed (1 per aminoacylation).

    Ordering note
    -------------
    In the original whole-cell model tRNAAminoacylation runs *before* Translation
    each time step, so charged tRNAs are available to the ribosome within the same
    tick. When composing with a translation Process, schedule this process first.
    """

    description = (
        "tRNA aminoacylation — reduced reproduction of Karr 2012 tRNAAminoacylation.\n"
        "Aminoacyl-tRNA synthetases conjugate amino acids to free tRNAs at a cost of 1 ATP each\n"
        "(ATP -> AMP + PPi). Reduced mechanism: the total number charged this step is\n"
        "    N = floor( min( sum(free_tRNA), amino_acid, atp, synthetase · kcat · Δt ) )\n"
        "distributed across species by multinomial draw proportional to free-tRNA counts, then\n"
        "free_trna -> aminoacylated_trna and atp decremented by N. Mirrors evolveState's greedy\n"
        "reactionLimits loop (min over synthetase/AA/ATP availability, weighted stochastic pick),\n"
        "collapsed to one aggregate limit + multinomial partition.\n"
        "Contract — in: free_trna (count map), amino_acid, atp, synthetase (pools/counts). "
        "out: free_trna (negative Δ), aminoacylated_trna (positive Δ), atp (negative Δ).\n"
        "ORDERING: runs BEFORE Translation (charged tRNAs feed the ribosome same-tick).\n"
        "Fidelity: REDUCED — one lumped synthetase pool + generic amino-acid/ATP pools, not 20 AAs\n"
        "× 37 tRNA/tmRNA reactions with per-synthetase kcats and glutamyl/methionyl transferases."
    )

    config_schema = {
        "synthetase_kcat": {"_type": "float", "_default": 20.0},  # charges / synthetase / second
        "atp_per_charge": {"_type": "float", "_default": 1.0},    # 1 ATP per aminoacylation
        "seed": {"_type": "integer", "_default": 3},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._kcat = float(self.config["synthetase_kcat"])
        self._atp_cost = float(self.config["atp_per_charge"])
        self._rng = np.random.default_rng(int(self.config["seed"]))

    def inputs(self):
        return {
            "free_trna": "map[float]",
            "amino_acid": "float",
            "atp": "float",
            "synthetase": "float",
        }

    def outputs(self):
        return {
            "free_trna": "map[float]",
            "aminoacylated_trna": "map[float]",
            "atp": "float",
        }

    def initial_state(self):
        return {"free_trna": {}, "amino_acid": 1e5, "atp": 1e6, "synthetase": 100.0}

    def update(self, state, interval):
        free = {k: float(v) for k, v in (state.get("free_trna", {}) or {}).items() if float(v) > 0}
        total_free = sum(free.values())
        if total_free <= 0:
            return {"free_trna": {}, "aminoacylated_trna": {}, "atp": 0.0}

        amino_acid = float(state.get("amino_acid", 0.0))
        atp = float(state.get("atp", 0.0))
        synthetase = float(state.get("synthetase", 0.0))

        # Aggregate reaction limit = min over the co-substrate availabilities
        # (free tRNA, amino acid, ATP, synthetase·kcat·Δt) — the reduced form of
        # the MATLAB greedy loop's reactionLimits (min across species columns).
        enzyme_limit = synthetase * self._kcat * interval
        atp_limit = atp / self._atp_cost if self._atp_cost > 0 else atp
        n_total = int(np.floor(min(total_free, amino_acid, atp_limit, enzyme_limit)))
        if n_total <= 0:
            return {"free_trna": {}, "aminoacylated_trna": {}, "atp": 0.0}

        # Partition the charged total across species proportional to free-tRNA
        # counts (multinomial) — the stochastic weighted pick over species.
        species = list(free.keys())
        weights = np.array([free[s] for s in species], dtype=float)
        weights /= weights.sum()
        draws = self._rng.multinomial(n_total, weights)

        consumed = {}
        produced = {}
        charged = 0
        for sp, d in zip(species, draws):
            d = min(float(d), free[sp])  # never charge more than are free
            if d > 0:
                consumed[sp] = -d
                produced[sp] = d
                charged += d

        return {
            "free_trna": consumed,
            "aminoacylated_trna": produced,
            "atp": -charged * self._atp_cost,
        }
