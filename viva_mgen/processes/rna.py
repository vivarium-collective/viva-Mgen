"""RNA maturation & regulation submodels — clean-room reproduction of Karr 2012.

Clean-room viva-native reproductions of four RNA-handling submodels of the Karr
2012 *M. genitalium* whole-cell model, each written as a
:class:`process_bigraph.Process`. They mirror the essential mechanism of the
corresponding MATLAB submodel (``evolveState``). Where the knowledge base
carries kinetic constants they are used (RNAProcessing's per-RNase specific
rates); where it does not (RNAModification, tRNAAminoacylation,
TranscriptionalRegulation all have empty parameter dicts), the mechanism is
faithful but rate constants are order-of-magnitude values. Each class's
``description`` states which case it is. The per-species reaction/stoichiometry
matrices and metabolite byproduct accounting are delegated to the metabolite
pools.

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
from ..kb import karr_process_params


def _rnase_rates() -> dict:
    """Real Karr KB per-RNase specific rates (1/s) for RNAProcessing:
    RNase III, RNase P, RNase J, DeaD, RsgA."""
    p = karr_process_params("RNAProcessing")
    return {
        "RNAseIII": p.get("enzymeSpecificRate_RNAseIII", 7.7),
        "RNAseP": p.get("enzymeSpecificRate_RNAseP", 6.0),
        "RNAseJ": p.get("enzymeSpecificRate_RNAseJ", 0.37),
        "DeaD": p.get("enzymeSpecificRate_DeaD", 1.48),
        "RsgA": p.get("enzymeSpecificRate_RsgA", 0.2917),
    }


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
        "Transcriptional regulation — reproduction of Karr 2012 TranscriptionalRegulation.\n"
        "TFs bind accessible promoter sites and multiply each transcription unit's RNA-pol binding\n"
        "probability by a fold change. Per gene g:\n"
        "    fold_change_g = prod over regulating TF f of  (fc_{f,g} ** activity_f)\n"
        "mirroring calcBindingProbabilityFoldChange (foldChange = prod of tfActivities;\n"
        "otherFoldChanges = prod(otherActivities .^ (enzymes>0))). An activity of 0 leaves the gene\n"
        "unregulated (fold change 1); an activated TF pushes fc_{f,g}>1, a repressor fc_{f,g}<1.\n"
        "The default network is the REAL KB regulatory network: 52 genes regulated by the 5 genuine\n"
        "M. genitalium TFs (MG_127/MG_236/MG_101 monomers, MG_205/MG_428 dimers) with the true\n"
        "per-edge fold-changes decoded from the TranscriptionUnit objects.\n"
        "Contract — config: gene_regulation (default = real KB network). in: tf_activity (per-TF "
        "fractional activity map). out: fold_change (per-gene multiplicative fold change, an "
        "overwrite sensor map).\n"
        "Fidelity: FAITHFUL regulatory network (real KB TF→gene edges + real fold-changes). The\n"
        "reduction is the promoter-occupancy / chromosome-accessibility state machine, collapsed to\n"
        "the fold-change law above driven by continuous TF activity."
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
        """The REAL Karr KB transcription-factor network: 52 genes regulated by 5
        M. genitalium TFs (MG_127, MG_236, MG_101 monomers; MG_205, MG_428 dimers)
        with the genuine per-edge fold-changes. Keyed by the expression-panel key
        (gene symbol when known, else gene id) so it composes with transcription."""
        from ..kb import load_genes, load_tf_regulation
        reg_by_id = load_tf_regulation()
        if not reg_by_id:  # fall back only if the dataset is unavailable
            return {g: {"rpoD_sigma": 1.0} for g in DEFAULT_GENES}
        # map gene_id -> panel key (symbol or id), matching expression_defaults
        panel_key = {}
        for g in load_genes():
            panel_key[g["gene_id"]] = (g.get("symbol") or "").strip() or g["gene_id"]
        reg = {}
        for gid, edges in reg_by_id.items():
            key = panel_key.get(gid, gid)
            reg[key] = dict(edges)
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
        "RNA processing — reproduction of Karr 2012 RNAProcessing.\n"
        "Ribonucleases (RNase III/P/J, RsgA, DeaD) cleave nascent polycistronic r/t/tmRNA\n"
        "transcripts into mature species. Per nascent species s:\n"
        "    matured_s ~ min(n_s, Poisson(n_s · k_proc · Δt))\n"
        "counts moving nascent_rna → mature_rna. Mirrors evolveState_Helper, which stochastically\n"
        "matures unprocessedRNAs → processedRNAs up to enzyme/substrate limits. The pooled rate\n"
        "k_proc is the mean of the REAL KB per-RNase specific rates (RNase III 7.7, P 6.0, J 0.37,\n"
        "DeaD 1.48, RsgA 0.29 s⁻¹ → ~3.17 s⁻¹).\n"
        "Contract — config: processing_rate (default = KB per-RNase mean), rnase_rates (the real\n"
        "per-RNase map). in: nascent_rna (per-species count map). "
        "out: nascent_rna (negative Δ), mature_rna (positive Δ).\n"
        "Fidelity: FAITHFUL rate constants (real KB per-RNase kcats). The lumping is species→RNase\n"
        "assignment (one pooled first-order rate) — the panel does not carry which RNase cleaves\n"
        "which transcript; cofactor (Mg²⁺/Zn²⁺) / ATP-GTP and intergenic-fragment accounting is\n"
        "delegated to the metabolite pools."
    )

    config_schema = {
        # pooled RNase kcat = mean of the real KB per-RNase specific rates
        "processing_rate": {"_type": "float", "_default": sum(_rnase_rates().values()) / 5.0},
        "rnase_rates": {"_type": "map[float]", "_default": _rnase_rates()},
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
        "RNA modification — reproduction of Karr 2012 RNAModification.\n"
        "Modification enzymes formylate / lysidinate / methylate / pseudouridylate / thiolate\n"
        "specific r/tRNA bases. Per species s:\n"
        "    modified_s ~ min(n_s, Poisson(n_s · k_mod · Δt)),   sum_s modified_s ≤ enzyme · kcat · Δt\n"
        "moving counts unmodified_rna → modified_rna, with the total capped by enzyme availability.\n"
        "Mirrors evolveState's greedy loop that modifies unmodifiedRNAs → modifiedRNAs while enzyme\n"
        "and substrate limits allow (reactionLimits gated on min over enzymes/substrates).\n"
        "Contract — in: unmodified_rna (count map), modification_enzyme (available enzyme count). "
        "out: unmodified_rna (negative Δ), modified_rna (positive Δ).\n"
        "Fidelity: mechanism-faithful (enzyme-budget-limited stochastic modification). The Karr KB\n"
        "RNAModification parameter dict is empty, so k_mod / kcat are order-of-magnitude values (not\n"
        "fabricated KB constants); one enzyme pool stands in for the 13 modification enzymes over 86\n"
        "base modifications, and metabolite substrate/byproduct accounting is delegated to the pools."
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
        "tRNA aminoacylation — reproduction of Karr 2012 tRNAAminoacylation.\n"
        "Aminoacyl-tRNA synthetases conjugate amino acids to free tRNAs at a cost of 1 ATP each\n"
        "(ATP -> AMP + PPi). The total number charged this step is\n"
        "    N = floor( min( sum(free_tRNA), amino_acid, atp, synthetase · kcat · Δt ) )\n"
        "distributed across species by multinomial draw proportional to free-tRNA counts, then\n"
        "free_trna -> aminoacylated_trna and atp decremented by N. Mirrors evolveState's greedy\n"
        "reactionLimits loop (min over synthetase/AA/ATP availability, weighted stochastic pick),\n"
        "as one aggregate limit + multinomial partition.\n"
        "Contract — in: free_trna (count map), amino_acid, atp, synthetase (pools/counts). "
        "out: free_trna (negative Δ), aminoacylated_trna (positive Δ), atp (negative Δ).\n"
        "ORDERING: runs BEFORE Translation (charged tRNAs feed the ribosome same-tick).\n"
        "Fidelity: mechanism-faithful (co-substrate-limited charging, ATP-coupled). The Karr KB\n"
        "tRNAAminoacylation parameter dict is empty, so the synthetase kcat is an order-of-magnitude\n"
        "value; one lumped synthetase pool + aggregate AA/ATP pools stand in for the 20 AA × 37\n"
        "tRNA/tmRNA reactions with per-synthetase kcats and glutamyl/methionyl transferases."
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
