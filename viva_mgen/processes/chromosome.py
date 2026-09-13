"""Chromosome DNA–protein interaction submodel — clean-room reproduction.

Reproduces the coordinate-resolved chromosome dynamics behind Karr 2012 Figure 3:
DNA-binding proteins (RNA polymerase, DNA polymerase, DnaA, gyrase, SMC, SSB,
topoisomerase, transcription factors) occupy and move along the 580,070 bp
chromosome, and head-on / co-directional COLLISIONS between them are detected and
resolved. In the original this emerges from the ``Chromosome`` state
(``proteinBoundSites``, ``polymerizedRegions``) plus the RNA-polymerase and
replisome position tracking in the Transcription/Replication submodels.

This is a reduced but genuinely coordinate-resolved model (the earlier viva-Mgen
transcription/replication submodels are aspatial): the genome is binned; RNA
polymerases initiate at gene start sites and elongate; two replisomes advance
bidirectionally from oriC; structural proteins occupy sites; and when a moving
polymerase enters an occupied bin a collision is recorded by protein pair and the
weaker partner is displaced (RNA pol yields to DNA pol; SMC/SSB yield to pols) —
reproducing Fig 3's occupancy map (A), chromosome-exploration kinetics (B),
polymerase position traces (D), collision-frequency matrix (E), and the
collisions-vs-binding-density relationship (F).

Gene positions are assigned from M. genitalium locus-tag order (MG_### tags run
sequentially along the G37 chromosome), a genuine approximation of the true
coordinates (which live in the paper's genome annotation, not the pieces of the
knowledge base loadable here).
"""

from __future__ import annotations

import numpy as np
from process_bigraph import Process

from .. import constants as C
from ..kb import load_genes

# structural / DNA-binding proteins tracked for the Fig-3E collision matrix
_STRUCTURAL = ("SMC", "SSB", "GyrAB", "Topo IV", "DnaA", "Fur", "HrcA")


def _gene_bins(n_bins: int) -> np.ndarray:
    """Bin index of each modelled gene, from MG_### locus-tag order along G37."""
    nums = []
    for g in load_genes():
        gid = g["gene_id"]
        digits = "".join(ch for ch in gid if ch.isdigit())
        if digits:
            nums.append(int(digits))
    if not nums:
        return np.arange(n_bins)
    nums = np.array(sorted(nums), dtype=float)
    # locus number -> genome fraction -> bin
    frac = nums / nums.max()
    return np.clip((frac * (n_bins - 1)).astype(int), 0, n_bins - 1)


class ChromosomeDynamicsReproductionProcess(Process):
    """Coordinate-resolved chromosome occupancy, polymerase motion, and collisions.

    Inputs
    ------
    rna_polymerase : float
        Number of RNA polymerases available to be actively transcribing.
    replication_active : float
        1.0 while the replisomes elongate (drives DNA-pol chromosome exploration).

    Outputs (current-value snapshots)
    -------
    fraction_explored : overwrite[float]      fraction of bins ever protein-bound
    dna_binding_density : overwrite[float]     currently-bound bins / total (knt^-1-like)
    n_collisions : overwrite[float]            cumulative collision count
    percent_rnap_explored : overwrite[float]   fraction explored by RNA pol alone
    percent_dnap_explored : overwrite[float]   fraction explored by DNA pol alone
    occupancy : overwrite[map[float]]          cumulative bound-time per bin (bin->count)
    rna_pol_positions : overwrite[list[float]] current active RNA-pol bin positions
    dna_pol_positions : overwrite[list[float]] the two replisome bin positions
    collisions : overwrite[map[float]]         "Mover||Occupant" -> count
    """

    description = (
        "Chromosome DNA-protein interactions — reduced coordinate-resolved reproduction "
        "of the dynamics behind Karr 2012 Fig 3.\n"
        "The 580,070 bp chromosome is binned; RNA polymerases initiate at gene sites "
        "(locus-tag-ordered positions, rRNA hotspot weighted) and elongate at "
        "~50 nt/s; two replisomes advance bidirectionally from oriC at ~100 nt/s; "
        "structural proteins (SMC, SSB, gyrase, topoisomerase, DnaA, TFs) occupy sites; "
        "and when a moving polymerase enters an occupied bin a COLLISION is recorded by "
        "protein pair and the weaker partner displaced (RNA pol yields to the fork; "
        "SMC/SSB yield to polymerases).\n"
        "Contract — in: rna_polymerase (available RNA pols), replication_active (0/1). "
        "out (snapshots): occupancy (bin→bound-time), fraction_explored, "
        "percent_rnap/dnap_explored, rna_pol_positions, dna_pol_positions, collisions "
        "(pair→count), n_collisions, dna_binding_density.\n"
        "Fidelity: REDUCED — binned genome, approximate (locus-ordered) gene coordinates, "
        "representative protein counts; the aspatial transcription/replication submodels "
        "do not provide coordinates, so this process supplies the spatial layer."
    )

    config_schema = {
        "genome_length_bp": {"_type": "float", "_default": float(C.GENOME_LENGTH_BP)},
        "n_bins": {"_type": "integer", "_default": 580},
        "n_rna_pol": {"_type": "integer", "_default": 50},         # concurrently elongating pols
        "rna_pol_rate_bp": {"_type": "float", "_default": 50.0},   # nt/s
        "dna_pol_rate_bp": {"_type": "float", "_default": 20.0},   # nt/s per replisome → ~4 h to traverse the genome (paper's replication phase)
        "n_smc": {"_type": "integer", "_default": 280},   # condensin densely coats the chromosome
        "n_ssb": {"_type": "integer", "_default": 30},
        "n_gyrase": {"_type": "integer", "_default": 20},
        "n_topo": {"_type": "integer", "_default": 10},
        "n_tf": {"_type": "integer", "_default": 12},
        # fraction of each structural protein's sites that unbind + rebind to a
        # NEW random position per second (sets the chromosome-exploration rate).
        "struct_turnover_per_s": {"_type": "float", "_default": 0.008},
        # structural proteins bind PROGRESSIVELY (not all at t=0): the bound count
        # ramps as 1 − e^(−t/τ), so early exploration/density is low and the dense
        # SMC coating (and the RNA-pol↔SMC collisions it causes) builds over the
        # first ~20 min — matching the paper's gradual exploration + >30k collisions.
        "struct_bind_tau_s": {"_type": "float", "_default": 1200.0},
        # bins an RNA pol elongates across before release (longer span → the
        # polymerases sweep more of the genome per pass → 90%% coverage sooner).
        "rna_gene_span_bins": {"_type": "integer", "_default": 20},
        "seed": {"_type": "integer", "_default": 0},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rng = np.random.default_rng(int(self.config["seed"]))
        self.nb = int(self.config["n_bins"])
        self.bp_per_bin = self.config["genome_length_bp"] / self.nb
        self.gene_bins = _gene_bins(self.nb)
        # highly-transcribed rRNA-like hotspot: a small contiguous region gets
        # extra RNA-pol initiation weight (the paper's rRNA operon in Fig 3A)
        self._init_weight = np.ones(self.nb)
        rr = int(0.15 * self.nb)
        self._init_weight[rr:rr + 6] = 25.0
        self.oriC = 0
        self.terC = self.nb // 2
        # state
        self.occ = np.zeros(self.nb)               # cumulative bound-time per bin
        self.explored_any = np.zeros(self.nb, bool)
        self.explored_rnap = np.zeros(self.nb, bool)
        self.explored_dnap = np.zeros(self.nb, bool)
        # per-protein cumulative occupancy + exploration, read by the Fig-3 viz to
        # draw the circular per-protein binding-probability rings (A) and the
        # per-protein exploration curves (B). Not exposed as composite ports.
        self.occ_by = {}                           # label -> cumulative bound-time per bin
        self.explored_by = {}                      # label -> bool array (ever bound)
        # structural proteins bind PERSISTENTLY and turn over a small fraction of
        # their sites each second, so the chromosome is explored GRADUALLY (a
        # coupon-collector saturation) rather than being ~fully covered in the
        # first step — matching the paper's "50% bound by 6 min, 90% by 20 min".
        self._struct_sites = {}                    # label -> np.array of bound bins
        self._t_struct = 0.0                        # elapsed time (drives the binding ramp)
        self.rna_pols = []                          # list of [pos, end]
        self.dna_pos = None                         # [left, right] once replicating
        self.collisions = {}
        self.cur_bound = np.zeros(self.nb, bool)

    def inputs(self):
        return {"rna_polymerase": "float", "replication_active": "float"}

    def outputs(self):
        return {
            "fraction_explored": "overwrite[float]",
            "dna_binding_density": "overwrite[float]",
            "n_collisions": "overwrite[float]",
            "percent_rnap_explored": "overwrite[float]",
            "percent_dnap_explored": "overwrite[float]",
            "occupancy": "overwrite[map[float]]",
            "rna_pol_positions": "overwrite[list[float]]",
            "dna_pol_positions": "overwrite[list[float]]",
            "collisions": "overwrite[map[float]]",
        }

    def initial_state(self):
        return {"rna_polymerase": 120.0, "replication_active": 1.0}

    def _collide(self, mover, occupant):
        key = f"{mover}||{occupant}"
        self.collisions[key] = self.collisions.get(key, 0.0) + 1.0

    def update(self, state, interval):
        nb = self.nb
        self.cur_bound = np.zeros(nb, bool)
        occupant = {}  # bin -> protein label (last writer for collision checks)

        def place(bin_idx, label):
            bin_idx = int(bin_idx) % nb
            self.cur_bound[bin_idx] = True
            self.occ[bin_idx] += 1.0
            self.explored_any[bin_idx] = True
            occupant[bin_idx] = label
            if label not in self.occ_by:
                self.occ_by[label] = np.zeros(nb)
                self.explored_by[label] = np.zeros(nb, bool)
            self.occ_by[label][bin_idx] += 1.0
            self.explored_by[label][bin_idx] = True

        # 1. structural / DNA-binding proteins occupy sites PERSISTENTLY, turning
        #    over a fraction to new positions each step (gradual exploration).
        cfg = self.config
        place(self.oriC, "DnaA")
        self._t_struct += interval
        ramp = 1.0 - np.exp(-self._t_struct / max(float(cfg["struct_bind_tau_s"]), 1.0))
        turnover = float(cfg["struct_turnover_per_s"]) * interval
        for label, n_max in (("SMC", cfg["n_smc"]), ("SSB", cfg["n_ssb"]),
                             ("GyrAB", cfg["n_gyrase"]), ("Topo IV", cfg["n_topo"])):
            n_now = max(1, int(round(int(n_max) * ramp)))         # progressively-bound count
            sites = self._struct_sites.get(label)
            if sites is None:
                sites = self._rng.integers(0, nb, n_now)
            else:
                k = int(round(min(1.0, turnover) * len(sites)))  # rebind a fraction
                if k > 0:
                    idx = self._rng.choice(len(sites), k, replace=False)
                    sites = sites.copy()
                    sites[idx] = self._rng.integers(0, nb, k)
                if n_now > len(sites):                            # grow toward n_max
                    sites = np.concatenate([sites, self._rng.integers(0, nb, n_now - len(sites))])
            self._struct_sites[label] = sites
            for b in sites:
                place(b, label)
        # transcription factors at a few promoters
        for b in self._rng.choice(self.gene_bins, min(int(cfg["n_tf"]), len(self.gene_bins)), replace=False):
            place(b, "Fur")

        step_bins_rna = max(1, int(round(cfg["rna_pol_rate_bp"] * interval / self.bp_per_bin)))
        step_bins_dna = max(1, int(round(cfg["dna_pol_rate_bp"] * interval / self.bp_per_bin)))

        # 2. RNA polymerase: initiate up to target, elongate, release at gene end
        avail = int(min(state.get("rna_polymerase", cfg["n_rna_pol"]), cfg["n_rna_pol"]))
        while len(self.rna_pols) < avail:
            start = int(self._rng.choice(nb, p=self._init_weight / self._init_weight.sum()))
            length = int(self._rng.integers(1, max(2, int(cfg["rna_gene_span_bins"]))))  # gene spans a few bins
            self.rna_pols.append([start, min(nb - 1, start + length)])
        still = []
        for pol in self.rna_pols:
            pos, end = pol
            newpos = pos + step_bins_rna
            # collision check along the path
            for b in range(pos, min(newpos, end) + 1):
                bb = b % nb
                occ_label = occupant.get(bb)
                if occ_label and occ_label not in ("RNA Pol",):
                    if occ_label == "DNA Pol":
                        self._collide("DNA Pol", "RNA Pol")  # RNA pol displaced by fork
                        newpos = end + 1  # release this RNA pol
                        break
                    else:
                        self._collide("RNA Pol", occ_label)  # RNA pol displaces structural
                        self.cur_bound[bb] = False
            pol[0] = newpos
            if newpos <= end:
                place(newpos, "RNA Pol")
                self.explored_rnap[newpos % nb] = True
                still.append(pol)
            # else: reached end -> released (slot frees for a new initiation)
        self.rna_pols = still

        # 3. DNA polymerase: two replisomes from oriC, bidirectional, if replicating
        if state.get("replication_active", 0.0) >= 0.5:
            if self.dna_pos is None:
                self.dna_pos = [self.oriC, self.oriC]
            self.dna_pos[0] = self.dna_pos[0] + step_bins_dna         # rightward toward terC
            self.dna_pos[1] = self.dna_pos[1] - step_bins_dna         # leftward toward terC
            for i, dp in enumerate(self.dna_pos):
                bb = int(dp) % nb
                occ_label = occupant.get(bb)
                if occ_label and occ_label not in ("DNA Pol",):
                    self._collide("DNA Pol", occ_label)              # fork displaces anything
                place(bb, "DNA Pol")
                self.explored_dnap[bb] = True

        n_coll = float(sum(self.collisions.values()))
        return {
            "fraction_explored": float(self.explored_any.mean()),
            "dna_binding_density": float(self.cur_bound.mean()),
            "n_collisions": n_coll,
            "percent_rnap_explored": float(self.explored_rnap.mean()),
            "percent_dnap_explored": float(self.explored_dnap.mean()),
            "occupancy": {str(i): float(v) for i, v in enumerate(self.occ) if v > 0},
            "rna_pol_positions": [float(p[0] % nb) for p in self.rna_pols],
            "dna_pol_positions": [float(x % nb) for x in (self.dna_pos or [])],
            "collisions": dict(self.collisions),
        }
