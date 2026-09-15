"""Chromosome DNA–protein interaction submodel — clean-room reproduction.

Reproduces the coordinate-resolved chromosome dynamics behind Karr 2012 Figure 3:
DNA-binding proteins (RNA polymerase, DNA polymerase, DnaA, gyrase, SMC, SSB,
topoisomerase, transcription factors) occupy and move along the 580,070 bp
chromosome, and head-on / co-directional COLLISIONS between them are detected and
resolved. In the original this emerges from the ``Chromosome`` state
(``proteinBoundSites``, ``polymerizedRegions``) plus the RNA-polymerase and
replisome position tracking in the Transcription/Replication submodels.

This is a genuinely coordinate-resolved model (the earlier viva-Mgen
transcription/replication submodels are aspatial): the genome is binned; RNA
polymerases initiate at gene start sites and elongate; two replisomes advance
bidirectionally from oriC; structural proteins occupy sites; and when a moving
polymerase enters an occupied bin a collision is recorded by protein pair and the
weaker partner is displaced (RNA pol yields to DNA pol; SMC/SSB yield to pols) —
reproducing Fig 3's occupancy map (A), chromosome-exploration kinetics (B),
polymerase position traces (D), collision-frequency matrix (E), and the
collisions-vs-binding-density relationship (F).

Gene positions are the REAL per-gene start coordinates on the 580,070 bp G37
chromosome, decoded natively from the knowledge base into
``datasets/karr_gene_expression.csv`` (kb.load_gene_expression); the rRNA
initiation hotspot of Fig 3A is placed at the real 16S/23S/5S operon locus
(~170–175 kb). Binning at n_bins resolution is the only spatial approximation.
"""

from __future__ import annotations

import numpy as np
from process_bigraph import Process

from .. import constants as C
from ..kb import load_gene_expression, load_genes

# structural / DNA-binding proteins tracked for the Fig-3E collision matrix
_STRUCTURAL = ("SMC", "SSB", "GyrAB", "Topo IV", "DnaA", "Fur", "HrcA")


def _gene_bins(n_bins: int, genome_length_bp: float) -> np.ndarray:
    """Bin index of each modelled gene from its REAL start coordinate on the
    580,070 bp G37 chromosome (decoded from the KB into karr_gene_expression.csv).
    Falls back to locus-tag order only if coordinates are unavailable."""
    expr = load_gene_expression()
    starts = [d["start_coordinate"] for d in expr.values()
              if d.get("start_coordinate") is not None]
    if starts:
        frac = np.array(sorted(starts), dtype=float) / float(genome_length_bp)
        return np.clip((frac * (n_bins - 1)).astype(int), 0, n_bins - 1)
    # fallback: MG_### locus-tag order along G37
    nums = []
    for g in load_genes():
        digits = "".join(ch for ch in g["gene_id"] if ch.isdigit())
        if digits:
            nums.append(int(digits))
    if not nums:
        return np.arange(n_bins)
    nums = np.array(sorted(nums), dtype=float)
    return np.clip((nums / nums.max() * (n_bins - 1)).astype(int), 0, n_bins - 1)


def _rrna_bins(n_bins: int, genome_length_bp: float) -> list:
    """Bin indices covering the real rRNA operon (16S/23S/5S) from the KB
    coordinates — the highly-transcribed RNA-pol initiation hotspot of Fig 3A."""
    expr = load_gene_expression()
    bins = set()
    for d in expr.values():
        if d.get("rna_type") == "rRNA" and d.get("start_coordinate") is not None:
            s = int(d["start_coordinate"] / genome_length_bp * (n_bins - 1))
            e = int((d.get("end_coordinate") or d["start_coordinate"])
                    / genome_length_bp * (n_bins - 1))
            bins.update(range(max(0, s), min(n_bins, e + 1)))
    return sorted(bins)


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
    rna_pol_positions : overwrite[list[float]] active RNA-pol positions, SIGNED bins from oriC
    dna_pol_positions : overwrite[list[float]] the two replisome positions, SIGNED bins from oriC
                                               ([+off, -off]: forks diverge from oriC to ±terC)
    collisions : overwrite[map[float]]         "Binder||Resident" -> count (Fig 4E panel)
    """

    description = (
        "Chromosome DNA-protein interactions — coordinate-resolved reproduction "
        "of the dynamics behind Karr 2012 Fig 3.\n"
        "The 580,070 bp chromosome is binned; RNA polymerases initiate at REAL per-gene "
        "start coordinates (KB-decoded, rRNA operon hotspot at its real ~170 kb locus) and "
        "elongate at ~50 nt/s; two replisomes initiate at oriC and advance OUTWARD, "
        "symmetrically, to terC (reported as signed positions ±off from oriC — a V, as in "
        "Fig 4D); the full DNA-binding-protein panel of Fig 4E occupies sites (SMC, SSB, "
        "GyrAB, Topo IV, DnaB, DnaN, DnaA, and the Fur/GntR/HrcA/LuxR transcription factors); "
        "and whenever a protein binds an already-occupied bin a binding×unbinding "
        "COLLISION is recorded by protein pair (RNA pol yields to the fork).\n"
        "Contract — in: rna_polymerase (available RNA pols), replication_active (0/1), "
        "lesion_map (the SHARED per-site chromosome — damaged bins from DNADamage). "
        "out (snapshots): occupancy (bin→bound-time), fraction_explored, "
        "percent_rnap/dnap_explored, rna_pol_positions, dna_pol_positions, collisions "
        "(pair→count, incl. *||lesion stalls), n_collisions, dna_binding_density.\n"
        "Fidelity: FAITHFUL gene coordinates (real KB per-gene loci + real rRNA operon "
        "position) and bidirectional oriC→terC fork geometry; the spatial approximation is the "
        "genome binning (n_bins), and structural-protein counts are representative. This process "
        "supplies the spatial layer the aspatial transcription/replication submodels lack, and "
        "reads the SHARED chromosome structure (gap #3): a DNA lesion on a polymerase's path "
        "stalls it (damage blocks elongation / the fork), so damage and the coordinate-resolved "
        "dynamics operate on ONE chromosome. Baseline damaging_agent=0 → no lesions → unchanged."
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
        "n_dnab": {"_type": "integer", "_default": 4},     # replicative helicase (DnaB)
        "n_dnan": {"_type": "integer", "_default": 4},     # sliding clamp (DnaN / beta-clamp)
        "n_tf": {"_type": "integer", "_default": 12},      # transcription factors, split across Fur/GntR/HrcA/LuxR
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
        # fraction of structural REBINDS that are treated as displacement
        # collisions. The high turnover above is an exploration device; most
        # rebinding lands on free DNA, so only a minority displaces a bound
        # protein — this keeps the moving polymerases the dominant collision
        # cause (Karr Fig 3F: ~84%% by RNA pol) while still populating the Fig-4E
        # structural rows/columns.
        "struct_collision_frac": {"_type": "float", "_default": 0.12},
        "seed": {"_type": "integer", "_default": 0},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._rng = np.random.default_rng(int(self.config["seed"]))
        self.nb = int(self.config["n_bins"])
        self.bp_per_bin = self.config["genome_length_bp"] / self.nb
        self.gene_bins = _gene_bins(self.nb, self.config["genome_length_bp"])
        # highly-transcribed rRNA hotspot: the REAL rRNA operon bins (16S/23S/5S,
        # KB coordinates ~170–175 kb) get extra RNA-pol initiation weight (Fig 3A)
        self._init_weight = np.ones(self.nb)
        rrna = _rrna_bins(self.nb, self.config["genome_length_bp"])
        if not rrna:  # fallback if coordinates unavailable
            rr = int(0.15 * self.nb)
            rrna = list(range(rr, rr + 6))
        for b in rrna:
            self._init_weight[b] = 25.0
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
        self.dna_off = None                         # fork offset from oriC (bins), once replicating
        self.collisions = {}
        self.cur_bound = np.zeros(self.nb, bool)
        # the full DNA-binding-protein panel of Karr 2012 Fig 4E (binding ×
        # unbinding collision matrix): the persistently-bound structural set +
        # the two replisome accessory factors + the transcription factors.
        self._tf_labels = ("Fur", "GntR", "HrcA", "LuxR")

    def inputs(self):
        # lesion_map is the SHARED per-site chromosome structure (gap #3): damaged
        # sites (from DNADamage) that this process reads so damage and the
        # coordinate-resolved dynamics operate on ONE chromosome. Baseline
        # damaging_agent=0 → empty map → no obstacles → behaviour unchanged.
        return {"rna_polymerase": "float", "replication_active": "float",
                "lesion_map": "map[float]"}

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
        return {"rna_polymerase": 120.0, "replication_active": 1.0, "lesion_map": {}}

    def _collide(self, mover, occupant):
        key = f"{mover}||{occupant}"
        self.collisions[key] = self.collisions.get(key, 0.0) + 1.0

    def _signed(self, bin_idx):
        """Position as a SIGNED offset from oriC (bins): oriC = 0 in the centre,
        terC = ±nb/2 at both ends. This is the Karr 2012 Fig 4D convention — the
        two replication forks diverge from oriC (centre) outward to terC (a V),
        rather than wrapping the [0, nb) axis into a crossing X."""
        b = int(bin_idx) % self.nb
        return b if b <= self.terC else b - self.nb

    def update(self, state, interval):
        nb = self.nb
        self.cur_bound = np.zeros(nb, bool)
        occupant = {}  # bin -> RESIDENT protein label (first binder this step)
        # SHARED chromosome (gap #3): bins carrying a DNA lesion act as physical
        # obstacles that stall an advancing polymerase (damage blocks elongation /
        # the replication fork). Read from the same per-site structure DNADamage/
        # DNARepair use. Empty at baseline (damaging_agent=0) → this set is empty
        # and every `in lesioned` test below is False → behaviour is unchanged.
        lesioned = {int(b) % nb for b, v in (state.get("lesion_map", {}) or {}).items()
                    if float(v) > 0.0}

        def place(bin_idx, label, record=True):
            bin_idx = int(bin_idx) % nb
            resident = occupant.get(bin_idx)
            if record and resident is not None and resident != label:
                # `label` binds a site already held by `resident` → a
                # binding × unbinding collision (Fig 4E). The resident stays the
                # site's occupant; the collision is what the matrix counts.
                self._collide(label, resident)
            self.cur_bound[bin_idx] = True
            self.occ[bin_idx] += 1.0
            self.explored_any[bin_idx] = True
            occupant.setdefault(bin_idx, label)
            if label not in self.occ_by:
                self.occ_by[label] = np.zeros(nb)
                self.explored_by[label] = np.zeros(nb, bool)
            self.occ_by[label][bin_idx] += 1.0
            self.explored_by[label][bin_idx] = True

        # 1. structural / DNA-binding proteins occupy sites PERSISTENTLY, turning
        #    over a fraction to new positions each step (gradual exploration). The
        #    full Fig-4E panel: condensin (SMC), SSB, the topoisomerases (GyrAB,
        #    Topo IV), and the replicative helicase/clamp (DnaB, DnaN).
        cfg = self.config
        place(self.oriC, "DnaA")
        self._t_struct += interval
        ramp = 1.0 - np.exp(-self._t_struct / max(float(cfg["struct_bind_tau_s"]), 1.0))
        turnover = float(cfg["struct_turnover_per_s"]) * interval
        # A collision is a BINDING event — a protein landing on an already-occupied
        # site — not persistent occupancy. So only NEWLY-bound sites (this step's
        # rebinds + growth) record collisions; sites that merely persist do not
        # (else the dense SMC coat would flood the matrix every step and RNA pol
        # would stop being the dominant cause, contra Fig 3F's ~84%).
        for label, n_max in (("SMC", cfg["n_smc"]), ("SSB", cfg["n_ssb"]),
                             ("GyrAB", cfg["n_gyrase"]), ("Topo IV", cfg["n_topo"]),
                             ("DnaB", cfg["n_dnab"]), ("DnaN", cfg["n_dnan"])):
            n_now = max(1, int(round(int(n_max) * ramp)))         # progressively-bound count
            sites = self._struct_sites.get(label)
            new_idx = set()
            if sites is None:
                sites = self._rng.integers(0, nb, n_now)
                new_idx = set(range(len(sites)))                  # initial binding
            else:
                k = int(round(min(1.0, turnover) * len(sites)))  # rebind a fraction
                if k > 0:
                    idx = self._rng.choice(len(sites), k, replace=False)
                    sites = sites.copy()
                    sites[idx] = self._rng.integers(0, nb, k)
                    new_idx.update(int(i) for i in idx)
                if n_now > len(sites):                            # grow toward n_max
                    old = len(sites)
                    sites = np.concatenate([sites, self._rng.integers(0, nb, n_now - old)])
                    new_idx.update(range(old, len(sites)))
            self._struct_sites[label] = sites
            frac = float(cfg["struct_collision_frac"])
            rec_idx = {j for j in new_idx if self._rng.random() < frac}
            for j, b in enumerate(sites):
                place(b, label, record=(j in rec_idx))
        # transcription factors at a few promoters (Fur / GntR / HrcA / LuxR) —
        # re-placed (a binding event) each step, so they record collisions.
        n_tf = min(int(cfg["n_tf"]), len(self.gene_bins))
        if n_tf > 0:
            for j, b in enumerate(self._rng.choice(self.gene_bins, n_tf, replace=False)):
                place(b, self._tf_labels[j % len(self._tf_labels)])

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
            released = False
            # a moving RNA pol displaces every bound protein along the bins it
            # sweeps this step (the dominant collision source, Fig 3F) — except
            # the replication fork, which displaces the RNA pol instead.
            stalled_at = None
            for b in range(pos + 1, min(newpos, end) + 1):
                bb = b % nb
                if lesioned and bb in lesioned:
                    # a DNA lesion stalls the elongating polymerase (transcription
                    # blocked at the damaged site) — a real damage↔dynamics coupling.
                    self._collide("RNA Pol", "lesion")
                    stalled_at = b - 1  # halts just before the lesion
                    break
                occ_label = occupant.get(bb)
                if occ_label is None or occ_label == "RNA Pol":
                    continue
                if occ_label == "DNA Pol":
                    self._collide("DNA Pol", "RNA Pol")  # RNA pol displaced by the fork
                    released = True
                    break
                self._collide("RNA Pol", occ_label)      # RNA pol displaces the bound protein
            if stalled_at is not None:
                # blocked by a lesion: hold position just before the damaged site
                # and remain bound (keeps elongating once the site is repaired).
                pol[0] = stalled_at
                place(stalled_at % nb, "RNA Pol", record=False)
                self.explored_rnap[stalled_at % nb] = True
                still.append(pol)
                continue
            pol[0] = newpos
            if not released and newpos <= end:
                place(newpos, "RNA Pol", record=False)   # collisions already counted along the path
                self.explored_rnap[newpos % nb] = True
                still.append(pol)
            # else: reached end / displaced -> released (slot frees for a new initiation)
        self.rna_pols = still

        # 3. DNA polymerase: two replisomes initiate at oriC and advance OUTWARD,
        #    symmetrically, to terC (Fig 4D). Positions are reported as signed
        #    offsets from oriC so the space–time plot is a V (oriC-centred), and
        #    each fork clamps once it reaches terC (that arm is fully replicated).
        dna_report = []
        if state.get("replication_active", 0.0) >= 0.5:
            advanced = min(self.terC, (self.dna_off or 0) + step_bins_dna)
            # a lesion on either fork's path stalls replication until it is repaired
            # (the fork cannot pass a damaged site) — empty lesion set at baseline
            # leaves this unchanged.
            if lesioned:
                for step in range((self.dna_off or 0) + 1, advanced + 1):
                    if (step % nb) in lesioned or ((-step) % nb) in lesioned:
                        self._collide("DNA Pol", "lesion")
                        advanced = step - 1
                        break
            self.dna_off = advanced
            for sign in (+1, -1):
                bb = int(sign * self.dna_off) % nb
                place(bb, "DNA Pol")                     # fork collides with whatever it meets
                self.explored_dnap[bb] = True
                dna_report.append(float(sign * self.dna_off))

        n_coll = float(sum(self.collisions.values()))
        return {
            "fraction_explored": float(self.explored_any.mean()),
            "dna_binding_density": float(self.cur_bound.mean()),
            "n_collisions": n_coll,
            "percent_rnap_explored": float(self.explored_rnap.mean()),
            "percent_dnap_explored": float(self.explored_dnap.mean()),
            "occupancy": {str(i): float(v) for i, v in enumerate(self.occ) if v > 0},
            # signed offsets from oriC (bins): oriC=0 centre, terC=±nb/2
            "rna_pol_positions": [float(self._signed(p[0])) for p in self.rna_pols],
            "dna_pol_positions": dna_report,
            "collisions": dict(self.collisions),
        }
