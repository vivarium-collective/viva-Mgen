"""Cell mass & growth submodel — clean-room reproduction.

Reproduces the mass/geometry accounting of the Karr 2012 whole-cell model
(``CellMass`` / ``CellGeometry``): it integrates the metabolism submodel's
calibrated growth into total dry mass, partitions mass into the fitted
dry-weight fractions (DNA/RNA/protein/lipid/…), derives cell volume, and flags
division when mass has doubled.

Growth is calibrated so that, unperturbed (``growth_fraction`` = 1.0), the cell
doubles in the model's fitted cell-cycle length (``CELL_CYCLE_LENGTH_S`` ≈ 9 h),
matching the measured *M. genitalium* doubling time reproduced in Fig 2A/B.
"""

from __future__ import annotations

import math

from process_bigraph import Process

from .. import constants as C
from ..mass_composition import emergent_composition


class MassGrowthReproductionProcess(Process):
    """Integrate calibrated growth into cell mass, fractions, volume, division.

    Inputs
    ------
    growth_fraction : float
        Unit-free growth relative to unperturbed (from metabolism). 1.0 == the
        calibrated wild-type rate that doubles the cell in ~9 h.
    mass : float
        Current total dry mass (fg).
    rna_counts, protein_counts : map[float]
        Per-gene RNA/protein copy numbers (read-only sensor for the emergent
        composition below).
    chromosome_copy : float
        Chromosome copy number (read-only sensor).
    atp, gtp, ntp, amino_acid : float
        Free-metabolite pool sizes (read-only sensor, diagnostic only).

    Outputs
    -------
    mass : float
        Dry-mass delta this interval (fg) — additive, composes with any other
        mass contributor.
    mass_fractions : overwrite[map[float]]
        Current dry mass per component (fg), snapshot.
    volume : overwrite[float]
        Current cell volume (fL), snapshot.
    division : overwrite[float]
        1.0 once mass has reached twice the initial dry mass, else 0.0.
    emergent_mass : overwrite[float]
        Σ species×MW over RNA + protein + DNA only (fg) — validated against
        Karr's fitted macromolecular dry-weight fractions.
    emergent_mass_fractions : overwrite[map[float]]
        {RNA, protein, DNA} fractions of emergent_mass, summing to 1 (or all
        zero if emergent_mass is 0). Does NOT include metabolite.
    metabolite_mass : overwrite[float]
        Free-metabolite pool mass (fg) — non-validated diagnostic, kept
        separate from emergent_mass/emergent_mass_fractions.
    """

    description = (
        "Cell mass & growth accounting — reproduction of Karr 2012 CellMass/CellGeometry.\n"
        "Integrates the metabolism submodel's calibrated growth into total dry mass by "
        "exponential growth over each interval:\n"
        "    Δmass = mass · (exp(μ · growth_fraction · Δt) − 1),   μ = ln2 / T_cycle\n"
        "with T_cycle = 32400 s (≈ 9 h) the model's fitted cell-cycle length, so unperturbed "
        "growth (growth_fraction = 1) doubles the cell in ~9 h. Dry mass is partitioned into "
        "the fitted dry-weight fractions (protein 0.620, DNA 0.169, RNA 0.093, lipid 0.057, …); "
        "volume follows from wet mass and cell density; division fires at 2× birth mass.\n"
        "Contract — in: growth_fraction (unit-free growth from metabolism), mass (current dry "
        "fg). out: mass (additive Δ, fg), mass_fractions (per-component fg snapshot), volume "
        "(fL snapshot), division ∈ {0,1}.\n"
        "Fidelity: FAITHFUL growth/composition phenotype — real fitted constants (dry-weight "
        "3.93e-15 g, T_cycle 32400 s, component fractions, cell density). This is the calibrated "
        "growth-law form: mass follows metabolism-driven exponential growth partitioned by the "
        "fitted fractions, rather than being summed emergently from every species count × MW. "
        "In addition, the emergent macromolecular dry-mass composition (Σ species×MW over "
        "RNA/protein/DNA, from the real molecular inventory) is now computed and emitted "
        "(emergent_mass/emergent_mass_fractions) for validation against Karr's fitted fractions; "
        "the free-metabolite pool mass is emitted separately as a non-validated diagnostic "
        "(metabolite_mass), since current pool sizes are allocation-dynamics placeholders rather "
        "than calibrated copy numbers; the growth-law total above still drives division, pending "
        "synthesis calibration that would let the emergent total itself close the loop."
    )

    config_schema = {
        "initial_mass_fg": {"_type": "float", "_default": C.CELL_INITIAL_DRY_WEIGHT_FG},
        "cell_cycle_length_s": {"_type": "float", "_default": C.CELL_CYCLE_LENGTH_S},
        # g/mL == 1e-12 fg-per-fL... using 1.1 g/mL cell density
        "density_g_per_ml": {"_type": "float", "_default": 1.1},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        # calibrated specific growth rate (per second): ln2 / cell-cycle-length
        self._mu = math.log(2.0) / float(self.config["cell_cycle_length_s"])

    def inputs(self):
        return {
            "growth_fraction": "float", "mass": "float",
            "rna_counts": "map[float]", "protein_counts": "map[float]",
            "chromosome_copy": "float",
            "atp": "float", "gtp": "float", "ntp": "float", "amino_acid": "float",
        }

    def outputs(self):
        return {
            "mass": "float",
            "mass_fractions": "overwrite[map[float]]",
            "volume": "overwrite[float]",
            "division": "overwrite[float]",
            "emergent_mass": "overwrite[float]",
            "emergent_mass_fractions": "overwrite[map[float]]",
            "metabolite_mass": "overwrite[float]",
        }

    def initial_state(self):
        return {
            "growth_fraction": 1.0, "mass": float(self.config["initial_mass_fg"]),
            "rna_counts": {}, "protein_counts": {}, "chromosome_copy": 1.0,
            "atp": 0.0, "gtp": 0.0, "ntp": 0.0, "amino_acid": 0.0,
        }

    def _volume_fl(self, dry_mass_fg: float) -> float:
        wet_mass_fg = dry_mass_fg / (1.0 - C.FRACTION_WET_WEIGHT)
        # fg -> g (1e-15), density g/mL, volume mL -> fL (1e12)
        wet_mass_g = wet_mass_fg * 1e-15
        vol_ml = wet_mass_g / float(self.config["density_g_per_ml"])
        return vol_ml * 1e12  # fL

    def update(self, state, interval):
        mass = float(state.get("mass", self.config["initial_mass_fg"]))
        gf = float(state.get("growth_fraction", 1.0))
        # exponential growth over the interval
        d_mass = mass * (math.exp(self._mu * gf * interval) - 1.0)
        new_mass = mass + d_mass
        fractions = {k: new_mass * v for k, v in C.DRY_MASS_FRACTIONS.items()}
        divided = 1.0 if new_mass >= 2.0 * float(self.config["initial_mass_fg"]) else 0.0

        # Emergent macromolecular dry-mass composition (RNA + protein + DNA),
        # computed from the real molecular inventory for validation against
        # Karr's fitted dry-weight fractions. The free-metabolite pool is
        # EXCLUDED from the validated total/fractions — the pool sizes are
        # allocation-dynamics placeholders, not calibrated copy numbers, and
        # would distort the macromolecule-dominated comparison — and is
        # instead emitted separately as a non-validated diagnostic.
        comp = emergent_composition(
            state.get("rna_counts", {}), state.get("protein_counts", {}),
            float(state.get("chromosome_copy", 1.0) or 1.0),
            {k: float(state.get(k, 0.0) or 0.0) for k in ("atp", "gtp", "ntp", "amino_acid")},
        )
        macro_total = comp["RNA"] + comp["protein"] + comp["DNA"]
        macro_fractions = (
            {k: comp[k] / macro_total for k in ("RNA", "protein", "DNA")}
            if macro_total > 0 else {k: 0.0 for k in ("RNA", "protein", "DNA")}
        )

        return {
            "mass": d_mass,
            "mass_fractions": fractions,
            "volume": self._volume_fl(new_mass),
            "division": divided,
            "emergent_mass": macro_total * 1e15,
            "emergent_mass_fractions": macro_fractions,
            "metabolite_mass": comp["metabolite"] * 1e15,
        }
