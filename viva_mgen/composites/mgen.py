"""The single, reusable *Mycoplasma genitalium* whole-cell composite.

One integrated cell — FBA metabolism + mass/growth + stochastic
transcription/translation/decay + dNTP-buffered replication — wired through
shared cell-variable stores. Every figure-study reuses THIS composite, differing
only in the parameters it passes (gene disruption, kinetic scaling, initial DnaA/
dNTP, nutrient level) and the observables it measures. This mirrors the Karr 2012
model's own design: one whole-cell model, many analyses.
"""

from __future__ import annotations

from process_bigraph.composite_generator import composite_generator

from ..expression_defaults import DEFAULT_GENES


def _stores(nutrient_scale, initial_dnaA, initial_dntp):
    # rna_counts / protein_counts pre-seeded with the panel gene keys (a
    # map[float] store accumulates deltas only into existing keys).
    return {
        # metabolism / growth
        "nutrient_scale": nutrient_scale,
        "growth_rate": 0.0, "growth_fraction": 1.0,
        "atp_production": 0.0, "gtp_production": 0.0, "feasible": 1.0,
        # mass / geometry
        "mass": 3.93, "mass_fractions": {}, "volume": 0.0, "division": 0.0,
        # expression
        "ntp": 1e8, "gtp": 1e8, "rna_pol": 100.0,
        "rna_counts": {g: 0.0 for g in DEFAULT_GENES},
        "protein_counts": {g: 0.0 for g in DEFAULT_GENES},
        # replication / cell cycle
        "dntp_synthesis_scale": 1.0,
        "replicated_fraction": 0.0, "dntp_pool": 0.0, "dnaA_complex": 0.0,
        "phase_code": 0.0, "chromosome_copy": 1.0,
        "initiation_duration": 0.0, "replication_duration": 0.0,
        "dntp_at_replication_start": 0.0,
    }


def build_mgen(core=None, *, nutrient_scale=1.0, disrupted_genes=None,
               reaction_bound_scale=None, initial_dnaA=0.0, initial_dntp=0.0,
               seed=0, interval=1.0):
    """Return the integrated M. genitalium composite document.

    Parameters drive the figure-specific behavior; every study reuses this one
    composite. ``disrupted_genes`` (Fig 6) and ``reaction_bound_scale`` (Fig 7)
    tune metabolism; ``initial_dnaA`` / ``initial_dntp`` (Fig 4) tune the cell
    cycle; ``nutrient_scale`` (Fig 2/5) tunes growth.
    """
    disrupted_genes = list(disrupted_genes or [])
    reaction_bound_scale = dict(reaction_bound_scale or {})

    def proc(address, config, inputs, outputs):
        return {"_type": "process", "address": f"local:{address}", "config": config,
                "interval": interval, "inputs": inputs, "outputs": outputs}

    doc = {
        "stores": _stores(nutrient_scale, initial_dnaA, initial_dntp),
        "metabolism": proc("MetabolismFbaReproductionProcess",
            {"disrupted_genes": disrupted_genes, "reaction_bound_scale": reaction_bound_scale},
            {"nutrient_scale": ["stores", "nutrient_scale"]},
            {"growth_rate": ["stores", "growth_rate"],
             "growth_fraction": ["stores", "growth_fraction"],
             "atp_production": ["stores", "atp_production"],
             "gtp_production": ["stores", "gtp_production"],
             "feasible": ["stores", "feasible"]}),
        "mass": proc("MassGrowthReproductionProcess", {},
            {"growth_fraction": ["stores", "growth_fraction"], "mass": ["stores", "mass"]},
            {"mass": ["stores", "mass"], "mass_fractions": ["stores", "mass_fractions"],
             "volume": ["stores", "volume"], "division": ["stores", "division"]}),
        "transcription": proc("TranscriptionReproductionProcess", {"seed": seed},
            {"ntp": ["stores", "ntp"], "rna_pol": ["stores", "rna_pol"]},
            {"rna_counts": ["stores", "rna_counts"], "ntp": ["stores", "ntp"]}),
        "translation": proc("TranslationReproductionProcess", {"seed": seed + 1},
            {"rna_counts": ["stores", "rna_counts"], "gtp": ["stores", "gtp"]},
            {"protein_counts": ["stores", "protein_counts"], "gtp": ["stores", "gtp"]}),
        "rna_decay": proc("RnaDecayReproductionProcess", {},
            {"rna_counts": ["stores", "rna_counts"]}, {"rna_counts": ["stores", "rna_counts"]}),
        "protein_decay": proc("ProteinDecayReproductionProcess", {},
            {"protein_counts": ["stores", "protein_counts"]},
            {"protein_counts": ["stores", "protein_counts"]}),
        "replication": proc("ReplicationReproductionProcess",
            {"initial_dnaA": initial_dnaA, "initial_dntp": initial_dntp, "seed": seed},
            {"dntp_synthesis_scale": ["stores", "dntp_synthesis_scale"]},
            {"replicated_fraction": ["stores", "replicated_fraction"],
             "dntp_pool": ["stores", "dntp_pool"], "dnaA_complex": ["stores", "dnaA_complex"],
             "phase_code": ["stores", "phase_code"], "chromosome_copy": ["stores", "chromosome_copy"],
             "initiation_duration": ["stores", "initiation_duration"],
             "replication_duration": ["stores", "replication_duration"],
             "dntp_at_replication_start": ["stores", "dntp_at_replication_start"]}),
    }
    emit = {"mass": "float", "mass_fractions": "map[float]",
            "growth_fraction": "float", "growth_rate": "float",
            "volume": "float", "division": "float",
            "atp_production": "float", "gtp_production": "float", "feasible": "float",
            "rna_counts": "map[float]", "protein_counts": "map[float]",
            "replicated_fraction": "float", "dntp_pool": "float",
            "dnaA_complex": "float", "phase_code": "float", "chromosome_copy": "float"}
    doc["emitter"] = {"_type": "step", "address": "local:RAMEmitter",
                      "config": {"emit": emit},
                      "inputs": {k: ["stores", k] for k in emit}}
    return doc


@composite_generator(
    name="mycoplasma_genitalium",
    description="The whole Mycoplasma genitalium cell: FBA metabolism + mass/growth + stochastic transcription/translation/decay + dNTP-buffered replication, wired through shared cell-variable stores. Reused by every figure-study.",
    parameters={
        "nutrient_scale": {"type": "float", "default": 1.0,
                           "description": "Carbon-uptake availability multiplier (Fig 2/5)"},
        "initial_dnaA": {"type": "float", "default": 0.0,
                         "description": "DnaA monomers at birth (Fig 4 cell cycle)"},
        "initial_dntp": {"type": "float", "default": 0.0,
                         "description": "dNTP molecules at birth (Fig 4)"},
        "seed": {"type": "integer", "default": 0, "description": "Stochastic seed"},
    },
    emitters=[{"address": "local:ParquetEmitter",
               "paths": ["stores/mass", "stores/growth_fraction", "stores/growth_rate",
                         "stores/volume", "stores/atp_production", "stores/gtp_production",
                         "stores/replicated_fraction", "stores/dntp_pool"]}],
)
def mycoplasma_genitalium(core=None, *, nutrient_scale=1.0, initial_dnaA=0.0,
                          initial_dntp=0.0, seed=0):
    return build_mgen(core, nutrient_scale=nutrient_scale, initial_dnaA=initial_dnaA,
                      initial_dntp=initial_dntp, seed=seed)
