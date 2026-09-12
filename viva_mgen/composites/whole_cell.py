"""Composite generators for the integrated whole-cell studies (Fig 1, 2, 3, 5).

All four share the same reduced integrated cell — metabolism (FBA) + mass/growth
+ stochastic transcription/translation/decay — wired through shared bigraph
stores, differing only in which observables they foreground.
"""

from __future__ import annotations

from process_bigraph.composite_generator import composite_generator


def _base_stores():
    return {
        "nutrient_scale": 1.0,
        "growth_fraction": 1.0,
        "growth_rate": 0.0,
        "atp_production": 0.0,
        "gtp_production": 0.0,
        "feasible": 1.0,
        "mass": 3.93,
        "mass_fractions": {},
        "volume": 0.0,
        "division": 0.0,
        "ntp": 1e8,
        "gtp": 1e8,
        "rna_pol": 100.0,
        "rna_counts": {},
        "protein_counts": {},
    }


def _metabolism_node():
    return {
        "_type": "process",
        "address": "local:MetabolismFbaReproductionProcess",
        "config": {},
        "interval": 1.0,
        "inputs": {"nutrient_scale": ["stores", "nutrient_scale"]},
        "outputs": {
            "growth_rate": ["stores", "growth_rate"],
            "growth_fraction": ["stores", "growth_fraction"],
            "atp_production": ["stores", "atp_production"],
            "gtp_production": ["stores", "gtp_production"],
            "feasible": ["stores", "feasible"],
        },
    }


def _mass_node():
    return {
        "_type": "process",
        "address": "local:MassGrowthReproductionProcess",
        "config": {},
        "interval": 1.0,
        "inputs": {"growth_fraction": ["stores", "growth_fraction"], "mass": ["stores", "mass"]},
        "outputs": {
            "mass": ["stores", "mass"],
            "mass_fractions": ["stores", "mass_fractions"],
            "volume": ["stores", "volume"],
            "division": ["stores", "division"],
        },
    }


def _expression_nodes():
    return {
        "transcription": {
            "_type": "process",
            "address": "local:TranscriptionReproductionProcess",
            "config": {},
            "interval": 1.0,
            "inputs": {"ntp": ["stores", "ntp"], "rna_pol": ["stores", "rna_pol"]},
            "outputs": {"rna_counts": ["stores", "rna_counts"], "ntp": ["stores", "ntp"]},
        },
        "translation": {
            "_type": "process",
            "address": "local:TranslationReproductionProcess",
            "config": {},
            "interval": 1.0,
            "inputs": {"rna_counts": ["stores", "rna_counts"], "gtp": ["stores", "gtp"]},
            "outputs": {"protein_counts": ["stores", "protein_counts"], "gtp": ["stores", "gtp"]},
        },
        "rna_decay": {
            "_type": "process",
            "address": "local:RnaDecayReproductionProcess",
            "config": {},
            "interval": 1.0,
            "inputs": {"rna_counts": ["stores", "rna_counts"]},
            "outputs": {"rna_counts": ["stores", "rna_counts"]},
        },
        "protein_decay": {
            "_type": "process",
            "address": "local:ProteinDecayReproductionProcess",
            "config": {},
            "interval": 1.0,
            "inputs": {"protein_counts": ["stores", "protein_counts"]},
            "outputs": {"protein_counts": ["stores", "protein_counts"]},
        },
    }


def _emitter(emit):
    return {
        "_type": "step",
        "address": "local:RAMEmitter",
        "config": {"emit": emit},
        "inputs": {k: ["stores", k] for k in emit},
    }


def _integrated_document(nutrient_scale=1.0):
    stores = _base_stores()
    stores["nutrient_scale"] = nutrient_scale
    doc = {"stores": stores, "metabolism": _metabolism_node(), "mass": _mass_node()}
    doc.update(_expression_nodes())
    doc["emitter"] = _emitter({
        "mass": "float", "growth_fraction": "float", "growth_rate": "float",
        "volume": "float", "division": "float",
        "atp_production": "float", "gtp_production": "float",
        "rna_counts": "map[float]", "protein_counts": "map[float]",
    })
    return doc


@composite_generator(
    name="fig1_architecture",
    description="Fig 1 — the integrated M. genitalium cell: FBA metabolism + mass + stochastic expression wired through shared cell-variable stores.",
    parameters={"nutrient_scale": {"type": "float", "default": 1.0,
                                   "description": "Carbon-uptake availability multiplier"}},
    emitters=[{"address": "local:ParquetEmitter",
               "paths": ["stores/mass", "stores/growth_fraction", "stores/rna_counts", "stores/protein_counts"]}],
)
def fig1_architecture(core=None, *, nutrient_scale=1.0):
    return _integrated_document(nutrient_scale)


@composite_generator(
    name="fig2_growth",
    description="Fig 2 — cell growth: mass accumulation, doubling time (~9 h), and dry-mass composition from FBA-driven growth.",
    parameters={"nutrient_scale": {"type": "float", "default": 1.0,
                                   "description": "Carbon-uptake availability multiplier"}},
    emitters=[{"address": "local:ParquetEmitter",
               "paths": ["stores/mass", "stores/growth_fraction", "stores/mass_fractions", "stores/volume"]}],
)
def fig2_growth(core=None, *, nutrient_scale=1.0):
    doc = {"stores": {**_base_stores(), "nutrient_scale": nutrient_scale},
           "metabolism": _metabolism_node(), "mass": _mass_node()}
    doc["emitter"] = _emitter({"mass": "float", "growth_fraction": "float",
                               "growth_rate": "float", "volume": "float",
                               "mass_fractions": "map[float]", "division": "float"})
    return doc


@composite_generator(
    name="fig3_expression",
    description="Fig 3 — gene expression dynamics: stochastic transcription/translation/decay producing bursty mRNA and accumulating protein across the representative gene panel.",
    parameters={"seed": {"type": "integer", "default": 0,
                         "description": "Random seed for the stochastic expression"}},
    emitters=[{"address": "local:ParquetEmitter",
               "paths": ["stores/rna_counts", "stores/protein_counts"]}],
)
def fig3_expression(core=None, *, seed=0):
    stores = _base_stores()
    doc = {"stores": stores}
    nodes = _expression_nodes()
    nodes["transcription"]["config"] = {"seed": seed}
    nodes["translation"]["config"] = {"seed": seed + 1}
    doc.update(nodes)
    doc["emitter"] = _emitter({"rna_counts": "map[float]", "protein_counts": "map[float]"})
    return doc


@composite_generator(
    name="fig5_energy",
    description="Fig 5 — global energy allocation: ATP/GTP production from FBA alongside the GTP/NTP consumption of transcription and translation.",
    parameters={"nutrient_scale": {"type": "float", "default": 1.0,
                                   "description": "Carbon-uptake availability multiplier"}},
    emitters=[{"address": "local:ParquetEmitter",
               "paths": ["stores/atp_production", "stores/gtp_production", "stores/growth_rate"]}],
)
def fig5_energy(core=None, *, nutrient_scale=1.0):
    return _integrated_document(nutrient_scale)
