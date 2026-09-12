"""Composite generators for the gene-disruption studies (Fig 6, 7).

Both drive the FBA metabolism process; the full single-gene-deletion scan
(Fig 6) and the kcat sweep (Fig 7) live in the studies' ``sims/run.py``.
"""

from __future__ import annotations

from process_bigraph.composite_generator import composite_generator


def _metabolism_doc(config):
    return {
        "stores": {"nutrient_scale": 1.0, "growth_rate": 0.0, "growth_fraction": 1.0,
                   "atp_production": 0.0, "gtp_production": 0.0, "feasible": 1.0},
        "metabolism": {
            "_type": "process",
            "address": "local:MetabolismFbaReproductionProcess",
            "config": config,
            "interval": 1.0,
            "inputs": {"nutrient_scale": ["stores", "nutrient_scale"]},
            "outputs": {
                "growth_rate": ["stores", "growth_rate"],
                "growth_fraction": ["stores", "growth_fraction"],
                "atp_production": ["stores", "atp_production"],
                "gtp_production": ["stores", "gtp_production"],
                "feasible": ["stores", "feasible"],
            },
        },
        "emitter": {
            "_type": "step",
            "address": "local:RAMEmitter",
            "config": {"emit": {"growth_rate": "float", "growth_fraction": "float",
                                "feasible": "float"}},
            "inputs": {k: ["stores", k] for k in ("growth_rate", "growth_fraction", "feasible")},
        },
    }


@composite_generator(
    name="fig6_gene_essentiality",
    description="Fig 6 — single-gene-disruption phenotypes: knock out a metabolic gene and measure growth via FBA (essential if growth collapses).",
    parameters={"disrupted_gene": {"type": "string", "default": "",
                                   "description": "Gene id to knock out (e.g. MG_023); empty = wild-type"}},
    emitters=[{"address": "local:ParquetEmitter",
               "paths": ["stores/growth_fraction", "stores/feasible"]}],
)
def fig6_gene_essentiality(core=None, *, disrupted_gene=""):
    genes = [disrupted_gene] if disrupted_gene else []
    return _metabolism_doc({"disrupted_genes": genes})


@composite_generator(
    name="fig7_kinetic_parameters",
    description="Fig 7 — kinetic-parameter characterization: scale a target reaction's flux bound (a kcat proxy) and measure the effect on growth.",
    parameters={
        "reaction_id": {"type": "string", "default": "PGK",
                        "description": "Reaction whose flux bound (kcat proxy) to scale"},
        "bound_scale": {"type": "float", "default": 1.0,
                        "description": "Multiplier on the reaction's flux bound"},
    },
    emitters=[{"address": "local:ParquetEmitter",
               "paths": ["stores/growth_fraction", "stores/growth_rate"]}],
)
def fig7_kinetic_parameters(core=None, *, reaction_id="PGK", bound_scale=1.0):
    return _metabolism_doc({"reaction_bound_scale": {reaction_id: bound_scale}})
