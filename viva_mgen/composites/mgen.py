"""The single, reusable *Mycoplasma genitalium* whole-cell composite.

One integrated cell wiring ALL ~28 Karr-2012 submodel reproductions
(metabolism, mass/growth, the full transcription→RNA-processing→modification→
aminoacylation chain, the translation→processing→translocation→folding→
modification→complexation→ribosome/terminal-organelle-assembly chain, the DNA
replication-initiation→replication→supercoiling→condensation→segregation→damage→
repair chain, and FtsZ→cytokinesis + host interaction) through shared
cell-variable stores. Every figure-study reuses THIS composite, differing only
in the parameters it passes and the observables it measures — mirroring the
Karr 2012 model's own "one whole-cell model, many analyses" design.
"""

from __future__ import annotations

from process_bigraph.composite_generator import composite_generator

from ..expression_defaults import DEFAULT_GENES
from ..processes import all_process_classes

# class name -> composite node name
_NODE_NAMES = {
    "MetabolismFbaReproductionProcess": "metabolism",
    "MassGrowthReproductionProcess": "mass",
    "TranscriptionReproductionProcess": "transcription",
    "TranslationReproductionProcess": "translation",
    "RnaDecayReproductionProcess": "rna_decay",
    "ProteinDecayReproductionProcess": "protein_decay",
    "ReplicationReproductionProcess": "replication",
    "ReplicationInitiationReproductionProcess": "replication_initiation",
    "DNASupercoilingReproductionProcess": "supercoiling",
    "ChromosomeCondensationReproductionProcess": "condensation",
    "ChromosomeSegregationReproductionProcess": "segregation",
    "DNADamageReproductionProcess": "dna_damage",
    "DNARepairReproductionProcess": "dna_repair",
    "TranscriptionalRegulationReproductionProcess": "transcriptional_regulation",
    "RNAProcessingReproductionProcess": "rna_processing",
    "RNAModificationReproductionProcess": "rna_modification",
    "TRNAAminoacylationReproductionProcess": "trna_aminoacylation",
    "ProteinProcessingIReproductionProcess": "protein_processing_i",
    "ProteinTranslocationReproductionProcess": "protein_translocation",
    "ProteinProcessingIIReproductionProcess": "protein_processing_ii",
    "ProteinFoldingReproductionProcess": "protein_folding",
    "ProteinModificationReproductionProcess": "protein_modification",
    "ProteinActivationReproductionProcess": "protein_activation",
    "MacromolecularComplexationReproductionProcess": "complexation",
    "RibosomeAssemblyReproductionProcess": "ribosome_assembly",
    "TerminalOrganelleAssemblyReproductionProcess": "terminal_organelle_assembly",
    "FtsZPolymerizationReproductionProcess": "ftsz",
    "CytokinesisReproductionProcess": "cytokinesis",
    "HostInteractionReproductionProcess": "host_interaction",
}

# port -> store overrides (disambiguate two processes' identically-named ports)
_PORT_STORE = {
    ("RNAModificationReproductionProcess", "modification_enzyme"): "rna_modification_enzyme",
    ("ProteinModificationReproductionProcess", "modification_enzyme"): "protein_modification_enzyme",
}

# store keys that hold map[float] (per-gene) values — pre-seeded so additive
# deltas land (a map[float] store drops deltas to keys it doesn't already have)
_MAP_STORES = {
    "rna_counts", "protein_counts", "mass_fractions", "fold_change", "nascent_rna",
    "mature_rna", "unmodified_rna", "modified_rna", "tf_activity", "nascent",
    "process_i_done", "translocated", "processed_ii", "unfolded", "folded",
    "unmodified", "modified", "monomers", "complexes", "active_fraction",
    "free_trna", "aminoacylated_trna", "rprotein_counts", "rrna_counts",
    "adhesin_proteins",
}

# initial float-store values (resources/enzymes/setpoints); unlisted floats -> 0.0
_FLOAT_INIT = {
    "atp": 1e9, "gtp": 1e9, "ntp": 1e8, "amino_acid": 1e8,
    "nutrient_scale": 1.0, "rna_pol": 100.0, "dntp_synthesis_scale": 1.0,
    "mass": 3.93, "growth_fraction": 1.0, "feasible": 1.0, "chromosome_copy": 1.0,
    "dnaA_free": 50.0, "gyrase": 20.0, "smc": 30.0, "damaging_agent": 0.0,
    "repair_enzyme": 20.0, "deformylase": 20.0, "translocase": 20.0,
    "signal_peptidase": 20.0, "chaperone_count": 50.0, "rna_modification_enzyme": 20.0,
    "protein_modification_enzyme": 20.0, "regulator": 1.0, "synthetase": 30.0,
    "assembly_factor": 20.0, "ftsz_monomer": 500.0, "septum_diameter": 200.0,
}

# observables the RAM emitter captures (study-measured + integration readouts)
_EMIT = {
    "mass": "float", "mass_fractions": "map[float]", "growth_fraction": "float",
    "growth_rate": "float", "volume": "float", "division": "float",
    "atp_production": "float", "gtp_production": "float", "feasible": "float",
    "rna_counts": "map[float]", "protein_counts": "map[float]",
    "replicated_fraction": "float", "dntp_pool": "float", "dnaA_complex": "float",
    "phase_code": "float", "chromosome_copy": "float",
    "condensed_fraction": "float", "segregated_fraction": "float",
    "superhelical_density": "float", "lesions": "float", "complex_size": "float",
    "initiation_ready": "float", "ftsz_ring_filaments": "float",
    "septum_diameter": "float", "contraction_progress": "float", "divided": "float",
    "terminal_organelle_fraction": "float", "adhesion_strength": "float",
    "ribosome_30S": "float", "ribosome_50S": "float",
}


def _store_key(cls_name, port):
    return _PORT_STORE.get((cls_name, port), port)


def build_mgen(core=None, *, nutrient_scale=1.0, disrupted_genes=None,
               reaction_bound_scale=None, initial_dnaA=0.0, initial_dntp=0.0,
               seed=0, interval=1.0):
    """Return the integrated 28-submodel M. genitalium composite document."""
    classes = all_process_classes()
    configs = {
        "metabolism": {"disrupted_genes": list(disrupted_genes or []),
                       "reaction_bound_scale": dict(reaction_bound_scale or {})},
        "transcription": {"seed": seed},
        "translation": {"seed": seed + 1},
        "replication": {"initial_dnaA": initial_dnaA, "initial_dntp": initial_dntp, "seed": seed},
    }

    stores, doc, used = {}, {}, set()

    def store(key):
        used.add(key)
        if key not in stores:
            stores[key] = ({g: 0.0 for g in DEFAULT_GENES} if key in _MAP_STORES
                           else _FLOAT_INIT.get(key, 0.0))
        return ["stores", key]

    for cls_name, node in _NODE_NAMES.items():
        cls = classes[cls_name]
        proc = cls(config={}, core=core)
        inputs = {p: store(_store_key(cls_name, p)) for p in proc.inputs()}
        outputs = {p: store(_store_key(cls_name, p)) for p in proc.outputs()}
        doc[node] = {"_type": "process", "address": f"local:{cls_name}",
                     "config": configs.get(node, {}), "interval": interval,
                     "inputs": inputs, "outputs": outputs}

    stores["nutrient_scale"] = nutrient_scale  # honor the param override
    doc["stores"] = stores
    emit = {k: v for k, v in _EMIT.items() if k in used}
    doc["emitter"] = {"_type": "step", "address": "local:RAMEmitter",
                      "config": {"emit": emit},
                      "inputs": {k: ["stores", k] for k in emit}}
    return doc


@composite_generator(
    name="mycoplasma_genitalium",
    description="The whole Mycoplasma genitalium cell: all ~28 Karr-2012 submodel reproductions (metabolism, mass, expression, RNA/protein maturation, DNA replication/repair, cytokinesis, host interaction) wired through shared cell-variable stores. Reused by every figure-study.",
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
