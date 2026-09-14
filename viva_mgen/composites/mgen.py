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
    "ChromosomeDynamicsReproductionProcess": "chromosome",
}

# port -> store overrides (disambiguate two processes' identically-named ports)
_PORT_STORE = {
    ("RNAModificationReproductionProcess", "modification_enzyme"): "rna_modification_enzyme",
    ("ProteinModificationReproductionProcess", "modification_enzyme"): "protein_modification_enzyme",
    ("ChromosomeDynamicsReproductionProcess", "rna_polymerase"): "rna_pol",
}

# store keys that hold map[float] (per-gene) values — pre-seeded so additive
# deltas land (a map[float] store drops deltas to keys it doesn't already have)
_MAP_STORES = {
    "rna_counts", "protein_counts", "mass_fractions", "fold_change", "nascent_rna",
    "mature_rna", "unmodified_rna", "modified_rna", "tf_activity", "nascent",
    "process_i_done", "translocated", "processed_ii", "unfolded", "folded",
    "unmodified", "modified", "monomers", "complexes", "active_fraction",
    "free_trna", "aminoacylated_trna", "rprotein_counts", "rrna_counts",
    "adhesin_proteins", "occupancy", "collisions",
}

# store keys that hold list[float] (chromosome polymerase positions)
_LIST_STORES = {"rna_pol_positions", "dna_pol_positions"}

# The whole cell is one root store tree named ``cell`` (not the generic
# "stores"), organised into biologically meaningful compartments/pools. Every
# cell variable lives under ``cell/<compartment>/<name>`` so the loom shows an
# intuitive, navigable hierarchy instead of one flat 87-node row. Each process
# still wires to the SAME leaf path everywhere, so shared state stays shared.
_ROOT = "cell"
_STORE_GROUP = {
    # whole-cell physiology / growth state
    "mass": "physiology", "volume": "physiology", "growth_rate": "physiology",
    "growth_fraction": "physiology", "phase_code": "physiology",
    # metabolism — energy carriers, precursor pools, biomass, FBA status
    "atp": "metabolism", "gtp": "metabolism", "ntp": "metabolism",
    "amino_acid": "metabolism", "dntp_pool": "metabolism",
    "dntp_at_replication_start": "metabolism", "dntp_synthesis_scale": "metabolism",
    "nutrient_scale": "metabolism", "atp_production": "metabolism",
    "gtp_production": "metabolism", "feasible": "metabolism",
    "mass_fractions": "metabolism",
    # genome — chromosome replication, structure, maintenance, DNA-binding (Fig 3)
    "chromosome_copy": "genome", "replicated_fraction": "genome",
    "replication_active": "genome", "replication_duration": "genome",
    "initiation_duration": "genome", "initiation_ready": "genome",
    "dnaA_free": "genome", "dnaA_complex": "genome", "complex_size": "genome",
    "superhelical_density": "genome", "condensed_fraction": "genome",
    "segregated_fraction": "genome", "lesions": "genome",
    "damaging_agent": "genome", "repair_enzyme": "genome", "gyrase": "genome",
    "smc": "genome", "occupancy": "genome", "collisions": "genome",
    "n_collisions": "genome", "fraction_explored": "genome",
    "dna_binding_density": "genome", "percent_rnap_explored": "genome",
    "percent_dnap_explored": "genome", "rna_pol_positions": "genome",
    "dna_pol_positions": "genome",
    # transcriptome — RNA synthesis/processing/modification, tRNA charging, regulation
    "rna_counts": "transcriptome", "nascent_rna": "transcriptome",
    "mature_rna": "transcriptome", "unmodified_rna": "transcriptome",
    "modified_rna": "transcriptome", "rna_pol": "transcriptome",
    "rna_modification_enzyme": "transcriptome", "tf_activity": "transcriptome",
    "fold_change": "transcriptome", "regulator": "transcriptome",
    "free_trna": "transcriptome", "aminoacylated_trna": "transcriptome",
    "synthetase": "transcriptome",
    # proteome — protein maturation chain and complexes
    "protein_counts": "proteome", "nascent": "proteome",
    "process_i_done": "proteome", "processed_ii": "proteome",
    "translocated": "proteome", "unfolded": "proteome", "folded": "proteome",
    "unmodified": "proteome", "modified": "proteome", "monomers": "proteome",
    "complexes": "proteome", "active_fraction": "proteome",
    "chaperone_count": "proteome", "deformylase": "proteome",
    "aminopeptidase": "proteome",
    "translocase": "proteome", "signal_peptidase": "proteome",
    "protein_modification_enzyme": "proteome",
    # ribosome — assembly of the translation machinery
    "ribosome_30S": "ribosome", "ribosome_50S": "ribosome",
    "rprotein_counts": "ribosome", "rrna_counts": "ribosome",
    "assembly_factor": "ribosome",
    # division — FtsZ ring, septum, cytokinesis
    "division": "division", "divided": "division", "ftsz_monomer": "division",
    "ftsz_ring_filaments": "division", "septum_diameter": "division",
    "contraction_progress": "division",
    # envelope — cell surface, terminal organelle, adhesion
    "terminal_organelle_fraction": "envelope", "adhesion_strength": "envelope",
    "adhesin_proteins": "envelope",
}


def _store_path(key):
    """``[cell, <compartment>, <key>]`` — fail loud on an unmapped store so a new
    process port can never silently land outside the biological hierarchy."""
    group = _STORE_GROUP.get(key)
    if group is None:
        raise KeyError(f"store {key!r} has no biological compartment in _STORE_GROUP")
    return [_ROOT, group, key]

# initial float-store values (resources/enzymes/setpoints); unlisted floats -> 0.0
_FLOAT_INIT = {
    "atp": 1e9, "gtp": 1e9, "ntp": 1e8, "amino_acid": 1e8,
    "nutrient_scale": 1.0, "rna_pol": 100.0, "dntp_synthesis_scale": 1.0,
    "mass": 3.93, "growth_fraction": 1.0, "feasible": 1.0, "chromosome_copy": 1.0,
    "dnaA_free": 50.0, "gyrase": 20.0, "smc": 30.0, "damaging_agent": 0.0,
    "repair_enzyme": 20.0, "deformylase": 20.0, "aminopeptidase": 20.0, "translocase": 20.0,
    "signal_peptidase": 20.0, "chaperone_count": 50.0, "rna_modification_enzyme": 20.0,
    "protein_modification_enzyme": 20.0, "regulator": 1.0, "synthetase": 30.0,
    "assembly_factor": 20.0, "ftsz_monomer": 500.0, "septum_diameter": 200.0,
    "replication_active": 1.0,
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
    "fraction_explored": "float", "dna_binding_density": "float",
    "n_collisions": "float", "percent_rnap_explored": "float",
    "percent_dnap_explored": "float",
}


def _store_key(cls_name, port):
    return _PORT_STORE.get((cls_name, port), port)


def build_mgen(core=None, *, nutrient_scale=1.0, disrupted_genes=None,
               reaction_bound_scale=None, initial_dnaA=0.0, initial_dntp=0.0,
               seed=0, interval=1.0):
    """Return the integrated 28-submodel M. genitalium composite document."""
    if core is None:
        # port introspection below instantiates each process, which needs a core;
        # the generator-build path calls us with core=None, so allocate one.
        from ..core import build_core
        core = build_core()
    classes = all_process_classes()
    configs = {
        "metabolism": {"disrupted_genes": list(disrupted_genes or []),
                       "reaction_bound_scale": dict(reaction_bound_scale or {})},
        "transcription": {"seed": seed},
        "translation": {"seed": seed + 1},
        "replication": {"initial_dnaA": initial_dnaA, "initial_dntp": initial_dntp, "seed": seed},
        "chromosome": {"seed": seed},
    }

    cell, doc, used = {}, {}, set()

    def _init_value(key):
        if key in _LIST_STORES:
            return []
        if key in _MAP_STORES:
            return ({} if key in ("occupancy", "collisions", "mass_fractions",
                                  "fold_change", "active_fraction")
                    else {g: 0.0 for g in DEFAULT_GENES})
        return _FLOAT_INIT.get(key, 0.0)

    def store(key):
        used.add(key)
        path = _store_path(key)
        group, leaf = path[1], path[2]
        grp = cell.setdefault(group, {})
        if leaf not in grp:
            grp[leaf] = _init_value(key)
        return path

    for cls_name, node in _NODE_NAMES.items():
        cls = classes[cls_name]
        proc = cls(config={}, core=core)
        inputs = {p: store(_store_key(cls_name, p)) for p in proc.inputs()}
        outputs = {p: store(_store_key(cls_name, p)) for p in proc.outputs()}
        # Dotted address so the dashboard/loom can import the class and show its
        # describe() contract + docstring (bare local:<Class> can't be imported).
        doc[node] = {"_type": "process", "address": f"local:{cls.__module__}.{cls_name}",
                     "config": configs.get(node, {}), "interval": interval,
                     "inputs": inputs, "outputs": outputs}

    # honor the param override (nutrient_scale lives under cell/metabolism)
    cell.setdefault("metabolism", {})["nutrient_scale"] = nutrient_scale
    doc[_ROOT] = cell
    emit = {k: v for k, v in _EMIT.items() if k in used}
    doc["emitter"] = {"_type": "step", "address": "local:RAMEmitter",
                      "config": {"emit": emit},
                      "inputs": {k: _store_path(k) for k in emit}}
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
               "paths": ["cell/physiology/mass", "cell/physiology/growth_fraction",
                         "cell/physiology/growth_rate", "cell/physiology/volume",
                         "cell/metabolism/atp_production", "cell/metabolism/gtp_production",
                         "cell/genome/replicated_fraction", "cell/metabolism/dntp_pool"]}],
)
def mycoplasma_genitalium(core=None, *, nutrient_scale=1.0, initial_dnaA=0.0,
                          initial_dntp=0.0, seed=0):
    return build_mgen(core, nutrient_scale=nutrient_scale, initial_dnaA=initial_dnaA,
                      initial_dntp=initial_dntp, seed=seed)


def build_parca(core=None, *, seed=0, interval=1.0):
    """The parameter-calculator (ParCa) composite: runs the Karr FitConstants
    reproduction and emits its fitted-parameter summary + metabolic-feasibility
    closure. It is the study that PRECEDES the figure studies, which consume the
    per-gene parameters it computes (viva_mgen.parca)."""
    if core is None:
        from ..core import build_core
        core = build_core()
    from ..processes.parameter_calculator import ParameterCalculatorReproductionProcess as _PC
    proc = _PC(config={"seed": seed}, core=core)
    ports = list(proc.outputs())
    doc = {
        "parameter_calculator": {
            "_type": "process",
            "address": f"local:{_PC.__module__}.{_PC.__name__}",
            "config": {"seed": seed}, "interval": interval,
            "inputs": {}, "outputs": {p: ["parca", p] for p in ports},
        },
        "parca": {p: 0.0 for p in ports},
        "emitter": {"_type": "step", "address": "local:RAMEmitter",
                    "config": {"emit": {"parca": {p: "float" for p in ports}}},
                    "inputs": {"parca": ["parca"]}},
    }
    return doc


@composite_generator(
    name="mycoplasma_parca",
    description="Parameter calculator (ParCa) — the native reproduction of Karr 2012 FitConstants. Computes the per-gene expression/decay/synthesis panel from the real observed knowledge-base data, fits it under the RNA-mass / DnaA-FtsZ-held / net-supercoiling constraints, and verifies the expression↔metabolism feasibility closure. The upstream study whose fitted parameters the seven figure studies consume.",
    parameters={"seed": {"type": "integer", "default": 0, "description": "Stochastic seed"}},
    emitters=[{"address": "local:ParquetEmitter",
               "paths": ["parca/n_genes", "parca/median_mrna_synthesis_rate",
                         "parca/mrna_halflife_min_mean", "parca/closed_loop_feasible",
                         "parca/nmp_supply_flux", "parca/aa_supply_flux"]}],
)
def mycoplasma_parca(core=None, *, seed=0):
    return build_parca(core, seed=seed)
