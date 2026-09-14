"""Constant provenance audit: classify every reduced-process kinetic constant by
source tier (real_kb / real_supplement / order_of_magnitude / irreducible) so the
reproduction's parameterization fidelity is fully transparent. See
docs/superpowers/specs/2026-09-14-constant-provenance-design.md and gap #5.
"""
from __future__ import annotations

import math

from .processes import all_process_classes

SOURCE_TIERS = ("real_kb", "real_supplement", "order_of_magnitude", "irreducible")

# Explicit per-constant -> KB-key map: (class_name, constant) -> (kb_param_key,
# expected_value, note). A constant is classified real_kb ONLY if it appears here
# AND its config default equals expected_value (within tolerance) — this replaces
# the earlier value-coincidence match (any config default that happened to equal
# ANY numeric value in the process's whole KB dict), which produced false
# positives such as MetabolismFbaReproductionProcess.enzyme_coupling_cap=1.0
# coinciding with an unrelated KB 1.0, DNASupercoilingReproductionProcess.initial_sigma
# =0.0 matching any zero in the dict, ReplicationReproductionProcess.dnaA_threshold=30.0
# coinciding with ssbComplexSpacing, and DNARepairReproductionProcess.atp_per_repair=4.0
# coinciding with an unrelated NER incision-margin value. expected_value may be a
# documented transform of the raw KB value (see dna_pol_rate below); the transform
# is spelled out in the note. Entries are only added when the process source has a
# clear code-comment citation tying the constant to a specific real KB value.
_REAL_KB = {
    ("DNASupercoilingReproductionProcess", "gyrase_rate"): ("gyraseActivityRate", 1.2, "KB gyraseActivityRate"),
    ("DNASupercoilingReproductionProcess", "atp_per_act"): ("gyraseATPCost", 2.0, "KB gyraseATPCost"),
    ("ProteinProcessingIReproductionProcess", "deformylase_specific_rate"): ("deformylaseSpecificRate", 38.0, "KB deformylaseSpecificRate"),
    ("ProteinProcessingIReproductionProcess", "aminopeptidase_specific_rate"): ("methionineAminoPeptidaseSpecificRate", 6.0, "KB methionineAminoPeptidaseSpecificRate"),
    ("ProteinProcessingIIReproductionProcess", "peptidase_specific_rate"): ("lipoproteinSignalPeptidaseSpecificRate", 11.0, "KB lipoproteinSignalPeptidaseSpecificRate"),
    ("ProteinProcessingIIReproductionProcess", "transferase_specific_rate"): ("lipoproteinDiacylglycerylTransferaseSpecificRate", 0.0165, "KB lipoproteinDiacylglycerylTransferaseSpecificRate"),
    ("ProteinTranslocationReproductionProcess", "translocase_specific_rate"): ("translocaseSpecificRate", 2.71e12, "KB translocaseSpecificRate"),
    ("ProteinTranslocationReproductionProcess", "gtp_per_monomer"): ("SRP_GTPUsedPerMonomer", 2.0, "KB SRP_GTPUsedPerMonomer"),
    ("FtsZPolymerizationReproductionProcess", "activation_rate"): ("activationFwd", 1.1, "KB activationFwd"),
    ("FtsZPolymerizationReproductionProcess", "dissociation_rate"): ("activationRev", 0.01, "KB activationRev"),
    ("ReplicationReproductionProcess", "dna_pol_rate"): ("dnaPolymeraseElongationRate", 200.0, "2x KB dnaPolymeraseElongationRate (100)"),
}

# Curated overrides: {(class_name, constant): {"source_tier": ..., "note": ...}}.
#
# Provenance of this curation (Task 2, gap #5): the only supplement file present in
# this workspace is workspace/references/papers/Karr2012_supplementary_information.pdf
# (16 pages). Read in full: it is the Cell main-text article plus two SI figures
# (S1, S2) and the "Extended Experimental Procedures" page — NOT the ~1,900-parameter
# appendix. That page states explicitly: "Chapters 2 and 3 of Data S1 discuss the
# mathematics of each cellular variable and process submodel" — i.e. the actual
# per-constant kinetic-parameter tables live in a separate document (Data S1) that is
# not included in this repository and was not fetched (out of the bounded scope of
# this audit; see docs/CONSTANT_PROVENANCE.md). Consequently almost no order_of_magnitude
# constant could be matched to an explicit supplement value — the one exception found
# by exhaustive read is CELL_CYCLE_LENGTH_S (Figure 2A explicitly states "Mean tau =
# 9.0 h", matching the existing default exactly).
#
# For the remaining order_of_magnitude constants in the empty-KB-dict processes named
# in the Task 2 brief (ProteinFolding, ProteinModification, RNAModification,
# tRNAAminoacylation, RibosomeAssembly, DNADamage, ProteinActivation,
# TerminalOrganelleAssembly, HostInteraction — TranscriptionalRegulation and
# MacromolecularComplexation have no numeric config constants to audit), each constant
# was classified by asking: does it stand in for a per-species/per-enzyme rate matrix
# that literally does not reduce to a single scalar (the original evolveState loop
# draws separate rates per protein/enzyme/species from a table this reduced model
# doesn't carry, and the Karr KB parameter dict is confirmed empty via
# kb.karr_process_params) -> "irreducible"; or is it a well-established stoichiometric
# cost/count (N ATP or GTP per reaction event, a Hill exponent, a discretization
# resolution) that is a defensible textbook/engineering estimate rather than a missing
# measurement -> stays "order_of_magnitude" (curated here only to record the
# justification, not to change tier).
_CURATED: dict = {
    # --- real_supplement: found and cited in the Karr 2012 SI main text ---
    ("MassGrowthReproductionProcess", "cell_cycle_length_s"): {
        "source_tier": "real_supplement",
        "note": "Karr 2012 SI, Figure 2A: mean predicted/observed doubling time tau = 9.0 h "
                "(32,400 s) — matches the existing default exactly; see constants.py "
                "CELL_CYCLE_LENGTH_S.",
    },

    # --- ProteinFolding: KB ProteinFolding param dict confirmed empty; SI silent ---
    ("ProteinFoldingReproductionProcess", "spontaneous_rate"): {
        "source_tier": "irreducible",
        "note": "Karr KB ProteinFolding parameter dict is empty (kb.karr_process_params); the SI "
                "(16-page main text + Fig S1/S2, no Data S1 parameter appendix in this repo) states "
                "no per-protein folding rate. The original assigns each protein its own folding "
                "rate via the notFolding-set Gillespie loop — a per-protein table, not a scalar. "
                "Would need the full Data S1 folding-kinetics table (unavailable here).",
    },
    ("ProteinFoldingReproductionProcess", "chaperone_rate"): {
        "source_tier": "irreducible",
        "note": "Same as spontaneous_rate: KB ProteinFolding dict empty, SI silent, and the "
                "original's chaperone-substrate specificity is per-protein, not a single rate. "
                "Would need the full Data S1 chaperone-kinetics table (unavailable here).",
    },
    ("ProteinFoldingReproductionProcess", "atp_per_fold"): {
        "source_tier": "order_of_magnitude",
        "note": "~7 ATP per GroEL/GroES folding cycle is an established chaperonin-ATPase "
                "turnover estimate from the general chaperonin literature, not a Karr KB or SI "
                "value; kept as a defensible physiological estimate.",
    },

    # --- ProteinModification: KB dict confirmed empty; SI silent ---
    ("ProteinModificationReproductionProcess", "modification_specific_rate"): {
        "source_tier": "irreducible",
        "note": "Karr KB ProteinModification parameter dict is empty; SI silent. The original "
                "runs an enzyme+substrate Gillespie loop over specific phosphorylation/lipoylation "
                "reactions with per-enzyme specific rates — a matrix, not a scalar. Would need the "
                "full Data S1 per-enzyme modification-rate table (unavailable here).",
    },
    ("ProteinModificationReproductionProcess", "atp_per_modification"): {
        "source_tier": "order_of_magnitude",
        "note": "1 ATP per phosphoryl-transfer reaction is standard kinase stoichiometry, not a "
                "Karr KB or SI value; kept as a defensible physiological estimate.",
    },

    # --- RNAModification: KB dict confirmed empty; SI silent ---
    ("RNAModificationReproductionProcess", "modification_rate"): {
        "source_tier": "irreducible",
        "note": "Karr KB RNAModification parameter dict is empty; SI silent. The original covers "
                "13 modification enzymes over 86 base modifications with per-species rates — a "
                "matrix this reduced single-pool model cannot recover from a scalar. Would need "
                "the full Data S1 RNA-modification kinetics table (unavailable here).",
    },
    ("RNAModificationReproductionProcess", "enzyme_kcat"): {
        "source_tier": "irreducible",
        "note": "Same matrix gap as modification_rate: one lumped enzyme kcat stands in for "
                "13 distinct modification enzymes' kcats, which the Karr KB dict (empty) and this "
                "SI (silent) do not provide.",
    },

    # --- tRNAAminoacylation: KB dict confirmed empty; SI silent ---
    ("TRNAAminoacylationReproductionProcess", "synthetase_kcat"): {
        "source_tier": "irreducible",
        "note": "Karr KB tRNAAminoacylation parameter dict is empty; SI silent. One lumped "
                "synthetase kcat stands in for 20 amino acids x 37 tRNA/tmRNA reactions with "
                "distinct per-synthetase kcats (e.g. GlnRS vs glutamyl/methionyl transferases) — "
                "a matrix, not a scalar. Would need the full Data S1 synthetase-kinetics table "
                "(unavailable here).",
    },
    ("TRNAAminoacylationReproductionProcess", "atp_per_charge"): {
        "source_tier": "order_of_magnitude",
        "note": "1 ATP per aminoacylation (ATP -> AMP + PPi) is textbook aminoacyl-tRNA-synthetase "
                "stoichiometry, not a Karr KB or SI value; kept as a defensible physiological "
                "estimate.",
    },

    # --- RibosomeAssembly: KB dict confirmed empty; SI silent ---
    ("RibosomeAssemblyReproductionProcess", "gtp_per_complex"): {
        "source_tier": "order_of_magnitude",
        "note": "~2 GTP per subunit-assembly event (one GTPase assembly-factor turnover) is a "
                "defensible physiological estimate; not in the Karr KB (empty dict) or this SI.",
    },

    # --- DNADamage: KB dict confirmed empty; SI silent ---
    ("DNADamageReproductionProcess", "base_rate"): {
        "source_tier": "irreducible",
        "note": "Karr KB DNADamage parameter dict is empty; SI silent. The original places "
                "damage at vulnerable-motif sites with per-site/per-motif probabilities, not one "
                "lumped spontaneous rate. Would need the full Data S1 site-specific lesion-rate "
                "table (unavailable here).",
    },
    ("DNADamageReproductionProcess", "agent_rate"): {
        "source_tier": "irreducible",
        "note": "Same matrix gap as base_rate: a per-agent, per-site dose-response damage-rate "
                "table (radiation/reactive-species) is needed, which the Karr KB (empty) and this "
                "SI (silent) do not provide.",
    },
    ("DNADamageReproductionProcess", "n_bins"): {
        "source_tier": "order_of_magnitude",
        "note": "Chromosome-lesion-map spatial resolution (bin count), an engineering "
                "discretization choice rather than a measured kinetic constant; not applicable to "
                "KB/SI sourcing.",
    },

    # --- ProteinActivation: KB dict confirmed empty; SI silent ---
    ("ProteinActivationReproductionProcess", "default_k"): {
        "source_tier": "irreducible",
        "note": "Karr KB ProteinActivation parameter dict is empty; SI silent. The original's "
                "evaluateActivationRules is a per-protein Boolean/metabolite-threshold rule table "
                "(regulator_k is meant to carry it per-protein but defaults empty); this scalar "
                "fallback cannot recover that per-rule table. Would need the full Data S1 "
                "activation-rule table (unavailable here).",
    },
    ("ProteinActivationReproductionProcess", "hill_n"): {
        "source_tier": "order_of_magnitude",
        "note": "Hill coefficient of 2 is a conventional default for a cooperative-binding "
                "approximation of a Boolean rule; not a Karr KB or SI value.",
    },

    # --- TerminalOrganelleAssembly: KB dict confirmed empty; SI silent ---
    ("TerminalOrganelleAssemblyReproductionProcess", "threshold"): {
        "source_tier": "irreducible",
        "note": "Karr KB TerminalOrganelleAssembly parameter dict is empty; SI silent. The "
                "presence-at-threshold gate stands in for the full per-reaction localization "
                "stoichiometry matrix (protein copy-number thresholds per localization step), "
                "which neither the KB nor this SI provides.",
    },

    # --- HostInteraction: KB dict confirmed empty; SI silent; original rule is Boolean ---
    ("HostInteractionReproductionProcess", "cooperativity"): {
        "source_tier": "order_of_magnitude",
        "note": "The original HostInteraction rule is a Boolean AND (no continuous parameter "
                "exists to source); this Hill exponent is a free modeling choice softening that "
                "rule, not a value the KB or SI could supply.",
    },
    ("HostInteractionReproductionProcess", "half_max_fraction"): {
        "source_tier": "order_of_magnitude",
        "note": "Same as cooperativity: the original Boolean rule has no continuous half-max "
                "parameter to source from KB/SI; this is a free modeling choice.",
    },
}


def _is_number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def audit_constants() -> dict:
    classes = all_process_classes()
    out: dict = {}
    for cname, cls in classes.items():
        schema = getattr(cls, "config_schema", {}) or {}
        consts = {}
        for key, spec in schema.items():
            if not isinstance(spec, dict):
                continue
            default = spec.get("_default")
            if not _is_number(default):
                continue
            if key in ("seed",):
                continue
            default = float(default)
            entry = _REAL_KB.get((cname, key))
            match = entry is not None and abs(default - entry[1]) < 1e-9 * (abs(entry[1]) or 1)
            cur = _CURATED.get((cname, key))
            if match:
                tier, note = "real_kb", entry[2]
            elif cur:
                tier, note = cur["source_tier"], cur["note"]
            else:
                tier, note = "order_of_magnitude", "physiological estimate; not in KB parameters.json"
            consts[key] = {"value": default, "source_tier": tier,
                           "kb_match": bool(match), "note": note}
        if consts:
            out[cname] = consts
    return out
