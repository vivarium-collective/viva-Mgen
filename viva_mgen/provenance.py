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
# Task 2 fills real_supplement / irreducible from the supplement read; a few known
# ones are seeded here. Only NON-real_kb constants should be curated.
_CURATED: dict = {}


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
