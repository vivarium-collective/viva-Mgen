"""Constant provenance audit: classify every reduced-process kinetic constant by
source tier (real_kb / real_supplement / order_of_magnitude / irreducible) so the
reproduction's parameterization fidelity is fully transparent. See
docs/superpowers/specs/2026-09-14-constant-provenance-design.md and gap #5.
"""
from __future__ import annotations

import math

from .processes import all_process_classes
from . import kb

SOURCE_TIERS = ("real_kb", "real_supplement", "order_of_magnitude", "irreducible")

# Map process CLASS name -> its karr_parameters.json process key.
_KB_KEY = {
    "MetabolismFbaReproductionProcess": "Metabolism",
    "MassGrowthReproductionProcess": "",  # mass constants live in states, handled as order_of_magnitude/real via constants.py
    "TranscriptionReproductionProcess": "Transcription",
    "TranslationReproductionProcess": "Translation",
    "RnaDecayReproductionProcess": "RNADecay",
    "ProteinDecayReproductionProcess": "ProteinDecay",
    "ReplicationReproductionProcess": "Replication",
    "ReplicationInitiationReproductionProcess": "ReplicationInitiation",
    "DNASupercoilingReproductionProcess": "DNASupercoiling",
    "ChromosomeCondensationReproductionProcess": "ChromosomeCondensation",
    "ChromosomeSegregationReproductionProcess": "ChromosomeSegregation",
    "DNADamageReproductionProcess": "DNADamage",
    "DNARepairReproductionProcess": "DNARepair",
    "TranscriptionalRegulationReproductionProcess": "TranscriptionalRegulation",
    "RNAProcessingReproductionProcess": "RNAProcessing",
    "RNAModificationReproductionProcess": "RNAModification",
    "TRNAAminoacylationReproductionProcess": "tRNAAminoacylation",
    "ProteinProcessingIReproductionProcess": "ProteinProcessingI",
    "ProteinTranslocationReproductionProcess": "ProteinTranslocation",
    "ProteinProcessingIIReproductionProcess": "ProteinProcessingII",
    "ProteinFoldingReproductionProcess": "ProteinFolding",
    "ProteinModificationReproductionProcess": "ProteinModification",
    "ProteinActivationReproductionProcess": "ProteinActivation",
    "MacromolecularComplexationReproductionProcess": "MacromolecularComplexation",
    "RibosomeAssemblyReproductionProcess": "RibosomeAssembly",
    "TerminalOrganelleAssemblyReproductionProcess": "TerminalOrganelleAssembly",
    "FtsZPolymerizationReproductionProcess": "FtsZPolymerization",
    "CytokinesisReproductionProcess": "Cytokinesis",
    "HostInteractionReproductionProcess": "HostInteraction",
    "ChromosomeDynamicsReproductionProcess": "",
}

# Curated overrides: {(class_name, constant): {"source_tier": ..., "note": ...}}.
# Task 2 fills real_supplement / irreducible from the supplement read; a few known
# ones are seeded here. Only NON-real_kb constants should be curated.
_CURATED: dict = {}


def _is_number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def _kb_values(proc_key):
    if not proc_key:
        return set()
    d = kb.karr_process_params(proc_key) or {}
    out = set()
    for v in d.values():
        if _is_number(v):
            out.add(round(float(v), 12))
    return out


def audit_constants() -> dict:
    classes = all_process_classes()
    out: dict = {}
    for cname, cls in classes.items():
        schema = getattr(cls, "config_schema", {}) or {}
        kb_vals = _kb_values(_KB_KEY.get(cname, ""))
        consts = {}
        for key, spec in schema.items():
            if not isinstance(spec, dict):
                continue
            default = spec.get("_default")
            if not _is_number(default):
                continue
            if key in ("seed",):
                continue
            match = round(float(default), 12) in kb_vals
            cur = _CURATED.get((cname, key))
            if match:
                tier, note = "real_kb", f"value present in karr_parameters.json[{_KB_KEY.get(cname)}]"
            elif cur:
                tier, note = cur["source_tier"], cur["note"]
            else:
                tier, note = "order_of_magnitude", "physiological estimate; not in KB parameters.json"
            consts[key] = {"value": float(default), "source_tier": tier,
                           "kb_match": bool(match), "note": note}
        if consts:
            out[cname] = consts
    return out
