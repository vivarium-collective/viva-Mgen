# Constant Provenance Implementation Plan (gap #5)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans. Steps use checkbox (`- [ ]`).

**Goal:** Give every reduced-process kinetic constant an explicit, auditable provenance (real_kb / real_supplement / order_of_magnitude / irreducible), source what the accessible Karr supplement provides, and document the rest as irreducible — closing gap #5.

**Architecture:** `viva_mgen/provenance.py` audits each process's `config_schema` numeric constants and cross-references `karr_parameters.json` (→ real_kb) and a curated overrides table (→ real_supplement / irreducible / justified order_of_magnitude). `scripts/extract_constant_provenance.py` writes `datasets/constant_provenance.json`; `docs/CONSTANT_PROVENANCE.md` summarizes by tier.

**Tech Stack:** Python 3.12, process-bigraph, pytest. Run with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--constants` and venv `/Users/eranagmon/code/viva-mGen/.venv`.

**Spec:** `docs/superpowers/specs/2026-09-14-constant-provenance-design.md`

## Global Constraints

- Work only in worktree `/Users/eranagmon/code/viva-mGen--constants` (branch `feat/constant-provenance`).
- Run everything with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--constants /Users/eranagmon/code/viva-mGen/.venv/bin/python`.
- No AI-attribution trailers.
- Do NOT fabricate literature kinetic values. A value is only `real_supplement` if actually found in `workspace/references/papers/Karr2012_supplementary_information.pdf` (cite where). Otherwise `order_of_magnitude` (physiological estimate) or `irreducible` (KB dict empty AND supplement silent — note what external source is needed).
- The only behavior change permitted is wiring a genuinely-sourced supplement value into a process default (with a citing comment); everything else is documentation.

## Facts: 11 empty KB param dicts — DNADamage, HostInteraction, MacromolecularComplexation, ProteinActivation, ProteinFolding, ProteinModification, RNAModification, RibosomeAssembly, TerminalOrganelleAssembly, TranscriptionalRegulation, tRNAAminoacylation. `kb.karr_process_params(name)` returns a process's KB dict; `all_process_classes()` enumerates process classes.

---

### Task 1: provenance audit module + dataset + tests

**Files:**
- Create: `viva_mgen/provenance.py`
- Create: `scripts/extract_constant_provenance.py`
- Create: `datasets/constant_provenance.json` (generated)
- Test: `tests/test_provenance.py`

**Interfaces:**
- Produces: `audit_constants() -> {process: {constant: {value, source_tier, kb_match, note}}}`; module constant `SOURCE_TIERS`; curated overrides `_CURATED` (initially with the confirmed real_kb/order_of_magnitude cases; Task 2 fills real_supplement/irreducible).

- [ ] **Step 1: Write failing tests**

```python
# tests/test_provenance.py
from viva_mgen.provenance import audit_constants, SOURCE_TIERS

def test_audit_covers_processes_and_shape():
    a = audit_constants()
    assert len(a) >= 10
    for proc, consts in a.items():
        for name, rec in consts.items():
            assert set(rec) >= {"value", "source_tier", "kb_match", "note"}
            assert rec["source_tier"] in SOURCE_TIERS

def test_real_kb_classification():
    a = audit_constants()
    # DNASupercoiling gyrase_rate 1.2 is in karr_parameters.json -> real_kb
    assert a["DNASupercoilingReproductionProcess"]["gyrase_rate"]["source_tier"] == "real_kb"
    # ProteinFolding spontaneous_rate has no KB entry -> not real_kb
    assert a["ProteinFoldingReproductionProcess"]["spontaneous_rate"]["source_tier"] != "real_kb"

def test_no_contradiction():
    a = audit_constants()
    for consts in a.values():
        for rec in consts.values():
            if rec["source_tier"] == "real_kb":
                assert rec["kb_match"] is True
            if rec["source_tier"] == "irreducible":
                assert rec["kb_match"] is False
```

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--constants /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_provenance.py -q`
Expected: FAIL — module missing.

- [ ] **Step 3: Implement `viva_mgen/provenance.py`**

```python
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
```

- [ ] **Step 4: Run to verify pass**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--constants /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_provenance.py -q`
Expected: PASS (3 tests). If `test_real_kb_classification` fails on `gyrase_rate`, inspect `karr_parameters.json[DNASupercoiling]` — confirm 1.2 is `gyraseActivityRate` (it is) and that the value-set match works.

- [ ] **Step 5: Extraction script + generate dataset**

Create `scripts/extract_constant_provenance.py`:
```python
#!/usr/bin/env python
"""Write datasets/constant_provenance.json from viva_mgen.provenance.audit_constants()."""
from __future__ import annotations
import json
from collections import Counter
from pathlib import Path
from viva_mgen.provenance import audit_constants, SOURCE_TIERS

OUT = Path(__file__).resolve().parents[1] / "datasets" / "constant_provenance.json"

def main():
    a = audit_constants()
    tier = Counter()
    n = 0
    for consts in a.values():
        for rec in consts.values():
            tier[rec["source_tier"]] += 1
            n += 1
    payload = {"source": "viva_mgen.provenance.audit_constants() — kinetic-constant provenance vs Karr 2012 KB + supplement",
               "n_constants": n, "by_tier": {t: tier.get(t, 0) for t in SOURCE_TIERS},
               "constants": a}
    OUT.write_text(json.dumps(payload, indent=2))
    print(f"wrote {OUT}: {n} constants, by_tier={dict(tier)}")

if __name__ == "__main__":
    main()
```
Run it; confirm it reports a plausible breakdown (a chunk `real_kb`, several `order_of_magnitude`).

- [ ] **Step 6: Commit**

```bash
git add viva_mgen/provenance.py scripts/extract_constant_provenance.py datasets/constant_provenance.json tests/test_provenance.py
git commit -m "Add kinetic-constant provenance audit (real_kb vs order_of_magnitude vs irreducible)"
```

---

### Task 2: supplement sourcing + curated overrides + doc + resources

**Files:**
- Modify: `viva_mgen/provenance.py` (`_CURATED`)
- Possibly modify: process `config_schema` defaults (only if a real supplement value is found — with a citing comment)
- Create: `docs/CONSTANT_PROVENANCE.md`
- Modify: `viva_mgen/data_sources.py` (Resources row); `docs/FIDELITY_GAPS.md`
- Regenerate: `datasets/constant_provenance.json`
- Test: `tests/test_provenance.py` (append)

**Interfaces:**
- Consumes: Task 1 audit.

- [ ] **Step 1: Read the Karr supplement for the empty-dict processes**

Read `workspace/references/papers/Karr2012_supplementary_information.pdf` (use the Read tool with page ranges) focusing on kinetic-parameter tables for the processes with empty KB dicts and order_of_magnitude constants: ProteinFolding, ProteinModification, RNAModification, tRNAAminoacylation, RibosomeAssembly, TranscriptionalRegulation, MacromolecularComplexation, DNADamage, ProteinActivation, TerminalOrganelleAssembly, HostInteraction. For EACH order_of_magnitude constant in the Task-1 audit, determine: does the supplement state a value?
  - If YES: record the value + the supplement location; set `_CURATED[(class, const)] = {"source_tier":"real_supplement","note":"Karr 2012 SI, <section/table>: <value>"}`, and WIRE the value into that process's `config_schema` default (edit the process file, add a comment `# Karr 2012 SI <ref>`). Only do this when the supplement value is unambiguous.
  - If NO: set `_CURATED[(class, const)] = {"source_tier":"irreducible" or "order_of_magnitude", "note": "<why: KB dict empty; SI gives no value; would need <external source>>"}`. Use `irreducible` when the KB dict is empty AND the SI is silent AND it's a genuine per-species/matrix constant the reduced model can't recover; use `order_of_magnitude` when a defensible physiological estimate stands (e.g. ATP-per-reaction ~1).

Be conservative and honest: do NOT invent values. If the SI is ambiguous, leave it order_of_magnitude with a note.

- [ ] **Step 2: Regenerate dataset + write the doc**

Re-run `scripts/extract_constant_provenance.py`. Then write `docs/CONSTANT_PROVENANCE.md`: an intro paragraph + a per-tier table (tier, count, and the `process.constant` list), generated from `datasets/constant_provenance.json` (you may compute it inline in the doc text). This is the transparency artifact — a reader sees how much of the parameterization is real vs estimated.

- [ ] **Step 3: Resources row + fidelity notes + backlog**

- `data_sources.py`: add a `parameters`-category row for `datasets/constant_provenance.json` (match the existing `_FILES` tuple format).
- In each process whose constants are `order_of_magnitude`/`irreducible`, append a short clause to its `description` fidelity line: "constant provenance recorded in the constant-provenance audit (gap #5)." (One clause per process; don't restate per-constant. Keep first lines unchanged.)
- `docs/FIDELITY_GAPS.md`: mark gap #5 done (provenance audit + supplement sourcing; irreducible constants documented).

- [ ] **Step 4: Tests (append)**

```python
def test_provenance_dataset_consistent():
    import json
    from pathlib import Path
    d = json.loads(Path("datasets/constant_provenance.json").read_text())
    assert sum(d["by_tier"].values()) == d["n_constants"]
    # at least some constants are real_kb and the audit ran over many processes
    assert d["by_tier"]["real_kb"] > 0 and d["n_constants"] > 20

def test_curated_overrides_valid():
    from viva_mgen.provenance import _CURATED, SOURCE_TIERS, audit_constants
    a = audit_constants()
    for (cls, const), rec in _CURATED.items():
        assert rec["source_tier"] in SOURCE_TIERS and rec["note"]
        # curated entries must reference real audited constants
        assert cls in a and const in a[cls]
```

- [ ] **Step 5: Run full suite + commit**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--constants /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/ -q` (all pass). If any process default was changed to a supplement value, note it (behavior change) in the report. Regenerate the composite-state snapshot ONLY if a config default changed (`scripts/regen_composite_state.py`).
```bash
git add -A
git commit -m "Constant provenance: source from Karr SI, document irreducible constants, add audit doc"
```

---

## Self-Review

**Spec coverage:** audit module + dataset (T1); supplement sourcing + curated tiers (T2 S1); doc (T2 S2); resources + fidelity notes + backlog (T2 S3); tests incl. consistency + curated validity (T1+T2). Covered.

**Placeholder scan:** T1 is concrete code; T2's supplement-read is inherently a judgment step but bounded by explicit rules (found→wire+cite; not-found→irreducible/order_of_magnitude with note; never invent).

**Type consistency:** `audit_constants() -> {class: {const: {value,source_tier,kb_match,note}}}`; `SOURCE_TIERS` tuple; `_CURATED` keyed by (class,const). Consistent across module, script, tests.
