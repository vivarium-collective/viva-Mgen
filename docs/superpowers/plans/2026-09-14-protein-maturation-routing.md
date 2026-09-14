# Per-protein Maturation Routing Implementation Plan (gap #4 phase 1)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans. Steps use checkbox (`- [ ]`).

**Goal:** Route the protein-maturation steps to the real per-protein subsets decoded from the KB — lipoprotein-only Lgt transferase (ProcessingII), N-terminal Met cleavage subset (ProcessingI), secretory/lipoprotein-only translocation — replacing "applies to all".

**Architecture:** `kb_decode.decode_protein_maturation` reads each ProteinMonomer's `signalSequenceType`/`nTerminalMethionineCleavage`/`signalSequenceLength` into `datasets/karr_protein_maturation.json`; loaders in `kb.py` key it by the expression-panel key; the three protein processes gain a `*_genes` config (defaulting to the decoded subset) that gates the relevant step.

**Tech Stack:** Python 3.12, process-bigraph, pytest. Run with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--submodels` and venv `/Users/eranagmon/code/viva-mGen/.venv`.

**Spec:** `docs/superpowers/specs/2026-09-14-protein-maturation-routing-design.md`

## Global Constraints

- Work only in worktree `/Users/eranagmon/code/viva-mGen--submodels` (branch `feat/submodel-dereduction`).
- Run everything with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--submodels /Users/eranagmon/code/viva-mGen/.venv/bin/python`.
- No AI-attribution trailers.
- The maturation-form stores are a parallel representation; routing must PASS EVERY protein through the pipeline (never drop one), only GATING which step's enzyme limit / metabolite cost applies — so deformylation still advances all, translocation passes non-substrates through without consuming GTP, etc. Verify the composite still builds + runs.
- Real KB rates already in place (deformylase 38, aminopeptidase 6, peptidase 11, transferase 0.0165, SRP 2 GTP/monomer) must be unchanged.

## KB facts (verified): signalSequenceType = {None:448, secretory:20, lipoprotein:14}; nTerminalMethionineCleavage = {0:447, 1:35}; MG_067_MONOMER is a lipoprotein. ProteinMonomer ids are `MG_###_MONOMER` → gene `MG_###`.

---

### Task 1: decode maturation classification + dataset + loaders

**Files:**
- Modify: `viva_mgen/kb_decode.py` (add `decode_protein_maturation`)
- Create: `scripts/extract_kb_maturation.py`
- Create: `datasets/karr_protein_maturation.json` (generated)
- Modify: `viva_mgen/kb.py` (add `load_protein_maturation`, `maturation_by_panel_key`)
- Modify: `viva_mgen/data_sources.py` (Resources row)
- Test: `tests/test_protein_maturation.py`

**Interfaces:**
- Produces: `decode_protein_maturation(mat) -> {gene_id: {"signal_type": str, "met_cleavage": bool, "signal_length": float}}`; `kb.load_protein_maturation() -> dict`; `kb.maturation_by_panel_key() -> dict`.

- [ ] **Step 1: Add `decode_protein_maturation` to kb_decode.py**

```python
def decode_protein_maturation(mat_path):
    """Per-protein maturation classification from the KB ProteinMonomer objects:
    {gene_id: {"signal_type": None|"secretory"|"lipoprotein",
               "met_cleavage": bool, "signal_length": float}}.
    Gene id is the monomer's MG_###_MONOMER stripped of the _MONOMER suffix."""
    cells = _read_subsystem_cells(mat_path)
    by = _objects_by_class(cells)
    out = {}
    for p in by.get("ProteinMonomer", []):
        wid = _scalar(p.get("wholeCellModelID"))
        if wid is None:
            continue
        gid = str(wid)
        if gid.endswith("_MONOMER"):
            gid = gid[: -len("_MONOMER")]
        sst = _scalar(p.get("signalSequenceType"))
        sst = str(sst) if sst is not None and str(sst) != "None" else None
        met = _scalar(p.get("nTerminalMethionineCleavage"))
        try:
            slen = float(_scalar(p.get("signalSequenceLength")) or 0.0)
        except (TypeError, ValueError):
            slen = 0.0
        out[gid] = {"signal_type": sst,
                    "met_cleavage": bool(met and float(met) > 0),
                    "signal_length": slen}
    return out
```

- [ ] **Step 2: Extraction script `scripts/extract_kb_maturation.py`**

```python
#!/usr/bin/env python
"""Extract per-protein maturation classification from knowledgeBase.mat into
datasets/karr_protein_maturation.json. Usage: python scripts/extract_kb_maturation.py [mat]"""
from __future__ import annotations
import json, sys
from pathlib import Path
from viva_mgen import kb_decode as k

DEFAULT_MAT = Path.home() / "code/WholeCell-reference/data/knowledgeBase.mat"
OUT = Path(__file__).resolve().parents[1] / "datasets" / "karr_protein_maturation.json"

def main():
    mat = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_MAT
    m = k.decode_protein_maturation(mat)
    from collections import Counter
    st = Counter(v["signal_type"] for v in m.values())
    payload = {"source": "Karr 2012 knowledgeBase.mat (ProteinMonomer signal/Met classification)",
               "n_proteins": len(m),
               "n_lipoprotein": sum(1 for v in m.values() if v["signal_type"] == "lipoprotein"),
               "n_secretory": sum(1 for v in m.values() if v["signal_type"] == "secretory"),
               "n_met_cleavage": sum(1 for v in m.values() if v["met_cleavage"]),
               "maturation": {g: v for g, v in sorted(m.items())}}
    OUT.write_text(json.dumps(payload, indent=2))
    print(f"wrote {OUT}: {len(m)} proteins, signal_type={dict(st)}")

if __name__ == "__main__":
    main()
```
Run it: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--submodels /Users/eranagmon/code/viva-mGen/.venv/bin/python scripts/extract_kb_maturation.py` and confirm it reports 14 lipoprotein / 20 secretory / 35 met_cleavage.

- [ ] **Step 3: Loaders in kb.py**

```python
@functools.lru_cache(maxsize=1)
def load_protein_maturation() -> dict:
    """gene_id -> {signal_type, met_cleavage, signal_length} from
    datasets/karr_protein_maturation.json (see scripts/extract_kb_maturation.py)."""
    path = dataset_path("karr_protein_maturation.json")
    if not path.is_file():
        return {}
    with open(path) as f:
        return json.load(f).get("maturation", {})


@functools.lru_cache(maxsize=1)
def maturation_by_panel_key() -> dict:
    """Same, re-keyed by the expression-panel key (symbol-or-id) so the protein
    maturation processes (which use panel-keyed maps) can classify each protein."""
    mat = load_protein_maturation()
    out = {}
    for g in load_genes():
        gid = g["gene_id"]
        if gid in mat:
            key = (g.get("symbol") or "").strip() or gid
            out[key] = mat[gid]
    return out
```

- [ ] **Step 4: Write failing tests, then run**

```python
# tests/test_protein_maturation.py
from viva_mgen.kb import load_protein_maturation, maturation_by_panel_key

def test_counts():
    m = load_protein_maturation()
    assert sum(1 for v in m.values() if v["signal_type"] == "lipoprotein") == 14
    assert sum(1 for v in m.values() if v["signal_type"] == "secretory") == 20
    assert sum(1 for v in m.values() if v["met_cleavage"]) == 35
    assert m["MG_067"]["signal_type"] == "lipoprotein"

def test_panel_key_map():
    pk = maturation_by_panel_key()
    assert len(pk) > 100 and any(v["signal_type"] == "lipoprotein" for v in pk.values())
```

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--submodels /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_protein_maturation.py -q` → PASS after the dataset is generated.

- [ ] **Step 5: Resources row (data_sources.py)**

Add to the `_FILES` list (knowledge-base category): `("protein-maturation", "datasets/karr_protein_maturation.json", "knowledge-base", "Per-protein maturation classification decoded from the Karr 2012 KB: signal-sequence type (14 lipoproteins, 20 secretory), N-terminal Met-cleavage flag (35 proteins), signal length. Routes the protein-maturation submodels to the correct protein subsets.")`.

- [ ] **Step 6: Commit**

```bash
git add viva_mgen/kb_decode.py scripts/extract_kb_maturation.py datasets/karr_protein_maturation.json viva_mgen/kb.py viva_mgen/data_sources.py tests/test_protein_maturation.py
git commit -m "Decode per-protein maturation classification (signal type, Met cleavage) from the KB"
```

---

### Task 2: route ProcessingI / ProcessingII / Translocation to the real subsets

**Files:**
- Modify: `viva_mgen/processes/protein.py` (ProcessingI, ProcessingII, Translocation)
- Test: `tests/test_protein_maturation.py` (append)
- Regenerate: `reports/composite-state/*.json`

**Interfaces:**
- Consumes: `kb.maturation_by_panel_key` (Task 1).

Read the current `ProteinProcessingIReproductionProcess`, `ProteinProcessingIIReproductionProcess`, `ProteinTranslocationReproductionProcess` in `viva_mgen/processes/protein.py` first — mirror their existing idioms.

- [ ] **Step 1: ProcessingII — lipoprotein-only Lgt transferase**

Add config `lipoprotein_genes` (`{"_type":"list[string]","_default":[]}`); in `__init__`, if empty, default to `[k for k,v in maturation_by_panel_key().items() if v["signal_type"]=="lipoprotein"]` (import the loader). Add config `transferase_specific_rate` (`0.0165`). In `update()`: the signal peptidase (rate 11) applies to ALL translocated as today; ADDITIONALLY, for proteins in the lipoprotein set, throughput is further limited by the transferase (`transferase_count * 0.0165 * interval` shared across the lipoprotein subset). Simplest correct form: compute the peptidase-limited move for all; then for the lipoprotein subset, cap their move by the transferase limit. Emit as before. Add a `diacylglyceryl_transferase` input (float, default 40.0) for the transferase enzyme count. Update the description fidelity line (lipoprotein subset now gets the Lgt step, from real KB).

- [ ] **Step 2: ProcessingI — Met-cleavage subset**

Add config `met_cleavage_genes` (default = `[k for k,v in maturation_by_panel_key().items() if v["met_cleavage"]]`). Deformylase (rate 38) applies to all nascent as today; the aminopeptidase (rate 6) enzyme limit now gates ONLY the met-cleavage subset (proteins not in the set are deformylated and pass to process_i_done without consuming aminopeptidase capacity). Keep both subsets moving nascent→process_i_done. Update the fidelity line.

- [ ] **Step 3: Translocation — secretory+lipoprotein only**

Add config `translocated_genes` (default = `[k for k,v in maturation_by_panel_key().items() if v["signal_type"] in ("secretory","lipoprotein")]`). Only these consume translocase/GTP and are "translocated"; a non-substrate protein in `process_i_done` passes through to `translocated` WITHOUT consuming GTP (so the pipeline doesn't stall). Update the fidelity line.

- [ ] **Step 4: Tests (append)**

```python
def test_processing_ii_lipoprotein_throttled():
    from viva_mgen.processes.protein import ProteinProcessingIIReproductionProcess as P2
    from viva_mgen.core import build_core
    p = P2(config={"lipoprotein_genes": ["lipoA"]}, core=build_core())
    # a lipoprotein and a non-lipoprotein both present in large amounts
    out = p.update({"translocated": {"lipoA": 1000.0, "plainB": 1000.0},
                    "signal_peptidase": 100.0, "diacylglyceryl_transferase": 1.0}, 1.0)
    moved = out["processed_ii"]
    # lipoA limited by transferase (1 * 0.0165) -> ~0; plainB limited only by peptidase -> more
    assert moved.get("plainB", 0) > moved.get("lipoA", 0)

def test_translocation_only_substrates():
    from viva_mgen.processes.protein import ProteinTranslocationReproductionProcess as T
    from viva_mgen.core import build_core
    p = T(config={"translocated_genes": ["secA"]}, core=build_core())
    out = p.update({"process_i_done": {"secA": 100.0, "cytoC": 100.0},
                    "translocase": 1e12, "gtp": 1e9}, 1.0)
    # both reach 'translocated' (pipeline doesn't stall) but only secA consumed GTP
    assert out["gtp"] < 0     # GTP consumed for secA
    # cytoC passed through without GTP: its translocated delta present, but GTP cost only for secA
```
(Adjust assertions to the exact implementation; the intent: lipoprotein throttled vs not; only substrates consume GTP.)

- [ ] **Step 5: Composite still builds; add the new `diacylglyceryl_transferase` port to the composite store map**

In `viva_mgen/composites/mgen.py` `_STORE_GROUP`, add `"diacylglyceryl_transferase": "proteome"` and a `_FLOAT_INIT` default (e.g. 40.0). Build + run the composite 3 s to confirm no wiring error.

- [ ] **Step 6: Run full suite + regenerate snapshot + commit**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--submodels /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/ -q` (all pass).
Then: `... scripts/regen_composite_state.py`.
```bash
git add -A
git commit -m "Protein maturation: route Lgt/Met-cleavage/translocation to real KB protein subsets"
```

---

## Self-Review

**Spec coverage:** decoder+dataset+loaders (T1); ProcessingII lipoprotein (T2 S1); ProcessingI Met subset (T2 S2); Translocation targeting (T2 S3); Resources row (T1 S5); tests (both); snapshot (T2 S6). Covered.

**Placeholder scan:** decoder/script/loaders/tests are concrete; the T2 process edits are described against the real existing code (implementer reads protein.py) with concrete defaults + test intent — acceptable since exact line edits depend on current code, but the classification sets + rates + port names are all specified.

**Type consistency:** `decode_protein_maturation -> {gene_id: {signal_type,met_cleavage,signal_length}}`; loaders return same; process configs `lipoprotein_genes`/`met_cleavage_genes`/`translocated_genes` (list[string]); new input `diacylglyceryl_transferase` (float). Consistent.
