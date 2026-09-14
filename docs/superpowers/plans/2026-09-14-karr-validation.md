# Karr-comparison Gating Metrics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans. Steps use checkbox (`- [ ]`).

**Goal:** Add Karr-cited **gating metrics to the study report cards** — new `behavior_tests` on the *emergent* observables (emergent mass composition from gap #2, single-cell mRNA CV from gap #7), backed by a versioned Karr reference dataset + shared extraction helpers — so the comparison to the paper is evaluated where studies gate, and moves as remaining gaps land.

**Architecture:** `datasets/karr_reference_values.json` holds the paper's targets/bands. `viva_mgen/validation.py` provides pure extraction/scoring helpers. `fig2-growth/sims/run.py` computes the new derived scalars (calling the helpers + a small ensemble) and records them in the completion event's `observables`; `fig2-growth/study.yaml` adds matching `behavior_tests` + `expected_behavior` with `pass_if` bands citing Karr.

**Tech Stack:** Python 3.12, numpy, process-bigraph, pytest. Run with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--validation` and venv `/Users/eranagmon/code/viva-mGen/.venv`.

**Spec:** `docs/superpowers/specs/2026-09-14-karr-validation-design.md`

## Global Constraints

- Work only in worktree `/Users/eranagmon/code/viva-mGen--validation` (branch `feat/karr-validation`).
- Run everything with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--validation /Users/eranagmon/code/viva-mGen/.venv/bin/python`.
- No AI-attribution trailers.
- New emergent gates are `classification: secondary` (may fail at baseline — that's the honest improvement signal; don't fudge bands to force a pass). Don't touch the existing tautological `protein_fraction`/`rna_fraction` tests.
- Keep fig2's added ensemble modest (n≈4, ~1h) so the study run stays tractable.
- Report-card `behavior_tests` are evaluated in the dashboard, not pytest CI — so a failing new gate does NOT break CI.

## How the study system works (verified)
`run.py main()` computes derived scalars and records them in `append_run_event(..., "observables": {field: value})` on completion. A `behavior_test` with `measure: {kind: derived_scalar, field: X}` + `pass_if: {op: range, low, high, provenance: {kind: experiment, note}}` gates against `observables[X]`. `expected_behavior` mirrors it (name, en, observable, condition, rationale, measure, expect).

---

### Task 1: Karr reference dataset + validation helpers

**Files:**
- Create: `datasets/karr_reference_values.json`
- Create: `viva_mgen/validation.py`
- Test: `tests/test_validation.py`

**Interfaces:**
- Produces: `reference(id) -> dict`; `emergent_macro_fractions(row) -> dict`; `mrna_cv(final_totals) -> float`; `score(observed, ref) -> dict`.

- [ ] **Step 1: Write `datasets/karr_reference_values.json`**

```json
{
  "source": "Karr et al. 2012, Cell 150:389-401 — reported quantitative results",
  "values": {
    "mass_fraction_protein": {"description": "protein fraction of macromolecular dry mass", "target": 0.705, "band": [0.55, 0.80], "unit": "fraction", "source": "Karr 2012 cell composition (protein 0.62 of total; 0.705 renormalized over protein+DNA+RNA)", "observable": "emergent_protein_fraction"},
    "mass_fraction_dna": {"description": "DNA fraction of macromolecular dry mass", "target": 0.192, "band": [0.10, 0.35], "unit": "fraction", "source": "Karr 2012 (DNA 0.169 of total; 0.192 renormalized)", "observable": "emergent_dna_fraction"},
    "mass_fraction_rna": {"description": "RNA fraction of macromolecular dry mass", "target": 0.106, "band": [0.05, 0.25], "unit": "fraction", "source": "Karr 2012 (RNA 0.093 of total; 0.106 renormalized)", "observable": "emergent_rna_fraction"},
    "single_cell_mrna_cv": {"description": "coefficient of variation of per-cell total mRNA across independent cells", "target": 0.3, "band": [0.001, 1.5], "unit": "cv", "source": "Karr 2012 Fig 2 single-cell distributions (bursty low-copy mRNA)", "observable": "single_cell_mrna_cv"},
    "doubling_time_h": {"description": "cell-cycle length", "target": 9.0, "band": [7.0, 11.0], "unit": "h", "source": "Karr 2012 Fig 2A (mean tau = 9.0 h)", "observable": "doubling_time_h"}
  }
}
```

- [ ] **Step 2: Write failing tests**

```python
# tests/test_validation.py
import math
from viva_mgen.validation import reference, emergent_macro_fractions, mrna_cv, score

def test_reference_entries_valid():
    for vid in ("mass_fraction_protein", "single_cell_mrna_cv", "doubling_time_h"):
        r = reference(vid)
        assert set(r) >= {"description", "target", "band", "unit", "source", "observable"}
        assert r["band"][0] <= r["band"][1]

def test_emergent_macro_fractions_renormalizes():
    row = {"emergent_mass_fractions": {"protein": 0.6, "DNA": 0.2, "RNA": 0.1, "metabolite": 0.1}}
    f = emergent_macro_fractions(row)  # renormalize over protein+DNA+RNA (drop metabolite)
    assert abs(sum(f.values()) - 1.0) < 1e-9
    assert abs(f["protein"] - 0.6/0.9) < 1e-9 and "metabolite" not in f

def test_emergent_macro_fractions_empty():
    assert emergent_macro_fractions({"emergent_mass_fractions": {}}) == {"protein": 0.0, "DNA": 0.0, "RNA": 0.0}

def test_mrna_cv():
    assert mrna_cv([100.0, 100.0, 100.0]) == 0.0
    cv = mrna_cv([80.0, 100.0, 120.0])
    assert abs(cv - (20.0*math.sqrt(2/3)) / 100.0) < 1e-9   # population std / mean

def test_score():
    r = reference("doubling_time_h")  # target 9, band [7,11]
    assert score(9.0, r)["pass"] is True and score(9.0, r)["closeness"] == 1.0
    assert score(20.0, r)["pass"] is False and score(20.0, r)["closeness"] < 0.5
```

- [ ] **Step 3: Implement `viva_mgen/validation.py`**

```python
"""Karr-2012 comparison helpers: load reference targets and extract/score the
emergent observables the study report cards gate on. See
docs/superpowers/specs/2026-09-14-karr-validation-design.md.
"""
from __future__ import annotations

import functools
import json

from .kb import dataset_path

_MACRO = ("protein", "DNA", "RNA")


@functools.lru_cache(maxsize=1)
def _reference_values() -> dict:
    path = dataset_path("karr_reference_values.json")
    with open(path) as f:
        return json.load(f).get("values", {})


def reference(vid: str) -> dict:
    vals = _reference_values()
    if vid not in vals:
        raise KeyError(f"no Karr reference value {vid!r}")
    return vals[vid]


def emergent_macro_fractions(row) -> dict:
    """Renormalize a row's emergent_mass_fractions over the macromolecules
    {protein, DNA, RNA} (drop metabolite/other); all-zero → zeros."""
    emf = (row or {}).get("emergent_mass_fractions", {}) or {}
    vals = {k: float(emf.get(k, 0.0)) for k in _MACRO}
    tot = sum(vals.values())
    if tot <= 0:
        return {k: 0.0 for k in _MACRO}
    return {k: vals[k] / tot for k in _MACRO}


def mrna_cv(final_totals) -> float:
    """Population coefficient of variation of per-cell final total mRNA."""
    xs = [float(x) for x in (final_totals or []) if x is not None]
    if len(xs) < 2:
        return 0.0
    mean = sum(xs) / len(xs)
    if mean == 0:
        return 0.0
    var = sum((x - mean) ** 2 for x in xs) / len(xs)
    return (var ** 0.5) / mean


def score(observed, ref) -> dict:
    lo, hi = ref["band"]
    target = float(ref["target"])
    passed = (observed is not None) and (lo <= float(observed) <= hi)
    if observed is None:
        closeness = 0.0
    else:
        denom = max(abs(target), 1e-9)
        closeness = max(0.0, 1.0 - abs(float(observed) - target) / denom)
    return {"observed": (None if observed is None else float(observed)),
            "target": target, "band": [lo, hi], "pass": bool(passed),
            "closeness": float(closeness)}
```

- [ ] **Step 4: Run tests to pass**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--validation /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_validation.py -q`
Expected: PASS (5 tests). (`dataset_path` is an existing helper in `viva_mgen/kb.py`.)

- [ ] **Step 5: Resources row (data_sources.py) + commit**

Add a `reference`-category row for `datasets/karr_reference_values.json` (match the existing `_FILES` tuple format; note: "Curated Karr 2012 reported quantitative results — targets + tolerance bands the study report cards gate against.").
```bash
git add datasets/karr_reference_values.json viva_mgen/validation.py viva_mgen/data_sources.py tests/test_validation.py
git commit -m "Add Karr reference values + comparison-metric helpers for study gating"
```

---

### Task 2: wire emergent + variation gating metrics into fig2-growth

**Files:**
- Modify: `workspace/studies/fig2-growth/sims/run.py`
- Modify: `workspace/studies/fig2-growth/study.yaml`
- Modify: `docs/FIDELITY_GAPS.md` (note the comparison harness)
- Test: `tests/test_validation.py` (append)

**Interfaces:**
- Consumes: `viva_mgen.validation`, `viva_mgen.ensemble` (gap #7).

Read the current `workspace/studies/fig2-growth/sims/run.py` and `study.yaml` first.

- [ ] **Step 1: Compute the new derived scalars in run.py `main()`**

After the existing derived-scalar block, add (using the final emitter row + a small ensemble):
```python
from viva_mgen import validation as val
from viva_mgen.ensemble import run_ensemble

emf = val.emergent_macro_fractions(rows[-1])
emergent_protein_fraction = emf["protein"]
emergent_dna_fraction = emf["DNA"]
emergent_rna_fraction = emf["RNA"]

# single-cell mRNA variation across an independent-cell ensemble (Fig 2 distributions)
ens = run_ensemble(n_cells=4, duration=3600.0)  # 1 h is enough for mRNA CV; keep it modest
def _tot_mrna(cell_rows):
    rc = cell_rows[-1].get("rna_counts", {}) if cell_rows else {}
    return sum(rc.values()) if isinstance(rc, dict) else 0.0
single_cell_mrna_cv = val.mrna_cv([_tot_mrna(c) for c in ens["cells"]])
```
Add all four to the printed summary and to the completion event's `observables` dict (alongside the existing fields).

- [ ] **Step 2: Add behavior_tests + expected_behavior to study.yaml**

For each of the four metrics add BOTH a `behavior_tests` entry and a mirrored `expected_behavior` entry, using the band from `karr_reference_values.json`:
```yaml
# behavior_tests: (append)
- name: emergent-mass-protein-fraction
  classification: secondary
  description: 'Emergent protein fraction of macromolecular dry mass ≈ 0.70 (Karr 2012 composition).'
  measure: {kind: derived_scalar, field: emergent_protein_fraction}
  pass_if:
    op: range
    low: 0.55
    high: 0.80
    provenance: {kind: experiment, note: 'Karr 2012 cell composition: protein 0.62 of total dry mass (0.705 renormalized over protein+DNA+RNA); emergent Σ(species×MW).'}
  requires_simulation: baseline
- name: single-cell-mrna-variation
  classification: secondary
  description: 'Cell-to-cell CV of total mRNA is nonzero and plausible (Fig 2 single-cell distributions).'
  measure: {kind: derived_scalar, field: single_cell_mrna_cv}
  pass_if:
    op: range
    low: 0.001
    high: 1.5
    provenance: {kind: experiment, note: 'Karr 2012 Fig 2: bursty low-copy single-cell mRNA → nonzero cell-to-cell variation.'}
  requires_simulation: baseline
```
Add analogous `emergent-mass-dna-fraction` (field `emergent_dna_fraction`, band 0.10–0.35) and `emergent-mass-rna-fraction` (field `emergent_rna_fraction`, band 0.05–0.25). Mirror all four in `expected_behavior` (with `en`, `observable`, `condition`, `rationale`, `measure`, `expect: {op: range, low, high}`), matching the file's existing style.

- [ ] **Step 3: Append a unit test for the run.py metric computation**

In `tests/test_validation.py`, import fig2's metric helpers if factored out, OR add a focused test that constructs a tiny `rows` fixture with `emergent_mass_fractions` and asserts `emergent_macro_fractions(rows[-1])` gives the expected fractions (the run.py wiring reuses the same helper, so testing the helper + a smoke of run_ensemble suffices; do NOT run the full 9h fig2 sim in the unit test).

- [ ] **Step 4: Backlog note + full suite**

`docs/FIDELITY_GAPS.md`: add a short "Comparison harness" section — Karr reference values + report-card gating metrics landed (emergent composition + single-cell CV in fig2), extensible to other figures; the emergent-composition gates are the current fidelity target for the mass/synthesis work.
Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--validation /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/ -q` (all pass).

- [ ] **Step 5: (optional) capture a baseline** — if a full fig2 run is affordable, run `workspace/studies/fig2-growth/sims/run.py` once and record the four new observables' values in the commit message / report so the baseline fidelity is captured. If too slow, skip (the report card computes it on demand).

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "fig2: gate emergent mass composition + single-cell mRNA CV against Karr reference"
```

---

## Self-Review

**Spec coverage:** reference dataset (T1); helpers (T1); fig2 emergent + CV derived scalars (T2 S1); behavior_tests/expected_behavior gates (T2 S2); resources (T1 S5); backlog (T2 S4); tests (both). Covered.

**Placeholder scan:** T1 concrete; T2's study.yaml edits show the exact YAML for 2 of 4 tests + explicit instructions for the other 2 + the run.py block — bands sourced from the dataset. No TBD.

**Type consistency:** `reference(id)->dict`, `emergent_macro_fractions(row)->{protein,DNA,RNA}`, `mrna_cv(list)->float`, `score(obs,ref)->{observed,target,band,pass,closeness}`. Derived-scalar field names (`emergent_protein_fraction` etc., `single_cell_mrna_cv`) match between run.py observables and study.yaml `measure.field`. Consistent.
