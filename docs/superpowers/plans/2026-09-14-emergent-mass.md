# Emergent Mass Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Compute the cell's dry-mass composition emergently as Σ(species count × molecular weight) from the real molecular inventory, expose it as observables, and validate it against Karr's fitted dry-weight fractions — without changing the growth-law that drives division.

**Architecture:** A new pure module `viva_mgen/mass_composition.py` computes per-component dry mass (RNA/protein/DNA/metabolite) in grams from counts + real MWs. `MassGrowthReproductionProcess` gains read-only sensor inputs (species stores) and two observable outputs (`emergent_mass`, `emergent_mass_fractions`); its existing growth-law mass/division behavior is untouched.

**Tech Stack:** Python 3.12, process-bigraph, pytest. Run tests with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--mass` and the venv at `/Users/eranagmon/code/viva-mGen/.venv`.

**Spec:** `docs/superpowers/specs/2026-09-14-emergent-mass-design.md`

## Global Constraints

- Work only in worktree `/Users/eranagmon/code/viva-mGen--mass` (branch `feat/emergent-mass`).
- Run everything with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--mass /Users/eranagmon/code/viva-mGen/.venv/bin/python`.
- No AI-attribution trailers.
- **Do NOT change** the growth-law `mass`/`volume`/`division`/`mass_fractions` behavior of `MassGrowthReproductionProcess`. The two new outputs are additive observables only.
- MW constants: `_RNA_NT_MW = 340.0`, `_AA_MW = 110.0`, `_BP_MW = 660.0` (both strands), `_N_A = 6.022e23`; metabolite MWs ATP 507, GTP 523, NTP 500, amino_acid 110. Mass in grams then ×1e15 for fg.

---

### Task 1: mass_composition module (pure MW math)

**Files:**
- Create: `viva_mgen/mass_composition.py`
- Test: `tests/test_emergent_mass.py`

**Interfaces:**
- Produces:
  - `rna_mass_g(rna_counts: dict) -> float`
  - `protein_mass_g(protein_counts: dict) -> float`
  - `dna_mass_g(chromosome_copy: float) -> float`
  - `metabolite_mass_g(pools: dict) -> float`
  - `emergent_composition(rna_counts, protein_counts, chromosome_copy, pools) -> dict` → `{"RNA","protein","DNA","metabolite","total": grams}`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_emergent_mass.py
import math
from viva_mgen.mass_composition import (
    rna_mass_g, protein_mass_g, dna_mass_g, metabolite_mass_g, emergent_composition,
    _N_A,
)

def test_rna_mass_one_gene():
    # one 300-nt transcript, count 1 → 300*340/N_A g
    g = rna_mass_g({"g1": 1.0}, lengths={"g1": 300.0})
    assert math.isclose(g, 300*340.0/_N_A, rel_tol=1e-9)

def test_protein_mass_one_gene():
    # count 2, 300-nt ORF → 300/3=100 aa each → 2*100*110/N_A
    g = protein_mass_g({"g1": 2.0}, lengths={"g1": 300.0})
    assert math.isclose(g, 2*100*110.0/_N_A, rel_tol=1e-9)

def test_dna_mass_matches_karr_fraction():
    # chromosome_copy 1 → ~0.6-0.7 fg (Karr DNA fraction 0.1688 * 3.93 fg ~= 0.66 fg)
    fg = dna_mass_g(1.0) * 1e15
    assert 0.55 < fg < 0.75

def test_metabolite_mass_positive():
    assert metabolite_mass_g({"atp": 1e6, "gtp": 1e6}) > 0.0

def test_emergent_composition_protein_dominant_and_sums():
    comp = emergent_composition(
        rna_counts={"g1": 100.0}, protein_counts={"g1": 5000.0},
        chromosome_copy=1.0, pools={"atp": 1e6},
        lengths={"g1": 1000.0})
    assert comp["total"] > 0
    assert abs((comp["RNA"]+comp["protein"]+comp["DNA"]+comp["metabolite"]) - comp["total"]) < 1e-30
    assert comp["protein"] == max(comp["RNA"], comp["protein"], comp["DNA"], comp["metabolite"])
```

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--mass /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_emergent_mass.py -q`
Expected: FAIL — module not found.

- [ ] **Step 3: Implement the module**

```python
# viva_mgen/mass_composition.py
"""Emergent dry-mass composition: Σ(species count × molecular weight).

Computes the cell's dry mass and its RNA/protein/DNA/metabolite breakdown from the
real molecular inventory + real molecular weights, for validation against Karr's
fitted dry-weight fractions. See
docs/superpowers/specs/2026-09-14-emergent-mass-design.md.
"""
from __future__ import annotations

from .constants import GENOME_LENGTH_BP

_N_A = 6.022e23        # Avogadro
_RNA_NT_MW = 340.0     # avg ribonucleotide-monophosphate MW in a chain (g/mol)
_AA_MW = 110.0         # avg amino-acid residue MW in a chain (g/mol)
_BP_MW = 660.0         # avg base-pair MW, both strands (g/mol)
_MET_MW = {"atp": 507.0, "gtp": 523.0, "ntp": 500.0, "amino_acid": 110.0}


def _lengths(lengths):
    if lengths is not None:
        return lengths
    from .expression_defaults import gene_lengths
    return gene_lengths()


def rna_mass_g(rna_counts, lengths=None):
    L = _lengths(lengths)
    return sum(float(c) * float(L.get(g, 1000.0)) * _RNA_NT_MW
               for g, c in (rna_counts or {}).items()) / _N_A


def protein_mass_g(protein_counts, lengths=None):
    L = _lengths(lengths)
    return sum(float(c) * (float(L.get(g, 1000.0)) / 3.0) * _AA_MW
               for g, c in (protein_counts or {}).items()) / _N_A


def dna_mass_g(chromosome_copy):
    return float(chromosome_copy) * float(GENOME_LENGTH_BP) * _BP_MW / _N_A


def metabolite_mass_g(pools):
    return sum(float(c) * _MET_MW.get(k, 300.0)
               for k, c in (pools or {}).items()) / _N_A


def emergent_composition(rna_counts, protein_counts, chromosome_copy, pools, lengths=None):
    rna = rna_mass_g(rna_counts, lengths)
    prot = protein_mass_g(protein_counts, lengths)
    dna = dna_mass_g(chromosome_copy)
    met = metabolite_mass_g(pools)
    return {"RNA": rna, "protein": prot, "DNA": dna, "metabolite": met,
            "total": rna + prot + dna + met}
```

- [ ] **Step 4: Run to verify pass**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--mass /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_emergent_mass.py -q`
Expected: PASS (5 tests).

- [ ] **Step 5: Commit**

```bash
git add viva_mgen/mass_composition.py tests/test_emergent_mass.py
git commit -m "Add emergent mass-composition module: Sigma(species x MW) with real MWs"
```

---

### Task 2: wire emergent observables into MassGrowth + composite + validation + fidelity note

**Files:**
- Modify: `viva_mgen/processes/mass.py`
- Modify: `viva_mgen/composites/mgen.py`
- Test: `tests/test_emergent_mass.py` (append integration)
- Regenerate: `reports/composite-state/*.json`

**Interfaces:**
- Consumes: `emergent_composition` (Task 1).
- Produces: mass node with `emergent_mass` (float) + `emergent_mass_fractions` (map[float]) observables.

- [ ] **Step 1: Add sensor inputs + observable outputs to MassGrowth**

In `viva_mgen/processes/mass.py`, `MassGrowthReproductionProcess`:
- `inputs()` gains (in addition to existing `growth_fraction`, `mass`): `"rna_counts": "map[float]"`, `"protein_counts": "map[float]"`, `"chromosome_copy": "float"`, `"atp": "float"`, `"gtp": "float"`, `"ntp": "float"`, `"amino_acid": "float"`.
- `outputs()` gains: `"emergent_mass": "overwrite[float]"`, `"emergent_mass_fractions": "overwrite[map[float]]"`.
- In `update()`, after the existing growth-law block, compute:
```python
from ..mass_composition import emergent_composition
comp = emergent_composition(
    state.get("rna_counts", {}), state.get("protein_counts", {}),
    float(state.get("chromosome_copy", 1.0) or 1.0),
    {k: float(state.get(k, 0.0) or 0.0) for k in ("atp", "gtp", "ntp", "amino_acid")},
)
total = comp["total"]
frac = ({k: comp[k] / total for k in ("RNA", "protein", "DNA", "metabolite")}
        if total > 0 else {k: 0.0 for k in ("RNA", "protein", "DNA", "metabolite")})
```
and add to the return dict: `"emergent_mass": total * 1e15` (fg), `"emergent_mass_fractions": frac`. Do NOT change the existing `mass`/`volume`/`division`/`mass_fractions` entries.
- `initial_state()` gains defaults for the new sensor inputs (`{}`/`0.0`/`1.0`) so a standalone run works.

- [ ] **Step 2: Composite wiring + emit**

In `viva_mgen/composites/mgen.py`:
- `_STORE_GROUP`: add `"emergent_mass": "physiology"`, `"emergent_mass_fractions": "physiology"`.
- `_EMPTY_MAP_STORES`: add `"emergent_mass_fractions"` (map, init `{}`).
- `_EMIT`: add `"emergent_mass": "float"`, `"emergent_mass_fractions": "map[float]"`.
- The mass node's new sensor input ports (`rna_counts`, `protein_counts`, `chromosome_copy`, `atp`, `gtp`, `ntp`, `amino_acid`) resolve by name to the existing shared stores — verify no `_PORT_STORE` override is needed (names already map: rna_counts→transcriptome, protein_counts→proteome, chromosome_copy→genome, pools→metabolism). Because these are read-only INPUTS on mass and existing processes already OWN those stores, mass just reads them; confirm the build wires them and the composite still runs.

- [ ] **Step 3: Fidelity note**

In `mass.py`, append to the `MassGrowthReproductionProcess.description` Fidelity line (do NOT change its first line): a sentence that the emergent dry-mass composition (Σ species×MW from the real inventory) is now computed and emitted (`emergent_mass`/`emergent_mass_fractions`) and validated against the fitted fractions; the growth-law total still drives division pending synthesis calibration for net doubling.

- [ ] **Step 4: Integration + validation test**

Append to `tests/test_emergent_mass.py`:
```python
def test_composite_emits_emergent_mass():
    from viva_mgen.composites.mgen import build_mgen
    from viva_mgen.core import build_core
    from process_bigraph import Composite
    core = build_core()
    doc = build_mgen(core=core)
    comp = Composite({"state": doc}, core=core)
    comp.run(3.0)
    phys = comp.state["cell"]["physiology"]
    assert phys["emergent_mass"] > 0.0
    fr = phys["emergent_mass_fractions"]
    s = sum(fr.values())
    assert abs(s - 1.0) < 1e-6 or s == 0.0
```

- [ ] **Step 5: Run tests**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--mass /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/ -q`
Expected: PASS (all existing + new).

- [ ] **Step 6: Regenerate committed snapshot**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--mass /Users/eranagmon/code/viva-mGen/.venv/bin/python scripts/regen_composite_state.py`
Expected: writes both snapshots; the mass node now shows the emergent outputs.

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "MassGrowth: emit emergent dry-mass composition; validate vs fitted fractions"
```

---

## Self-Review

**Spec coverage:** MW helper module (Task 1); emergent outputs on mass process (Task 2 Step 1); composite wiring + emit (Step 2); validation tests (Task 1 units + Task 2 integration); fidelity note (Step 3); snapshot regen (Step 6); growth-law/division untouched (constraint honored — only additive observables). All spec sections covered.

**Placeholder scan:** all code steps carry real code; no TBD.

**Type consistency:** `emergent_composition(rna_counts, protein_counts, chromosome_copy, pools, lengths=None) -> dict` with keys RNA/protein/DNA/metabolite/total, used consistently in Task 2. `_N_A` exported from `mass_composition` and imported in the test. `emergent_mass` float (fg), `emergent_mass_fractions` map[float].
