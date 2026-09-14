# Unified Chromosome (Phase 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Introduce a shared per-site chromosome-state structure and migrate the DNADamage→DNARepair subsystem onto it (per-site lesions that damage adds and repair clears from the same sites), with zero baseline behavior change.

**Architecture:** A pure `viva_mgen/chromosome_state.py` operates on a bin-indexed **lesion map** (`{bin(str) -> count(float)}`). `DNADamageReproductionProcess` emits per-site lesion additions; `DNARepairReproductionProcess` clears specific sites; both share one `lesion_map` store. The `lesions` scalar becomes the Σ of that structure. Later phases migrate the other DNA submodels (staged).

**Tech Stack:** Python 3.12, numpy, process-bigraph, pytest. Run with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--chromosome` and venv `/Users/eranagmon/code/viva-mGen/.venv`.

**Spec:** `docs/superpowers/specs/2026-09-14-unified-chromosome-design.md`

## Global Constraints

- Work only in worktree `/Users/eranagmon/code/viva-mGen--chromosome` (branch `feat/unified-chromosome`).
- Run everything with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--chromosome /Users/eranagmon/code/viva-mGen/.venv/bin/python`.
- No AI-attribution trailers.
- **Baseline non-regression:** with `damaging_agent=0` (the default), no lesions arise ⇒ `lesions` stays 0 and every existing result is unchanged. A test must assert this.
- `lesion_map` is `map[float]` keyed by `str(int bin)`, PRE-SEEDED to all-bins-zero so additive deltas land (a `map[float]` store drops deltas to absent keys).
- Repair can never drive a bin below 0 or repair more than the lesions present.

---

### Task 1: chromosome_state library (pure per-site ops)

**Files:**
- Create: `viva_mgen/chromosome_state.py`
- Test: `tests/test_chromosome_state.py`

**Interfaces:**
- Produces: `empty_lesion_map(n_bins) -> dict`; `add_lesions(rng, n_bins, count) -> dict`; `repair_sites(lesion_map, capacity, rng) -> (dict, float)`; `n_lesions(lesion_map) -> float`; `lesion_positions(lesion_map, genome_length_bp, n_bins) -> list[float]`.

- [ ] **Step 1: Write failing tests**

```python
# tests/test_chromosome_state.py
import numpy as np
from viva_mgen.chromosome_state import (
    empty_lesion_map, add_lesions, repair_sites, n_lesions, lesion_positions,
)

def test_empty_lesion_map():
    m = empty_lesion_map(10)
    assert len(m) == 10 and set(m.values()) == {0.0} and all(isinstance(k, str) for k in m)

def test_add_lesions_total_weight():
    rng = np.random.default_rng(0)
    d = add_lesions(rng, n_bins=100, count=5)
    assert sum(d.values()) == 5.0
    assert all(0 <= int(b) < 100 for b in d)

def test_add_lesions_zero():
    assert add_lesions(np.random.default_rng(0), 100, 0) == {}

def test_repair_never_exceeds_present():
    m = {"3": 2.0, "7": 1.0}   # 3 lesions present
    delta, repaired = repair_sites(m, capacity=10, rng=np.random.default_rng(1))
    assert repaired == 3.0
    # applying delta zeroes the map, never negative
    for b, dv in delta.items():
        assert m[b] + dv >= -1e-9
    assert sum(delta.values()) == -3.0

def test_repair_partial():
    m = {"3": 5.0}
    delta, repaired = repair_sites(m, capacity=2, rng=np.random.default_rng(2))
    assert repaired == 2.0 and delta == {"3": -2.0}

def test_n_lesions():
    assert n_lesions({"1": 2.0, "2": 0.0, "5": 3.0}) == 5.0

def test_lesion_positions():
    pos = lesion_positions({"0": 1.0, "50": 2.0, "3": 0.0}, genome_length_bp=1000.0, n_bins=100)
    assert 0.0 in pos and 500.0 in pos and len(pos) == 2   # bin 3 has 0 count -> excluded
```

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--chromosome /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_chromosome_state.py -q`
Expected: FAIL — module missing.

- [ ] **Step 3: Implement the module**

```python
# viva_mgen/chromosome_state.py
"""Shared per-site chromosome state (bin-indexed) — the viva-native stand-in for
Karr 2012's CircularSparseMat. Phase 1 implements the LESION layer: a lesion map
{bin(str) -> count(float)} that DNADamage adds to and DNARepair clears from,
shared as one structure. Later phases add bound-protein footprints, per-region
linking number, and polymerized regions on the same bin index (see
docs/superpowers/specs/2026-09-14-unified-chromosome-design.md).
"""
from __future__ import annotations


def empty_lesion_map(n_bins: int) -> dict:
    """Pre-seeded all-zero lesion map (so additive map[float] deltas land)."""
    return {str(b): 0.0 for b in range(int(n_bins))}


def add_lesions(rng, n_bins: int, count) -> dict:
    """Additive delta placing ``count`` independent lesions at random bins."""
    count = int(count)
    if count <= 0:
        return {}
    delta: dict = {}
    for b in rng.integers(0, int(n_bins), size=count):
        delta[str(int(b))] = delta.get(str(int(b)), 0.0) + 1.0
    return delta


def n_lesions(lesion_map) -> float:
    return float(sum(float(v) for v in (lesion_map or {}).values()))


def repair_sites(lesion_map, capacity, rng):
    """Repair up to ``capacity`` lesion instances from bins with count>0, weighted
    by count. Returns (negative additive delta, number repaired). Never repairs
    more than present, never drives a bin below zero."""
    present = {b: float(v) for b, v in (lesion_map or {}).items() if float(v) > 0.0}
    total = sum(present.values())
    cap = int(min(float(capacity), total))
    if cap <= 0:
        return {}, 0.0
    # expand to per-instance bin list, sample `cap` without replacement
    bins = []
    for b, v in present.items():
        bins.extend([b] * int(v))
    chosen = rng.choice(len(bins), size=cap, replace=False)
    delta: dict = {}
    for i in chosen:
        b = bins[int(i)]
        delta[b] = delta.get(b, 0.0) - 1.0
    return delta, float(cap)


def lesion_positions(lesion_map, genome_length_bp: float, n_bins: int) -> list:
    """bp coordinates (bin centers) of currently-damaged bins."""
    bp_per_bin = float(genome_length_bp) / int(n_bins)
    return [int(b) * bp_per_bin for b, v in (lesion_map or {}).items() if float(v) > 0.0]
```

- [ ] **Step 4: Run to verify pass**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--chromosome /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_chromosome_state.py -q`
Expected: PASS (7 tests).

- [ ] **Step 5: Commit**

```bash
git add viva_mgen/chromosome_state.py tests/test_chromosome_state.py
git commit -m "Add shared per-site chromosome state (lesion layer): the CircularSparseMat stand-in"
```

---

### Task 2: migrate DNADamage/DNARepair onto the shared lesion structure

**Files:**
- Modify: `viva_mgen/processes/dna.py` (DNADamage, DNARepair)
- Modify: `viva_mgen/composites/mgen.py` (lesion_map store + pre-seed)
- Test: `tests/test_chromosome_state.py` (append process + integration)
- Regenerate: `reports/composite-state/*.json`

**Interfaces:**
- Consumes: `chromosome_state` (Task 1).

- [ ] **Step 1: DNADamage emits per-site lesions**

In `DNADamageReproductionProcess` (dna.py): add config `n_bins` (`{"_type":"integer","_default":580}`); add `lesion_map` to `inputs()` and `outputs()` (`map[float]`); in `__init__` keep the rng. In `update()`, after computing `new_lesions` (unchanged Poisson), add:
```python
from ..chromosome_state import add_lesions
...
        return {"lesions": new_lesions,
                "lesion_map": add_lesions(self._rng, int(self.config["n_bins"]), int(new_lesions))}
```
`initial_state()` gains `"lesion_map": {}`.

- [ ] **Step 2: DNARepair clears specific sites**

In `DNARepairReproductionProcess` (dna.py): add `lesion_map` to `inputs()` and `outputs()` (`map[float]`). In `update()`, replace the abstract `repaired = min(lesions, capacity, atp/atp_per)` clearing with site-based repair:
```python
from ..chromosome_state import repair_sites, n_lesions
...
        lesion_map = state.get("lesion_map", {}) or {}
        atp_cap = atp / atp_per if atp_per > 0 else n_lesions(lesion_map)
        capacity = min(self.config["repair_rate"] * enzyme * interval, atp_cap)
        delta, repaired = repair_sites(lesion_map, capacity, self._rng)
        return {"lesions": -repaired, "atp": -repaired * atp_per, "lesion_map": delta}
```
Add `self._rng = np.random.default_rng(int(self.config.get("seed", 5)))` in `__init__` (DNARepair currently has no rng — add a `seed` config default 5) and `lesion_map` to `initial_state()`.

- [ ] **Step 3: Composite wiring + pre-seed**

In `viva_mgen/composites/mgen.py`:
- `_STORE_GROUP`: add `"lesion_map": "genome"`.
- Pre-seed it to all-bins-zero: in the map-init logic (where `_EMPTY_MAP_STORES`/`_MAP_STORES` decide the initial value), special-case `lesion_map` to `empty_lesion_map(580)` (import `empty_lesion_map` from `..chromosome_state`) so additive deltas land. (Mirror how demand__<pool> was special-cased for the allocator.)
- Add `lesion_map` to `_MAP_STORES` so it's treated as a map.
- DNADamage/DNARepair `lesion_map` ports auto-wire by name to the store.

- [ ] **Step 4: Fidelity notes**

Append to DNADamage and DNARepair `description` (first line unchanged): a sentence that lesions are now tracked per-site on the shared chromosome structure (`lesion_map`) that damage adds to and repair clears from the same sites — the first consumer of the unified per-site chromosome (remaining DNA submodels staged).

- [ ] **Step 5: Process + integration tests (append to tests/test_chromosome_state.py)**

```python
def test_damage_populates_lesion_map():
    from viva_mgen.processes.dna import DNADamageReproductionProcess
    from viva_mgen.core import build_core
    p = DNADamageReproductionProcess(config={"agent_rate": 1.0, "n_bins": 100}, core=build_core())
    out = p.update({"damaging_agent": 50.0}, 1.0)
    from viva_mgen.chromosome_state import n_lesions
    assert n_lesions(out["lesion_map"]) == out["lesions"] and out["lesions"] > 0

def test_damage_zero_agent_no_lesions():
    from viva_mgen.processes.dna import DNADamageReproductionProcess
    from viva_mgen.core import build_core
    p = DNADamageReproductionProcess(config={"base_rate": 0.0, "n_bins": 100}, core=build_core())
    out = p.update({"damaging_agent": 0.0}, 1.0)
    assert out["lesions"] == 0.0 and out["lesion_map"] == {}

def test_repair_clears_sites():
    from viva_mgen.processes.dna import DNARepairReproductionProcess
    from viva_mgen.core import build_core
    p = DNARepairReproductionProcess(config={"repair_rate": 1.0}, core=build_core())
    out = p.update({"lesion_map": {"3": 2.0, "9": 1.0}, "repair_enzyme": 100.0, "atp": 1e9}, 1.0)
    assert out["lesions"] < 0 and sum(out["lesion_map"].values()) == out["lesions"]

def test_composite_baseline_no_lesions():
    from viva_mgen.composites.mgen import build_mgen
    from viva_mgen.core import build_core
    from process_bigraph import Composite
    core = build_core()
    comp = Composite({"state": build_mgen(core=core)}, core=core)
    comp.run(3.0)
    lm = comp.state["cell"]["genome"]["lesion_map"]
    assert sum(lm.values()) == 0.0   # damaging_agent defaults 0 -> no lesions (non-regression)
```

- [ ] **Step 6: Run full suite + regenerate snapshot**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--chromosome /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/ -q` (all pass).
Then: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--chromosome /Users/eranagmon/code/viva-mGen/.venv/bin/python scripts/regen_composite_state.py`.

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "DNADamage/DNARepair: per-site lesions on the shared chromosome structure"
```

---

## Self-Review

**Spec coverage:** shared structure lib (T1); damage per-site (T2 S1); repair per-site (T2 S2); composite pre-seed (T2 S3); fidelity notes (T2 S4); baseline non-regression test (T2 S5); snapshot (T2 S6). Covered.

**Placeholder scan:** all steps carry real code; no TBD.

**Type consistency:** `add_lesions(rng,n_bins,count)->dict`, `repair_sites(lesion_map,capacity,rng)->(dict,float)`, `n_lesions(map)->float`, `empty_lesion_map(n)->dict`; `lesion_map` is `map[float]` keyed by str(bin); DNARepair gains a `seed` config. Consistent across tasks.
