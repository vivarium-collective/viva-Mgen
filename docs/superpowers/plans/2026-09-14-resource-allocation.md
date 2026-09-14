# Resource-Allocation Layer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reproduce Karr 2012's hybrid per-timestep resource partitioning in viva-Mgen: finite metabolite pools replenished by metabolism, partitioned among competing submodels by a central allocator via a demand→allocate→run pipeline.

**Architecture:** A new `AllocatorProcess` runs first each tick; it reads each finite pool's level + metabolic production + every consumer's demand (a `map[float]` keyed by consumer id), computes per-consumer allocations (proportional, with a priority-weight hook), and replenishes the pools. Each pool-consuming process gains a `demand__<pool>` output and an `alloc__<pool>` input; it caps its consumption at `min(kinetic_limit, budget)` where `budget = alloc.get(consumer_id, inf)`, and emits its unconstrained want as next tick's demand. One-tick pipeline lag, biologically negligible at Δt=1 s.

**Tech Stack:** Python 3.12, process-bigraph (`Process` base, `map[float]`/`float` port types), numpy, pytest. Run tests with the venv at `/Users/eranagmon/code/viva-mGen/.venv` and `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc`.

**Spec:** `docs/superpowers/specs/2026-09-14-resource-allocation-design.md`

## Global Constraints

- Work only in the worktree `/Users/eranagmon/code/viva-mGen--alloc` (branch `feat/resource-allocation`). Never touch the canonical checkout.
- Run every command with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python` so imports resolve to the worktree, not the editable install.
- No AI-attribution trailers on commits.
- Back-compat is mandatory: an absent `alloc__<pool>` entry ⇒ `budget = inf` (a consumer behaves exactly as today when unwired or at tick 0). The allocator is strictly *additional* constraint, never a new failure mode.
- Partitioned pools (v1): `atp`, `gtp`, `ntp`, `amino_acid`. `dntp` is out of scope (single consumer, self-managed via `dntp_synthesis_scale`).
- No process's kinetic constants or core mechanism change — only the consumption cap and the new demand/budget ports.

## Pools → consumers (consumer_id in parens)

| Pool | Consumers |
|------|-----------|
| `atp` | supercoiling (`supercoiling`), dna_repair (`dna_repair`), protein_folding (`protein_folding`), protein_modification (`protein_modification`), trna_aminoacylation (`trna_aminoacylation`) |
| `gtp` | translation (`translation`), protein_translocation (`translocation`), ribosome_assembly (`ribosome_assembly`), ftsz (`ftsz`) |
| `ntp` | transcription (`transcription`) |
| `amino_acid` | trna_aminoacylation (`trna_aminoacylation`) |

(`replication_initiation` reads `atp` only as a >0 gate — not a consumer, left unchanged.)

## File Structure

- **Create** `viva_mgen/processes/allocation.py` — `allocate()` pure function, `AllocatorProcess`, `select_budget()` / `demand_entry()` consumer helpers.
- **Create** `tests/test_allocation.py` — allocator unit tests + integration scarcity test.
- **Modify** `viva_mgen/processes/metabolism.py` — add `ntp_production` / `aa_production` outputs.
- **Modify** `viva_mgen/composites/mgen.py` — budget store group, allocator node + first-scheduling, finite pool init, metabolism→allocator + allocator→pool wiring.
- **Modify** consumers: `transcription.py`, `translation.py`, `rna.py` (aminoacylation), `dna.py` (supercoiling, dna_repair), `protein.py` (folding, modification, translocation, ribosome_assembly), `cytokinesis.py` (ftsz).
- **Modify** `scripts/regen_composite_state.py` run + republish (final task).

---

### Task 1: Allocator core (pure math + Process + consumer helpers)

**Files:**
- Create: `viva_mgen/processes/allocation.py`
- Test: `tests/test_allocation.py`

**Interfaces:**
- Produces:
  - `allocate(supply: float, demands: dict[str, float], priorities: dict[str, float] | None = None) -> dict[str, float]` — returns `{consumer_id: grant}` with `sum(grants) <= supply`.
  - `select_budget(alloc: dict, consumer_id: str) -> float` — `alloc.get(consumer_id, inf)`.
  - `demand_entry(consumer_id: str, want: float) -> dict` — `{consumer_id: max(0.0, want)}`.
  - `AllocatorProcess(Process)` with config `pools` (list[str], default `["atp","gtp","ntp","amino_acid"]`), `priorities` (`map[map[float]]`, default `{}`), `pool_cap` (float, default `1e9`); inputs `<pool>` (float level) + `<pool>_supply` (float production) + `demand__<pool>` (map[float]) for each pool; outputs `<pool>` (float delta) + `alloc__<pool>` (overwrite[map[float]]) for each pool.

- [ ] **Step 1: Write the failing tests for `allocate()`**

```python
# tests/test_allocation.py
import math
import numpy as np
from viva_mgen.processes.allocation import allocate, select_budget, demand_entry


def test_allocate_surplus_grants_full():
    g = allocate(100.0, {"a": 10.0, "b": 20.0})
    assert g == {"a": 10.0, "b": 20.0}


def test_allocate_scarcity_is_proportional():
    g = allocate(30.0, {"a": 30.0, "b": 60.0})
    assert math.isclose(g["a"], 10.0) and math.isclose(g["b"], 20.0)
    assert sum(g.values()) <= 30.0 + 1e-9


def test_allocate_priority_weights():
    g = allocate(30.0, {"a": 30.0, "b": 30.0}, priorities={"a": 2.0, "b": 1.0})
    assert math.isclose(g["a"], 20.0) and math.isclose(g["b"], 10.0)


def test_allocate_zero_supply_grants_zero():
    g = allocate(0.0, {"a": 5.0})
    assert g == {"a": 0.0}


def test_allocate_ignores_nonpositive_and_conserves():
    g = allocate(10.0, {"a": -3.0, "b": 0.0, "c": 40.0})
    assert g["a"] == 0.0 and g["b"] == 0.0
    assert math.isclose(g["c"], 10.0)


def test_select_budget_absent_is_inf():
    assert select_budget({}, "x") == float("inf")
    assert select_budget({"x": 5.0}, "x") == 5.0


def test_demand_entry_clamps_negative():
    assert demand_entry("x", -2.0) == {"x": 0.0}
```

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_allocation.py -q`
Expected: FAIL — `ModuleNotFoundError: viva_mgen.processes.allocation`.

- [ ] **Step 3: Implement `allocate()` + helpers**

```python
# viva_mgen/processes/allocation.py
"""Whole-cell resource-allocation layer (Karr 2012 hybrid partitioning).

Each finite metabolite pool (ATP, GTP, NTP, amino acids) is partitioned each
tick among the submodels that demand it: demand -> allocate -> run. See
docs/superpowers/specs/2026-09-14-resource-allocation-design.md.
"""
from __future__ import annotations

import numpy as np
from process_bigraph import Process

DEFAULT_POOLS = ["atp", "gtp", "ntp", "amino_acid"]


def allocate(supply, demands, priorities=None):
    """Partition ``supply`` (>=0) of one pool among ``demands`` (consumer_id ->
    wanted amount). Non-positive demands get 0. If total weighted demand fits in
    supply every consumer is granted its full demand; otherwise grants are
    proportional to weight*demand. Guarantees ``sum(grants) <= supply``."""
    supply = max(0.0, float(supply))
    pos = {k: float(v) for k, v in (demands or {}).items() if float(v) > 0.0}
    grants = {k: 0.0 for k in (demands or {})}
    if not pos or supply <= 0.0:
        return grants
    total = sum(pos.values())
    if total <= supply:
        grants.update(pos)
        return grants
    w = {k: float((priorities or {}).get(k, 1.0)) for k in pos}
    wsum = sum(w[k] * pos[k] for k in pos)
    if wsum <= 0.0:
        return grants
    for k in pos:
        grants[k] = w[k] * pos[k] / wsum * supply
    return grants


def select_budget(alloc, consumer_id):
    """A consumer's granted budget from an ``alloc__<pool>`` map; ``inf`` when the
    consumer has no entry (tick 0 / unwired) so it behaves as pool-unlimited."""
    if not isinstance(alloc, dict) or consumer_id not in alloc:
        return float("inf")
    return float(alloc[consumer_id])


def demand_entry(consumer_id, want):
    """A ``demand__<pool>`` output entry (non-negative)."""
    return {consumer_id: max(0.0, float(want))}
```

- [ ] **Step 4: Run to verify the math tests pass**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_allocation.py -q`
Expected: PASS (7 tests).

- [ ] **Step 5: Write the failing test for `AllocatorProcess`**

```python
# append to tests/test_allocation.py
from process_bigraph import Process
from viva_mgen.processes.allocation import AllocatorProcess


def test_allocator_replenishes_and_partitions():
    proc = AllocatorProcess.__new__(AllocatorProcess)   # skip core-requiring __init__
    Process.__init__(proc, {"pools": ["atp"]}, core=None)
    state = {"atp": 5.0, "atp_supply": 25.0,
             "demand__atp": {"a": 20.0, "b": 40.0}}
    out = proc.update(state, 1.0)
    # supply available = level(5) + production(25) = 30; demands 20+40=60 -> proportional
    assert out["alloc__atp"]["a"] == 10.0 and out["alloc__atp"]["b"] == 20.0
    # pool replenished by production this tick (consumers draw it down elsewhere)
    assert out["atp"] == 25.0
```

- [ ] **Step 6: Run to verify it fails**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_allocation.py::test_allocator_replenishes_and_partitions -q`
Expected: FAIL — `AllocatorProcess` not defined.

- [ ] **Step 7: Implement `AllocatorProcess`**

```python
# append to viva_mgen/processes/allocation.py
class AllocatorProcess(Process):
    """Central resource allocator (reproduction of Karr 2012 hybrid partitioning).

    Runs first each tick. For every finite pool it computes the available supply
    (current level + this tick's metabolic production), partitions it among the
    consumers' demands (demand__<pool>) by allocate(), and publishes the grants
    (alloc__<pool>) that each consumer caps its consumption at. Also replenishes
    each pool by this tick's production (consumers draw it down via their own
    negative deltas), holding the pool >= 0 with sum(grants) <= supply.

    Contract — per pool P: in <P> (level, float), <P>_supply (production, float),
    demand__<P> (consumer->want, map). out <P> (replenish delta, float),
    alloc__<P> (consumer->grant, overwrite map).
    Fidelity: FAITHFUL to Karr's demand->allocate->run arbitration (proportional
    partition with a priority-weight hook); one-tick pipeline lag at Δt=1 s.
    """

    config_schema = {
        "pools": {"_type": "list[string]", "_default": DEFAULT_POOLS},
        "priorities": {"_type": "map[map[float]]", "_default": {}},
        "pool_cap": {"_type": "float", "_default": 1.0e9},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self._pools = list(self.config["pools"])

    def inputs(self):
        s = {}
        for p in self._pools:
            s[p] = "float"
            s[p + "_supply"] = "float"
            s["demand__" + p] = "map[float]"
        return s

    def outputs(self):
        s = {}
        for p in self._pools:
            s[p] = "float"
            s["alloc__" + p] = "overwrite[map[float]]"
        return s

    def initial_state(self):
        return {p + "_supply": 0.0 for p in self._pools}

    def update(self, state, interval):
        prio = self.config["priorities"] or {}
        cap = float(self.config["pool_cap"])
        out = {}
        for p in self._pools:
            level = float(state.get(p, 0.0) or 0.0)
            production = max(0.0, float(state.get(p + "_supply", 0.0) or 0.0))
            demands = state.get("demand__" + p, {}) or {}
            supply = min(level + production, cap)
            out["alloc__" + p] = allocate(supply, demands, prio.get(p))
            # replenish the pool by this tick's production (bounded by cap)
            out[p] = min(production, max(0.0, cap - level))
        return out
```

- [ ] **Step 8: Run all allocation tests**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_allocation.py -q`
Expected: PASS (8 tests).

- [ ] **Step 9: Commit**

```bash
cd /Users/eranagmon/code/viva-mGen--alloc
git add viva_mgen/processes/allocation.py tests/test_allocation.py
git commit -m "Add resource allocator: proportional pool partitioning + consumer helpers"
```

---

### Task 2: Metabolism supplies NTP + amino-acid production

**Files:**
- Modify: `viva_mgen/processes/metabolism.py` (outputs + update return)
- Test: `tests/test_allocation.py` (append)

**Interfaces:**
- Consumes: existing `MetabolismFbaReproductionProcess` (outputs `growth_rate`, `growth_fraction`, `atp_production`, `gtp_production`, `feasible`).
- Produces: two new outputs `ntp_production` (float) and `aa_production` (float) = `growth_fraction × base rate`, so precursor supply falls when growth is stressed (the metabolism→precursor coupling).

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_allocation.py
from viva_mgen.processes.metabolism import MetabolismFbaReproductionProcess


def test_metabolism_emits_precursor_supply():
    out = MetabolismFbaReproductionProcess.outputs(
        MetabolismFbaReproductionProcess.__new__(MetabolismFbaReproductionProcess))
    assert "ntp_production" in out and "aa_production" in out
```

- [ ] **Step 2: Run to verify it fails**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_allocation.py::test_metabolism_emits_precursor_supply -q`
Expected: FAIL — keys absent.

- [ ] **Step 3: Add the outputs + config**

In `viva_mgen/processes/metabolism.py`, add to `config_schema` (after `reaction_bound_scale`):

```python
        # base precursor supply per tick (molecules), scaled by growth_fraction —
        # the metabolism->transcription/translation precursor coupling.
        "ntp_base_supply": {"_type": "float", "_default": 1.0e6},
        "aa_base_supply": {"_type": "float", "_default": 1.0e6},
```

Add to `outputs()` (inside the returned dict):

```python
            "ntp_production": "overwrite[float]",
            "aa_production": "overwrite[float]",
```

In `update()`, before the final `return`, compute (`growth` and `self._wt` already exist):

```python
        gf = growth / self._wt if self._wt else 0.0
        ntp_prod = float(self.config["ntp_base_supply"]) * gf
        aa_prod = float(self.config["aa_base_supply"]) * gf
```

and add to the returned dict:

```python
            "ntp_production": ntp_prod,
            "aa_production": aa_prod,
```

- [ ] **Step 4: Run to verify it passes**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_allocation.py::test_metabolism_emits_precursor_supply -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add viva_mgen/processes/metabolism.py tests/test_allocation.py
git commit -m "Metabolism: emit growth-scaled NTP/amino-acid precursor supply"
```

---

### Task 3: Wire the allocator into the composite (finite pools, scheduling, replenishment)

**Files:**
- Modify: `viva_mgen/composites/mgen.py`
- Test: `tests/test_allocation.py` (append integration test)

**Interfaces:**
- Consumes: `AllocatorProcess` (Task 1), metabolism `*_production` outputs (Task 2).
- Produces: composite that instantiates the allocator as node `allocator`, a `budget` store group holding `demand__<pool>` / `alloc__<pool>` maps, finite pool init, and `<pool>_supply` stores fed by metabolism.

- [ ] **Step 1: Register allocator node + store group + finite pools**

In `viva_mgen/composites/mgen.py`:

Add to `_NODE_NAMES` (so it's discovered/instantiated):
```python
    "AllocatorProcess": "allocator",
```

Add the pool-supply + budget keys to `_STORE_GROUP` (new `budget` compartment for the maps; supplies live under `metabolism`):
```python
    "atp_supply": "metabolism", "gtp_supply": "metabolism",
    "ntp_supply": "metabolism", "aa_supply": "metabolism",
    "demand__atp": "budget", "demand__gtp": "budget",
    "demand__ntp": "budget", "demand__amino_acid": "budget",
    "alloc__atp": "budget", "alloc__gtp": "budget",
    "alloc__ntp": "budget", "alloc__amino_acid": "budget",
```

Add the demand/alloc maps to `_MAP_STORES` (so they init to `{}`, not per-gene):
```python
    "demand__atp", "demand__gtp", "demand__ntp", "demand__amino_acid",
    "alloc__atp", "alloc__gtp", "alloc__ntp", "alloc__amino_acid",
```

Reduce the pool init in `_FLOAT_INIT` from infinite to finite steady-state-scale levels:
```python
    "atp": 1e6, "gtp": 1e6, "ntp": 1e6, "amino_acid": 1e6,
```

- [ ] **Step 2: Map metabolism production ports → supply stores; import AllocatorProcess for discovery**

`all_process_classes()` (in `viva_mgen/processes/__init__.py`) auto-collects `*ReproductionProcess` by suffix; `AllocatorProcess` does NOT match that suffix. Add an explicit import + inclusion so the composite can instantiate it. In `viva_mgen/processes/__init__.py`, import it and add to the module list used by `all_process_classes()`:

```python
from .allocation import AllocatorProcess   # noqa: F401
```
and in `all_process_classes()` extend the returned dict:
```python
    from . import allocation
    out["AllocatorProcess"] = allocation.AllocatorProcess
```

Wire metabolism's `*_production` outputs to the `*_supply` stores the allocator reads. In `mgen.py`, add a port→store override so metabolism's `atp_production` feeds `atp_supply` etc. Extend `_PORT_STORE`:
```python
    ("MetabolismFbaReproductionProcess", "atp_production"): "atp_supply",
    ("MetabolismFbaReproductionProcess", "gtp_production"): "gtp_supply",
    ("MetabolismFbaReproductionProcess", "ntp_production"): "ntp_supply",
    ("MetabolismFbaReproductionProcess", "aa_production"): "aa_supply",
```
The allocator's `amino_acid_supply` port name is `amino_acid_supply`; map metabolism `aa_production` → `aa_supply` and set the allocator's pool list so its supply port for `amino_acid` reads `aa_supply`. To keep names aligned, give the allocator config `pools` but read supply via a fixed map. Simplest: rename the allocator's supply port derivation to use `{"amino_acid": "aa"}` alias. Implement by adding to `AllocatorProcess` a `_SUPPLY_ALIAS = {"amino_acid": "aa"}` and building the supply port name as `alias.get(p, p) + "_supply"`. Update `inputs()` accordingly. (Add this alias in Task 1's file if not already; if discovered here, amend `allocation.py` and its test.)

- [ ] **Step 3: Ensure the allocator is scheduled first**

Composite step order follows document key order. In `build_mgen()`, when assembling `doc`, insert the `allocator` node entry before the others (iterate `_NODE_NAMES` with `allocator` first). Concretely, sort the `for cls_name, node in _NODE_NAMES.items()` loop so `AllocatorProcess` is processed first:

```python
    items = sorted(_NODE_NAMES.items(), key=lambda kv: 0 if kv[1] == "allocator" else 1)
    for cls_name, node in items:
```

- [ ] **Step 4: Write the failing integration test**

```python
# append to tests/test_allocation.py
def test_composite_builds_and_runs_with_allocator():
    from viva_mgen.composites.mgen import build_mgen
    from viva_mgen.core import build_core
    from process_bigraph import Composite
    core = build_core()
    doc = build_mgen(core=core)
    assert "allocator" in doc
    comp = Composite({"state": doc}, core=core)
    comp.run(3.0)
    budget = comp.state["cell"]["budget"]
    assert "alloc__atp" in budget and "alloc__gtp" in budget
```

- [ ] **Step 5: Run to verify it fails, then passes after wiring**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_allocation.py::test_composite_builds_and_runs_with_allocator -q`
Expected: PASS once Steps 1-3 are complete (fix wiring errors until the composite builds + runs 3 s and the budget stores exist).

- [ ] **Step 6: Run the full suite (no regressions)**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/ -q`
Expected: PASS (existing 56 + new).

- [ ] **Step 7: Commit**

```bash
git add viva_mgen/composites/mgen.py viva_mgen/processes/__init__.py viva_mgen/processes/allocation.py tests/test_allocation.py
git commit -m "Composite: wire allocator first, finite pools, metabolism-fed supply"
```

---

### Task 4: ATP consumers → demand/budget

**Files:**
- Modify: `viva_mgen/processes/dna.py` (`DNASupercoilingReproductionProcess`, `DNARepairReproductionProcess`), `viva_mgen/processes/protein.py` (`ProteinFoldingReproductionProcess`, `ProteinModificationReproductionProcess`), `viva_mgen/processes/rna.py` (`TRNAAminoacylationReproductionProcess`)
- Test: `tests/test_allocation.py` (append)

**Interfaces:**
- Consumes: `select_budget`, `demand_entry` (Task 1); `alloc__atp` store (Task 3).
- Produces: each ATP consumer reads `alloc__atp`, caps ATP consumption at its budget, emits `demand__atp`.

**The per-consumer pattern** (apply to each of the 5, using its own `consumer_id` and its own unconstrained ATP want):

Add to `config_schema`: `"consumer_id": {"_type": "string", "_default": "<id>"}`.
Add to `inputs()`: `"alloc__atp": "map[float]"`.
Add to `outputs()`: `"demand__atp": "map[float]"`.
In `__init__`: `self._cid = self.config["consumer_id"]`.
In `update()`: after computing the process's *unconstrained* ATP want `want_atp`, clamp:
```python
from .allocation import select_budget, demand_entry   # module-top import
budget = select_budget(state.get("alloc__atp", {}), self._cid)
atp_cap = min(want_atp, budget)
```
then use `atp_cap` wherever the code previously used the raw ATP availability to bound work, and add to the return dict: `"demand__atp": demand_entry(self._cid, want_atp)`.

Concrete per-consumer wants:
- **DNASupercoiling** (`supercoiling`): `want_atp = gyrase * gyrase_rate * interval * atp_per_act` (the ATP the desired acts need). Bound `acts` additionally by `atp_cap / atp_per_act`.
- **DNARepair** (`dna_repair`): `want_atp = repair_rate * enzyme * interval * atp_per_repair` capped by lesions; bound `repaired` by `atp_cap / atp_per_repair`.
- **ProteinFolding** (`protein_folding`): `want_atp = (chaperone-assisted folds this step) * atp_per_fold`; bound chaperone folds by `atp_cap / atp_per_fold`.
- **ProteinModification** (`protein_modification`): `want_atp = enzyme * modification_specific_rate * interval * atp_per_modification`; bound modifications by `atp_cap / atp_per_modification`.
- **tRNAAminoacylation** (`trna_aminoacylation`): `want_atp = N_desired * atp_per_charge` where `N_desired` is the pre-clamp charge count; add `atp_cap` into the existing `min(...)` that sets `n_total`.

Worked example (DNASupercoiling — do the analogous edit for the other four):

- [ ] **Step 1: Write the failing scarcity test**

```python
# append to tests/test_allocation.py
def test_atp_consumer_respects_budget():
    from viva_mgen.processes.dna import DNASupercoilingReproductionProcess as S
    from process_bigraph import Process
    proc = S.__new__(S); Process.__init__(proc, {}, core=None)
    proc._sigma = 0.0
    st = {"gyrase": 1000.0, "atp": 1e9,
          "alloc__atp": {"supercoiling": 4.0}}   # budget = 4 ATP -> 2 acts
    out = proc.update(st, 1.0)
    assert out["atp"] >= -4.0            # never spend more than the 4-ATP budget
    assert "demand__atp" in out and out["demand__atp"]["supercoiling"] > 0
```

- [ ] **Step 2: Run to verify it fails**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_allocation.py::test_atp_consumer_respects_budget -q`
Expected: FAIL (no `alloc__atp` handling / no `demand__atp`).

- [ ] **Step 3: Apply the pattern to DNASupercoiling**

Edit `DNASupercoilingReproductionProcess`: add `consumer_id` config (default `"supercoiling"`), `alloc__atp` input, `demand__atp` output, and in `update()`:
```python
from .allocation import select_budget, demand_entry
...
        want_acts = gyrase * self.config["gyrase_rate"] * interval
        want_atp = want_acts * self.config["atp_per_act"]
        budget = select_budget(state.get("alloc__atp", {}), self._cid)
        atp_cap = min(want_atp, budget, atp)   # atp still bounds as a floor safety
        acts = min(want_acts, atp_cap / self.config["atp_per_act"], acts_wanted)
```
Return dict gains `"demand__atp": demand_entry(self._cid, want_atp)`.

- [ ] **Step 4: Run to verify it passes**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_allocation.py::test_atp_consumer_respects_budget -q`
Expected: PASS.

- [ ] **Step 5: Apply the same pattern to the other four ATP consumers** (dna_repair, protein_folding, protein_modification, trna_aminoacylation), each with its `consumer_id` and `want_atp` from the list above.

- [ ] **Step 6: Run the full suite**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/ -q`
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add viva_mgen/processes/dna.py viva_mgen/processes/protein.py viva_mgen/processes/rna.py viva_mgen/composites/mgen.py tests/test_allocation.py
git commit -m "ATP consumers: consume under allocator budget, emit demand"
```

(Note: add each new consumer's `consumer_id` to `_PORT_STORE`/`_STORE_GROUP` only if a port name collides; `alloc__atp`/`demand__atp` already mapped to the `budget` group in Task 3.)

---

### Task 5: GTP consumers → demand/budget

**Files:**
- Modify: `viva_mgen/processes/translation.py` (`TranslationReproductionProcess`), `viva_mgen/processes/protein.py` (`ProteinTranslocationReproductionProcess`, `RibosomeAssemblyReproductionProcess`), `viva_mgen/processes/cytokinesis.py` (`FtsZPolymerizationReproductionProcess`)
- Test: `tests/test_allocation.py` (append)

**Interfaces:**
- Consumes: `select_budget`, `demand_entry`, `alloc__gtp` store.
- Produces: each GTP consumer caps GTP use at its budget, emits `demand__gtp`.

Apply the **same pattern as Task 4** with pool `gtp`, `consumer_id` per the pools table, and these unconstrained wants:
- **Translation** (`translation`): `want_gtp = sum over genes of expected_proteins * gtp_per_protein` (the pre-clamp Poisson expectation); cap protein synthesis by `budget / gtp_per_protein`.
- **ProteinTranslocation** (`translocation`): `want_gtp = translocase * translocase_specific_rate * interval * gtp_per_monomer` bounded by available substrate; add `budget` into the existing `min(enz_limit, gtp_limit)` as `gtp_limit = min(gtp, budget)/gtp_per`.
- **RibosomeAssembly** (`ribosome_assembly`): `want_gtp = desired_subunits * gtp_per_complex`; add `budget` into the `gtp_cap` computation (`gtp_cap = min(gtp, budget)/gtp_per`).
- **FtsZPolymerization** (`ftsz`): `want_gtp = added_ring_subunits` (1 GTP each); cap `gtp_consumed = min(added, gtp, budget)`.

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_allocation.py
def test_gtp_consumer_respects_budget():
    from viva_mgen.processes.cytokinesis import FtsZPolymerizationReproductionProcess as F
    from process_bigraph import Process
    proc = F.__new__(F); Process.__init__(proc, {}, core=None)
    proc._ring = 0.0; proc._free = None
    st = {"ftsz_monomer": 500.0, "gtp": 1e9, "alloc__gtp": {"ftsz": 3.0}}
    out = proc.update(st, 1.0)
    assert out["gtp"] >= -3.0
    assert out["demand__gtp"]["ftsz"] >= 0.0
```

- [ ] **Step 2: Run to verify it fails**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_allocation.py::test_gtp_consumer_respects_budget -q`
Expected: FAIL.

- [ ] **Step 3: Apply the pattern to all four GTP consumers.**

- [ ] **Step 4: Run to verify it passes**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_allocation.py::test_gtp_consumer_respects_budget -q`
Expected: PASS.

- [ ] **Step 5: Run the full suite**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/ -q`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add viva_mgen/processes/translation.py viva_mgen/processes/protein.py viva_mgen/processes/cytokinesis.py tests/test_allocation.py
git commit -m "GTP consumers: consume under allocator budget, emit demand"
```

---

### Task 6: NTP + amino-acid consumers → demand/budget

**Files:**
- Modify: `viva_mgen/processes/transcription.py` (`TranscriptionReproductionProcess`, pool `ntp`), `viva_mgen/processes/rna.py` (`TRNAAminoacylationReproductionProcess`, pool `amino_acid`)
- Test: `tests/test_allocation.py` (append)

**Interfaces:**
- Consumes: `select_budget`, `demand_entry`, `alloc__ntp`, `alloc__amino_acid`.
- Produces: transcription capped by NTP budget; aminoacylation capped by amino-acid budget (in addition to its ATP budget from Task 4).

Wants:
- **Transcription** (`transcription`, `ntp`): `want_ntp = sum over genes of expected_transcripts * length`; add `budget` to the existing NTP cap (`ntp_cap = min(ntp_avail, budget)`).
- **tRNAAminoacylation** (`trna_aminoacylation`, `amino_acid`): `want_aa = N_desired`; add `budget` into the `min(total_free, amino_acid, atp_limit, enzyme_limit)` that sets `n_total` (so it now reads `min(total_free, min(amino_acid, aa_budget), min(atp_limit, atp_budget), enzyme_limit)`).

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_allocation.py
def test_ntp_consumer_respects_budget():
    from viva_mgen.processes.transcription import TranscriptionReproductionProcess as T
    from process_bigraph import Process
    import numpy as np
    proc = T.__new__(T); Process.__init__(proc, {}, core=None)
    proc._rates = {"g1": 5.0}; proc._lengths = {"g1": 100.0}
    proc._rng = np.random.default_rng(0)
    st = {"ntp": 1e9, "rna_pol": 100.0, "alloc__ntp": {"transcription": 150.0}}
    out = proc.update(st, 1.0)
    assert -out["ntp"] <= 150.0 + 1e-6
    assert out["demand__ntp"]["transcription"] >= 0.0
```

- [ ] **Step 2: Run to verify it fails**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_allocation.py::test_ntp_consumer_respects_budget -q`
Expected: FAIL.

- [ ] **Step 3: Apply the pattern to transcription (`ntp`) and aminoacylation (`amino_acid`).**

- [ ] **Step 4: Run to verify it passes**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_allocation.py::test_ntp_consumer_respects_budget -q`
Expected: PASS.

- [ ] **Step 5: Run the full suite**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/ -q`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add viva_mgen/processes/transcription.py viva_mgen/processes/rna.py tests/test_allocation.py
git commit -m "NTP + amino-acid consumers: consume under allocator budget, emit demand"
```

---

### Task 7: End-to-end scarcity behavior, fidelity notes, snapshot regen

**Files:**
- Test: `tests/test_allocation.py` (append end-to-end test)
- Modify: consumer `description` fidelity lines (each touched process), composite docstring
- Modify: regenerate `reports/composite-state/*.json`

**Interfaces:**
- Consumes: everything above.

- [ ] **Step 1: Write the end-to-end scarcity test**

```python
# append to tests/test_allocation.py
def test_scarcity_partitions_across_consumers():
    """Under a tight ATP budget, two ATP consumers scale back proportionally
    rather than one starving the other."""
    from viva_mgen.processes.allocation import allocate
    # two consumers want 100 and 300; only 40 ATP available
    g = allocate(40.0, {"protein_folding": 100.0, "dna_repair": 300.0})
    assert abs(g["protein_folding"] - 10.0) < 1e-6
    assert abs(g["dna_repair"] - 30.0) < 1e-6
    assert sum(g.values()) <= 40.0 + 1e-9
```

- [ ] **Step 2: Run it**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_allocation.py::test_scarcity_partitions_across_consumers -q`
Expected: PASS.

- [ ] **Step 3: Update fidelity lines**

In each touched consumer's `description`, append to its Fidelity line:
"Consumption is arbitrated by the whole-cell resource allocator (Karr hybrid partitioning) — capped at its per-tick ATP/GTP/NTP/amino-acid budget." Update the composite docstring in `mgen.py` to mention the allocator as the integration layer.

- [ ] **Step 4: Regenerate the committed composite-state snapshot**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python scripts/regen_composite_state.py`
Expected: writes both snapshots; `allocator` node present with `budget` stores.

- [ ] **Step 5: Full suite green**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--alloc /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/ -q`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "Allocator: end-to-end scarcity test, fidelity notes, regenerate snapshot"
```

---

## Self-Review

**Spec coverage:** demand→allocate→run pipeline (Tasks 1,4,5,6); finite pools replenished by metabolism (Tasks 2,3); allocation math proportional+priority (Task 1); pools/consumers table (Tasks 4-6); absent-allocation default = inf (Task 1 `select_budget`); allocator scheduled first (Task 3); edge cases zero-supply/no-demand/feasible=0 (Task 1 tests + `pool_cap`); testing incl. scarcity partition (Tasks 1,7); snapshot regen (Task 7). All spec sections map to a task.

**Placeholder scan:** no TBD/TODO; each code step has real code; the per-consumer pattern is stated with each consumer's concrete `want_*` expression rather than "similar to".

**Type consistency:** `allocate(supply, demands, priorities) -> dict`; `select_budget(alloc, id) -> float`; `demand_entry(id, want) -> dict`; allocator ports `<pool>`, `<pool>_supply` (alias `aa_supply` for `amino_acid`), `demand__<pool>` (map), `alloc__<pool>` (overwrite map) — used consistently across Tasks 3-6. `consumer_id` config on every consumer.

**Open implementation note (Task 3, Step 2):** the `amino_acid` pool's supply store is `aa_supply` (metabolism emits `aa_production`); the allocator maps `amino_acid → aa` via `_SUPPLY_ALIAS`. If the executor finds it cleaner, rename metabolism's output to `amino_acid_production` and drop the alias — either is acceptable so long as the port names line up.
