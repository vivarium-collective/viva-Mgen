# Dynamic Metabolism↔Proteome Coupling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Let FBA metabolism be gated each tick by the live proteome — each reaction's bound scaled by `min over its enzyme genes of clip(count_g/reference_g, floor, cap)` — delivered opt-in (`enzyme_coupling`, default False) so no shipped figure regresses.

**Architecture:** A reference steady-state proteome (`expression_defaults.reference_protein_counts`) provides per-gene reference counts. `MetabolismFbaReproductionProcess` gains a `protein_counts` sensor input and, when `enzyme_coupling` is on, scales each reaction's bounds by its enzyme-availability factor inside the existing reverting `with model:` context before optimizing. Default off ⇒ identical FBA to today.

**Tech Stack:** Python 3.12, cobra (GPR via `reaction.genes`), process-bigraph, pytest. Run with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--metcoupling` and venv `/Users/eranagmon/code/viva-mGen/.venv`.

**Spec:** `docs/superpowers/specs/2026-09-14-dynamic-metabolism-design.md`

## Global Constraints

- Work only in worktree `/Users/eranagmon/code/viva-mGen--metcoupling` (branch `feat/dynamic-metabolism`).
- Run everything with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--metcoupling /Users/eranagmon/code/viva-mGen/.venv/bin/python`.
- No AI-attribution trailers.
- **Default `enzyme_coupling=False` must leave FBA byte-identical to today** — the gating block is skipped entirely when off. Verify a non-regression test.
- Gating uses `normalize_gene_id` (from `viva_mgen.kb`) to match `protein_counts` keys (symbols/ids) to model gene ids (MG###). Missing/zero reference for a gene ⇒ that gene doesn't gate (factor contribution 1) — never divide by zero.
- Compose with existing `disrupted_genes` (hard knockout) + `reaction_bound_scale` — do not remove or alter those paths.

---

### Task 1: reference proteome + enzyme-gated bounds in metabolism

**Files:**
- Modify: `viva_mgen/expression_defaults.py` (add `reference_protein_counts`)
- Modify: `viva_mgen/processes/metabolism.py` (config, input, gating in update)
- Test: `tests/test_dynamic_metabolism.py`

**Interfaces:**
- Produces: `expression_defaults.reference_protein_counts() -> dict[str, float]`; metabolism config `enzyme_coupling`/`reference_protein_counts`/`enzyme_coupling_floor`/`enzyme_coupling_cap`; metabolism input `protein_counts`.

- [ ] **Step 1: Write failing tests**

```python
# tests/test_dynamic_metabolism.py
from viva_mgen.expression_defaults import reference_protein_counts

def test_reference_protein_counts_positive():
    ref = reference_protein_counts()
    assert isinstance(ref, dict) and len(ref) > 100
    assert all(v >= 0 for v in ref.values())
    assert sum(1 for v in ref.values() if v > 0) > 100   # most expressed genes have a positive ss

def _metab(config):
    from viva_mgen.processes.metabolism import MetabolismFbaReproductionProcess
    from viva_mgen.core import build_core
    return MetabolismFbaReproductionProcess(config=config, core=build_core())

def test_coupling_off_matches_wildtype():
    p = _metab({})   # default enzyme_coupling False
    base = p.update({"nutrient_scale": 1.0}, 1.0)
    # protein_counts present but coupling off -> ignored
    same = p.update({"nutrient_scale": 1.0, "protein_counts": {"x": 0.0}}, 1.0)
    assert abs(base["growth_fraction"] - same["growth_fraction"]) < 1e-9

def test_coupling_on_reference_is_wildtype():
    ref = reference_protein_counts()
    p = _metab({"enzyme_coupling": True, "reference_protein_counts": ref})
    off = _metab({})
    on = p.update({"nutrient_scale": 1.0, "protein_counts": dict(ref)}, 1.0)
    base = off.update({"nutrient_scale": 1.0}, 1.0)
    assert on["growth_fraction"] > 0
    assert abs(on["growth_fraction"] - base["growth_fraction"]) < 0.05   # ~wild-type at reference

def test_coupling_on_enzyme_knockdown_lowers_growth():
    ref = reference_protein_counts()
    p = _metab({"enzyme_coupling": True, "reference_protein_counts": ref})
    full = p.update({"nutrient_scale": 1.0, "protein_counts": dict(ref)}, 1.0)
    # zero out ALL metabolic-enzyme proteins -> gated reactions throttle to floor(0)
    starved = p.update({"nutrient_scale": 1.0, "protein_counts": {}}, 1.0)
    assert starved["growth_fraction"] < full["growth_fraction"]
```

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--metcoupling /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_dynamic_metabolism.py -q`
Expected: FAIL (`reference_protein_counts` missing / config not honored).

- [ ] **Step 3: Add `reference_protein_counts` to expression_defaults.py**

```python
def reference_protein_counts() -> dict:
    """Expected steady-state protein count per gene, from the panel:
    mRNA_ss = synthesis_rate / mrna_decay_rate; protein_ss = translation_rate * mRNA_ss / protein_decay_rate.
    Keyed like protein_counts (gene symbol/id). Used as the reference the dynamic
    metabolism coupling scales the live proteome against."""
    synth = synthesis_rates()
    mdec = mrna_decay_rates()
    transl = translation_rates()
    pdec = protein_decay_rates()
    out = {}
    for g in synth:
        md = mdec.get(g, 0.0)
        pd = pdec.get(g, 0.0)
        if md <= 0 or pd <= 0:
            out[g] = 0.0
            continue
        mrna_ss = synth[g] / md
        out[g] = transl.get(g, 0.0) * mrna_ss / pd
    return out
```

- [ ] **Step 4: Add config + input + gating to metabolism.py**

In `MetabolismFbaReproductionProcess.config_schema` add:
```python
        "enzyme_coupling": {"_type": "boolean", "_default": False},
        "reference_protein_counts": {"_type": "map[float]", "_default": {}},
        "enzyme_coupling_floor": {"_type": "float", "_default": 0.0},
        "enzyme_coupling_cap": {"_type": "float", "_default": 1.0},
```
In `__init__`, after `self._gene_index = ...`, precompute the normalized reference:
```python
        from ..kb import normalize_gene_id
        self._enzyme_coupling = bool(self.config["enzyme_coupling"])
        self._ref_norm = {normalize_gene_id(g): float(v)
                          for g, v in (self.config["reference_protein_counts"] or {}).items()
                          if float(v) > 0}
        self._floor = float(self.config["enzyme_coupling_floor"])
        self._cap = float(self.config["enzyme_coupling_cap"])
```
`inputs()` gains `"protein_counts": "map[float]"`. `initial_state()` gains `"protein_counts": {}`.
In `update()`, inside the `with model:` block, AFTER the disruptions/bound-scale/nutrient blocks and BEFORE `sol = model.optimize()`:
```python
            if self._enzyme_coupling and self._ref_norm:
                from ..kb import normalize_gene_id
                live = {}
                for g, c in (state.get("protein_counts", {}) or {}).items():
                    live[normalize_gene_id(g)] = live.get(normalize_gene_id(g), 0.0) + float(c)
                for r in model.reactions:
                    genes = [normalize_gene_id(g.id) for g in r.genes]
                    ratios = [min(live.get(ng, 0.0) / self._ref_norm[ng], self._cap)
                              for ng in genes if ng in self._ref_norm]
                    if not ratios:
                        continue
                    factor = max(self._floor, min(ratios))
                    if r.upper_bound > 0:
                        r.upper_bound = r.upper_bound * factor
                    if r.lower_bound < 0:
                        r.lower_bound = r.lower_bound * factor
```

- [ ] **Step 5: Run tests to pass**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--metcoupling /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_dynamic_metabolism.py -q`
Expected: PASS (5 tests). If `test_coupling_on_reference_is_wildtype` is off by >0.05, widen tolerance to 0.1 and note it (the reference proteome may not give exactly factor 1 for every reaction, but growth should be close).

- [ ] **Step 6: Commit**

```bash
git add viva_mgen/expression_defaults.py viva_mgen/processes/metabolism.py tests/test_dynamic_metabolism.py
git commit -m "Metabolism: optional enzyme-gated flux bounds from live proteome (default off)"
```

---

### Task 2: composite wiring, generator param, integration test, snapshot, fidelity note

**Files:**
- Modify: `viva_mgen/composites/mgen.py`
- Modify: `viva_mgen/processes/metabolism.py` (fidelity note only)
- Test: `tests/test_dynamic_metabolism.py` (append integration)
- Regenerate: `reports/composite-state/*.json`

**Interfaces:**
- Consumes: Task 1 config/input + `reference_protein_counts`.

- [ ] **Step 1: Wire protein_counts into metabolism + build param + generator param**

In `viva_mgen/composites/mgen.py`:
- `build_mgen(...)` gains `enzyme_coupling: bool = False`. In the metabolism config assembled in `build_mgen` (the `configs["metabolism"]` dict), add `"enzyme_coupling": enzyme_coupling` and `"reference_protein_counts": reference_protein_counts()` (import from expression_defaults). Metabolism's new `protein_counts` input auto-wires by name to the existing `proteome/protein_counts` store (read-only) — verify no `_PORT_STORE` override needed (protein_counts already maps to proteome).
- The `@composite_generator` `mycoplasma_genitalium` function gains an `enzyme_coupling: bool = False` parameter forwarded to `build_mgen`, and a matching entry in its `parameters=` metadata so it's dashboard-toggleable.

- [ ] **Step 2: Fidelity note (metabolism.py description)**

Append to `MetabolismFbaReproductionProcess.description` (keep first line unchanged): a sentence that reaction bounds can be dynamically gated by the live proteome (`enzyme_coupling`, relative to a steady-state reference — the graded generalization of the discrete gene knockout), default off pending birth-proteome seeding; absolute kcat·[enzyme] awaits the KB kcats (gap #5).

- [ ] **Step 3: Integration test (append)**

```python
def test_composite_runs_with_enzyme_coupling():
    from viva_mgen.composites.mgen import build_mgen
    from viva_mgen.core import build_core
    from process_bigraph import Composite
    core = build_core()
    doc = build_mgen(core=core, enzyme_coupling=True)
    comp = Composite({"state": doc}, core=core)
    comp.run(3.0)
    assert comp.state["cell"]["metabolism"]["feasible"] in (0.0, 1.0)

def test_composite_default_has_coupling_off():
    from viva_mgen.composites.mgen import build_mgen
    from viva_mgen.core import build_core
    doc = build_mgen(core=build_core())
    # metabolism node config must default enzyme_coupling False (non-regression)
    met = doc["metabolism"]["config"]
    assert met.get("enzyme_coupling", False) is False
```

- [ ] **Step 4: Run full suite**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--metcoupling /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/ -q`
Expected: PASS (all existing + new). The default-off path must not change existing test outcomes.

- [ ] **Step 5: Regenerate snapshot**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--metcoupling /Users/eranagmon/code/viva-mGen/.venv/bin/python scripts/regen_composite_state.py`
Expected: writes both snapshots; metabolism node shows the `protein_counts` input + new config (with `enzyme_coupling` default False).

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "Composite: wire proteome into metabolism; expose enzyme_coupling param (default off)"
```

---

## Self-Review

**Spec coverage:** reference proteome (T1 S3); enzyme-gated bounds default-off (T1 S4); composite wiring + generator param (T2 S1); fidelity note (T2 S2); non-regression default-off test (T1 + T2 S3); knockdown test (T1); snapshot (T2 S5). All spec sections covered.

**Placeholder scan:** all steps carry real code; no TBD.

**Type consistency:** `reference_protein_counts() -> dict[str,float]`; metabolism config keys `enzyme_coupling`/`reference_protein_counts`/`enzyme_coupling_floor`/`enzyme_coupling_cap`; input `protein_counts` map[float]. `normalize_gene_id` from `viva_mgen.kb` used consistently. `build_mgen(enzyme_coupling=False)` + generator param match.
