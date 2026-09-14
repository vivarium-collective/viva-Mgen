# Single-cell Ensembles Implementation Plan (gap #7)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans. Steps use checkbox (`- [ ]`).

**Goal:** Add a reusable ensemble layer that runs N independent seeded whole-cell composites and aggregates observables (mean ± spread + end-of-cycle distribution), the infrastructure for reproducing Karr's cell-to-cell variation. No model change.

**Architecture:** `viva_mgen/ensemble.py` runs `build_mgen(seed=…)` composites per seed, gathers emitter rows, and aggregates a scalar observable across cells; `scripts/run_ensemble.py` is a CLI demonstrating it.

**Tech Stack:** Python 3.12, process-bigraph (`Composite`, `gather_emitter_results`), numpy, pytest. Run with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--ensembles` and venv `/Users/eranagmon/code/viva-mGen/.venv`.

**Spec:** `docs/superpowers/specs/2026-09-14-ensembles-design.md`

## Global Constraints

- Work only in worktree `/Users/eranagmon/code/viva-mGen--ensembles` (branch `feat/ensembles`).
- Run everything with `PYTHONPATH=/Users/eranagmon/code/viva-mGen--ensembles /Users/eranagmon/code/viva-mGen/.venv/bin/python`.
- No AI-attribution trailers.
- No change to the single-cell composite or any process — ensembles only RUN it.
- Offline; keep test ensembles tiny (n≤3, short duration) so the suite stays fast.

## Facts: `build_mgen(core=None, *, ..., seed=0, interval=1.0, ...)` seeds transcription/translation/replication/chromosome. Results via `gather_emitter_results(sim)[("emitter",)]` → list of per-emit-step dicts. Build the core with `viva_mgen.core.build_core`.

---

### Task 1: ensemble module + tests

**Files:**
- Create: `viva_mgen/ensemble.py`
- Test: `tests/test_ensemble.py`

**Interfaces:**
- Produces: `run_ensemble(n_cells, duration, *, build_kwargs=None, seeds=None, core=None) -> dict`; `aggregate(cells, observable) -> dict`; `final_values(cells, observable) -> list`.

- [ ] **Step 1: Write failing tests**

```python
# tests/test_ensemble.py
import pytest
from viva_mgen.ensemble import run_ensemble, aggregate, final_values

def test_run_ensemble_independent_cells():
    ens = run_ensemble(n_cells=3, duration=5.0)
    assert len(ens["cells"]) == 3 and len(ens["seeds"]) == 3
    assert len(set(ens["seeds"])) == 3
    assert all(len(rows) > 0 for rows in ens["cells"])

def test_cells_diverge():
    ens = run_ensemble(n_cells=3, duration=10.0)
    # total mRNA at the last common step differs across at least two cells (stochastic)
    def total_rna(rows):
        rc = rows[-1].get("rna_counts", {})
        return sum(rc.values()) if isinstance(rc, dict) else 0.0
    totals = [total_rna(c) for c in ens["cells"]]
    assert len(set(totals)) > 1   # not all identical -> seeds genuinely diverge

def test_aggregate_shapes():
    ens = run_ensemble(n_cells=3, duration=10.0)
    agg = aggregate(ens["cells"], "mass")
    n = len(agg["mean"])
    assert n > 0 and len(agg["std"]) == n and len(agg["per_cell"]) == 3
    assert all(s >= 0.0 for s in agg["std"])

def test_final_values():
    ens = run_ensemble(n_cells=3, duration=5.0)
    fv = final_values(ens["cells"], "mass")
    assert len(fv) == 3

def test_missing_observable_raises():
    ens = run_ensemble(n_cells=2, duration=5.0)
    with pytest.raises((KeyError, ValueError)):
        aggregate(ens["cells"], "no_such_observable")
```

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--ensembles /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_ensemble.py -q`
Expected: FAIL — module missing.

- [ ] **Step 3: Implement `viva_mgen/ensemble.py`**

```python
"""Single-cell ensembles: run N independent seeded whole-cell composites and
aggregate observables to mean +/- spread + end-of-cycle distributions — the
infrastructure for reproducing Karr 2012's cell-to-cell variation. The submodels
are already stochastic; this exercises them as a population (no model change).
See docs/superpowers/specs/2026-09-14-ensembles-design.md.
"""
from __future__ import annotations

import numpy as np
from process_bigraph import Composite, gather_emitter_results

from .composites.mgen import build_mgen
from .core import build_core


def run_ensemble(n_cells, duration, *, build_kwargs=None, seeds=None, core=None):
    """Run ``n_cells`` independent composites (one per seed) for ``duration`` s,
    returning ``{"seeds": [...], "cells": [rows_per_cell]}`` where each entry is
    the cell's emitter rows (list of per-step observable dicts)."""
    if core is None:
        core = build_core()
    if seeds is None:
        seeds = list(range(int(n_cells)))
    seeds = [int(s) for s in seeds]
    base = dict(build_kwargs or {})
    cells = []
    for s in seeds:
        kw = dict(base)
        kw["seed"] = s
        doc = build_mgen(core=core, **kw)
        sim = Composite({"state": doc}, core=core)
        sim.run(float(duration))
        rows = gather_emitter_results(sim).get(("emitter",), [])
        cells.append(rows)
    return {"seeds": seeds, "cells": cells}


def _scalar(row, observable):
    v = row.get(observable)
    if isinstance(v, (int, float)):
        return float(v)
    return None


def aggregate(cells, observable):
    """Mean/std across cells at each timepoint for a SCALAR ``observable``.
    Aligns cells to the common (minimum) row count. Raises if no cell has it."""
    series = []
    for rows in cells:
        vals = [_scalar(r, observable) for r in rows]
        if any(v is not None for v in vals):
            series.append([v if v is not None else 0.0 for v in vals])
    if not series:
        raise ValueError(f"observable {observable!r} not found (or non-scalar) in any cell")
    n_steps = min(len(s) for s in series)
    if n_steps == 0:
        raise ValueError(f"no emitter rows for observable {observable!r}")
    mat = np.array([s[:n_steps] for s in series], dtype=float)  # (cells, steps)
    return {
        "observable": observable,
        "n": len(series),
        "t_index": list(range(n_steps)),
        "mean": mat.mean(axis=0).tolist(),
        "std": mat.std(axis=0).tolist(),
        "per_cell": mat.tolist(),
    }


def final_values(cells, observable):
    """Last value of a scalar ``observable`` per cell (end-of-cycle distribution)."""
    out = []
    for rows in cells:
        val = None
        for r in reversed(rows):
            v = _scalar(r, observable)
            if v is not None:
                val = v
                break
        out.append(val)
    return out
```

- [ ] **Step 4: Run to verify pass**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--ensembles /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/test_ensemble.py -q`
Expected: PASS (5 tests). If `test_cells_diverge` is flaky (cells identical), lengthen `duration` to 30.0 so stochastic events accumulate — but with distinct seeds they should differ quickly.

- [ ] **Step 5: Commit**

```bash
git add viva_mgen/ensemble.py tests/test_ensemble.py
git commit -m "Add single-cell ensemble runner + aggregation (cell-to-cell variation infra)"
```

---

### Task 2: demonstration CLI + backlog update

**Files:**
- Create: `scripts/run_ensemble.py`
- Modify: `docs/FIDELITY_GAPS.md`
- Test: `tests/test_ensemble.py` (append a CLI smoke test)

**Interfaces:**
- Consumes: `viva_mgen.ensemble` (Task 1).

- [ ] **Step 1: Implement `scripts/run_ensemble.py`**

```python
#!/usr/bin/env python
"""Run an ensemble of the whole-cell composite and write an aggregated observable
(mean +/- spread + end-of-cycle distribution) — a demonstration of viva_mgen's
single-cell variation infrastructure.

  python scripts/run_ensemble.py --n 8 --hours 1.0 --observable mass \
      --out reports/ensemble_mass.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from viva_mgen.ensemble import run_ensemble, aggregate, final_values


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--hours", type=float, default=1.0)
    ap.add_argument("--observable", default="mass")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    ens = run_ensemble(args.n, args.hours * 3600.0)
    agg = aggregate(ens["cells"], args.observable)
    fv = [v for v in final_values(ens["cells"], args.observable) if v is not None]
    payload = {**agg, "seeds": ens["seeds"], "final_values": fv,
               "hours": args.hours}
    out = Path(args.out) if args.out else Path(f"reports/ensemble_{args.observable}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2))
    if fv:
        print(f"{args.observable}: {len(fv)} cells, final mean={np.mean(fv):.4g} "
              f"std={np.std(fv):.4g} (CV={np.std(fv)/np.mean(fv):.2%})" if np.mean(fv) else
              f"{args.observable}: {len(fv)} cells, final mean=0")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: CLI smoke test (append to tests/test_ensemble.py)**

```python
def test_cli_writes_aggregate(tmp_path):
    import subprocess, sys, json, os
    out = tmp_path / "ens.json"
    env = dict(os.environ, PYTHONPATH="/Users/eranagmon/code/viva-mGen--ensembles")
    r = subprocess.run(
        [sys.executable, "scripts/run_ensemble.py", "--n", "2", "--hours", "0.01",
         "--observable", "mass", "--out", str(out)],
        cwd="/Users/eranagmon/code/viva-mGen--ensembles", env=env,
        capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr
    d = json.loads(out.read_text())
    assert d["observable"] == "mass" and d["n"] == 2 and len(d["final_values"]) == 2
```
(Use the venv python — the test uses `sys.executable`; run the suite with the venv so `sys.executable` is the venv interpreter.)

- [ ] **Step 3: Backlog update**

In `docs/FIDELITY_GAPS.md`, mark Gap 7 DONE: ensemble runner + aggregation + demo CLI landed (independent seeded cells → mean±spread + end-of-cycle distributions); a dashboard variation-figure study consuming it is staged.

- [ ] **Step 4: Run full suite**

Run: `PYTHONPATH=/Users/eranagmon/code/viva-mGen--ensembles /Users/eranagmon/code/viva-mGen/.venv/bin/python -m pytest tests/ -q`
Expected: PASS (all existing + new). (No snapshot regen needed — no composite/process change.)

- [ ] **Step 5: Commit**

```bash
git add scripts/run_ensemble.py tests/test_ensemble.py docs/FIDELITY_GAPS.md
git commit -m "Ensemble demo CLI + mark gap #7 done in backlog"
```

---

## Self-Review

**Spec coverage:** run_ensemble/aggregate/final_values (T1); CLI (T2 S1); tests incl. divergence + missing-observable + CLI smoke (T1+T2); backlog (T2 S3). Covered.

**Placeholder scan:** all code present; no TBD.

**Type consistency:** `run_ensemble(...) -> {"seeds","cells"}`; `aggregate(cells, observable) -> {observable,n,t_index,mean,std,per_cell}`; `final_values(cells, observable) -> list`. Consistent across module, CLI, tests.
