#!/usr/bin/env python3
"""Persist a Results-tab-readable zarr store per study.

The workbench Results tab reads a run's time-series via
``explorer_data.get_series``, which needs a FLAT xarray/zarr dataset (observables
as top-level data_vars, not nested ``stores/<name>`` groups). This script runs
each study's baseline composite with a flat-named XArrayEmitter over the study's
readout store-paths, writes a committed zarr under
``workspace/studies/<slug>/results/baseline.zarr``, and records a study-linked
run in ``.pbg/runs.jsonl`` (store_path → that zarr). ``build_study_results`` then
renders a per-store preview, so the Results tab is populated in the read-only
snapshot too (both the zarr and the run log are committed for CI).

Re-run any time to refresh: it rewrites the run log and the per-study stores.
"""
from __future__ import annotations

import json
import shutil
import time
import uuid
from pathlib import Path

WS = Path(__file__).resolve().parents[1]

# study -> (baseline params, {flat_name: store-path}, sim_time_s, interval_s)
STUDIES = {
    "fig1-architecture": ({}, {"mass": "stores/mass", "growth_fraction": "stores/growth_fraction",
                               "growth_rate": "stores/growth_rate"}, 32400.0, 300.0),
    "fig2-growth": ({}, {"mass": "stores/mass", "growth_fraction": "stores/growth_fraction",
                         "growth_rate": "stores/growth_rate", "volume": "stores/volume"}, 32400.0, 300.0),
    "fig3-expression": ({"seed": 0}, {"ntp": "stores/ntp", "gtp": "stores/gtp"}, 3600.0, 30.0),
    "fig4-cell-cycle": ({"initial_dnaA": 5, "initial_dntp": 0, "seed": 0},
                        {"replicated_fraction": "stores/replicated_fraction", "dntp_pool": "stores/dntp_pool",
                         "dnaA_complex": "stores/dnaA_complex", "phase_code": "stores/phase_code"}, 32400.0, 200.0),
    "fig5-energy": ({}, {"atp_production": "stores/atp_production", "gtp_production": "stores/gtp_production",
                         "growth_rate": "stores/growth_rate"}, 32400.0, 300.0),
    "fig6-gene-essentiality": ({}, {"growth_fraction": "stores/growth_fraction",
                                    "feasible": "stores/feasible"}, 10800.0, 100.0),
    "fig7-kinetic-parameters": ({}, {"growth_fraction": "stores/growth_fraction",
                                     "growth_rate": "stores/growth_rate"}, 10800.0, 100.0),
}
SPEC_ID = "viva_mgen.composites.mgen.mycoplasma_genitalium"


def _persist_one(core, slug, params, readouts, sim_time, interval):
    from process_bigraph import Composite
    import pbg_emitters
    from vivarium_workbench.lib.emitters import _xarray_emitter_config, _flush_step_emitters
    from viva_emitters.xarray_emitter.view import view_from_emit_paths
    from viva_mgen.composites import build_mgen

    core.register_link("XArrayEmitter", pbg_emitters.XArrayEmitter)
    run_id = f"{slug}-baseline"  # stable id == zarr experiment_id, referenced by study.yaml runs[]
    results_dir = WS / "workspace" / "studies" / slug / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    store = results_dir / "baseline.zarr"
    if store.exists():
        shutil.rmtree(store)

    doc = build_mgen(core, interval=interval, **params)
    flat = {name: ["stores", sp.split("/", 1)[1]] for name, sp in readouts.items()}
    view = view_from_emit_paths(sorted(flat), dtype="<f8")
    cfg = _xarray_emitter_config(str(store), view, None)
    cfg["metadata"] = {"experiment_id": run_id}
    doc["emitter"] = {"_type": "step", "address": "local:XArrayEmitter",
                      "config": {**cfg, "emit": {k: "node" for k in flat}},
                      "inputs": {**flat, "global_time": ["global_time"]}}
    sim = Composite({"state": doc}, core=core)
    sim.run(sim_time)
    _flush_step_emitters(sim)
    n = int(sim_time / interval)
    return run_id, str(store.relative_to(WS)), n


def _write_runs_block(slug, run_id, store_rel, n):
    """Insert/replace a clean ``runs:`` block in the study.yaml that points the
    Simulations DB + Results tab at the committed zarr store (study-linked)."""
    import re
    p = WS / "workspace" / "studies" / slug / "study.yaml"
    t = p.read_text()
    block = (
        "runs:\n"
        f"- name: baseline\n"
        f"  run_id: {run_id}\n"
        f"  status: completed\n"
        f"  emitter: xarray\n"
        f"  n_steps: {n}\n"
        f"  store_path: {store_rel}\n"
        f"  timestamp: {time.time():.1f}\n\n"
    )
    # drop ALL existing top-level runs: blocks (server flush may add one after
    # pipeline_gate / at EOF), then insert one clean block before pipeline_gate
    t = re.sub(r"(?ms)^runs:\n(?:[ \t-].*\n|\n)*?(?=^[A-Za-z_]+:|\Z)", "", t)
    i = t.find("\npipeline_gate:")
    t = t[:i + 1] + block + t[i + 1:]
    p.write_text(t)


def main() -> int:
    from viva_mgen.core import build_core
    from vivarium_workbench.lib.results_views import build_study_results

    # clean stale run bookkeeping so the committed study.yaml runs win the index
    for stale in [WS / ".pbg" / "composite-runs.db", WS / ".pbg" / "runs.jsonl"]:
        stale.unlink(missing_ok=True)
    for rundir in (WS / ".pbg" / "runs").glob("*"):
        shutil.rmtree(rundir, ignore_errors=True)
    for db in (WS / "workspace" / "studies").glob("*/runs.db"):
        db.unlink(missing_ok=True)

    core = build_core()
    for slug, (params, readouts, sim_time, interval) in STUDIES.items():
        run_id, store_rel, n = _persist_one(core, slug, params, readouts, sim_time, interval)
        _write_runs_block(slug, run_id, store_rel, n)
        print(f"{slug}: {run_id} -> {store_rel} ({n} steps)")

    print("\n=== Results preview per study ===")
    ok = True
    for slug in STUDIES:
        res, _ = build_study_results(WS, slug)
        n = len(res.get("stores", []))
        ok &= bool(res.get("present") and n > 0)
        print(f"  {slug}: present={res.get('present')} stores={n} "
              + (",".join(s["path"] for s in res.get("stores", [])[:4])))
    print("ALL RESULTS READABLE:", ok)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
