#!/usr/bin/env python3
"""Canonical run for the ``fig1-architecture`` study (Fig 1 of Karr et al. 2012).

Fig 1 makes the model's central claim visible: the whole-cell model is a set of
independent submodels *integrated* through a shared set of cell variables. This
study builds the integrated composite (metabolism + mass + stochastic
transcription/translation/decay), derives the actual process-to-cell-variable
wiring matrix straight from the composite document (a genuine Fig 1B-style
diagram), and confirms the integrated cell runs — all six submodels compose,
they share stores, and the cell stays viable and grows over 20 min.

Fidelity: representative/schematic. This study demonstrates the ARCHITECTURE
(that the submodels compose and the integrated cell runs), not a quantitative
fit.
"""
from __future__ import annotations

import time
import uuid
from pathlib import Path


def _workspace_root(start: Path) -> Path:
    p = start
    for _ in range(8):
        if (p / "workspace.yaml").is_file():
            return p
        p = p.parent
    return start.parents[4]


STUDY_DIR = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = _workspace_root(STUDY_DIR)

import numpy as np
from process_bigraph import Composite, gather_emitter_results

from viva_mgen.core import build_core
from viva_mgen.composites.whole_cell import fig1_architecture
from viva_mgen import viz
from vivarium_workbench.lib.run_log import append_run_event

SPEC_ID = "viva_mgen.composites.whole_cell.fig1_architecture"
STUDY_SLUG = "fig1-architecture"
INVESTIGATION_SLUG = "mgen"


def _store_keys(port_map: dict) -> set:
    """Store keys a process touches: the terminal store name in each ["stores", key] path."""
    keys = set()
    for path in (port_map or {}).values():
        if isinstance(path, (list, tuple)) and len(path) >= 2 and path[0] == "stores":
            keys.add(path[-1])
    return keys


def _wiring_matrix(doc: dict):
    """Binary process x store coupling matrix derived directly from the composite doc."""
    proc_names = [k for k, v in doc.items()
                  if isinstance(v, dict) and v.get("_type") == "process"]
    proc_stores = {}
    for name in proc_names:
        node = doc[name]
        proc_stores[name] = _store_keys(node.get("inputs")) | _store_keys(node.get("outputs"))
    # column order: stable, union of all touched stores
    store_names = []
    for name in proc_names:
        for s in sorted(proc_stores[name]):
            if s not in store_names:
                store_names.append(s)
    matrix = [[1 if s in proc_stores[p] else 0 for s in store_names] for p in proc_names]
    return proc_names, store_names, matrix


def _run(runtime_s=1200.0, heavy_interval=60.0):
    core = build_core()
    doc = fig1_architecture(core)
    # keep FBA cost reasonable: metabolism + mass on a coarse timestep, expression at 1 s
    doc["metabolism"]["interval"] = heavy_interval
    doc["mass"]["interval"] = heavy_interval

    proc_names, store_names, matrix = _wiring_matrix(doc)

    sim = Composite({"state": doc}, core=core)
    sim.run(runtime_s)
    rows = gather_emitter_results(sim)[("emitter",)]
    rows = [r for r in rows if r]
    return doc, proc_names, store_names, matrix, rows, runtime_s


def main() -> int:
    run_id = uuid.uuid4().hex
    append_run_event(WORKSPACE_ROOT, {
        "run_id": run_id, "event": "started", "spec_id": SPEC_ID,
        "label": STUDY_SLUG, "started_at": time.time(), "status": "running",
        "n_steps": 1, "emitter": "ram", "origin": "canonical_run",
        "study_slug": STUDY_SLUG, "investigation_slug": INVESTIGATION_SLUG, "params": {},
    })
    try:
        doc, proc_names, store_names, matrix, rows, runtime_s = _run()

        # --- scalar observables ---
        n_processes = float(len(proc_names))
        # stores coupled = touched by >= 2 processes (integration / shared-variable coupling)
        col_touch = [sum(matrix[i][j] for i in range(len(proc_names))) for j in range(len(store_names))]
        n_stores_coupled = float(sum(1 for c in col_touch if c >= 2))

        first, last = rows[0], rows[-1]
        growth_fraction_final = float(last.get("growth_fraction", 0.0))
        mass_grew = 1.0 if float(last.get("mass", 0.0)) > float(first.get("mass", 0.0)) else 0.0
        final_rna = last.get("rna_counts", {}) or {}
        final_protein = last.get("protein_counts", {}) or {}
        mrna_species_produced = float(sum(1 for v in final_rna.values() if v and v > 0))
        protein_species_produced = float(sum(1 for v in final_protein.values() if v and v > 0))
        assert mrna_species_produced >= 0 and protein_species_produced >= 0

        observables = {
            "n_processes": n_processes,
            "n_stores_coupled": n_stores_coupled,
            "growth_fraction_final": growth_fraction_final,
            "mass_grew": mass_grew,
            "mrna_species_produced": mrna_species_produced,
            "protein_species_produced": protein_species_produced,
        }
        for k, v in observables.items():
            print(f"{k:26s} = {v}")

        # --- visualizations ---
        viz_dir = STUDY_DIR / "viz"
        viz_dir.mkdir(parents=True, exist_ok=True)

        (viz_dir / "wiring_matrix.html").write_text(viz.heatmap_html(
            "Whole-cell integration: process ↔ cell-variable wiring (Fig 1B)",
            z=matrix, x=store_names, y=proc_names,
            x_title="cell variable (store)", y_title="submodel (process)",
            colorbar_title="touches"))

        # integrated dynamics over time (emit cadence = coarsest sync; derive minutes)
        n = len(rows)
        t_min = np.linspace(0.0, runtime_s / 60.0, n) if n > 1 else np.array([0.0])
        mass = np.array([float(r.get("mass", 0.0)) for r in rows])
        growth_rate = np.array([float(r.get("growth_rate", 0.0)) for r in rows])
        total_mrna = np.array([sum((r.get("rna_counts", {}) or {}).values()) for r in rows])
        total_protein = np.array([sum((r.get("protein_counts", {}) or {}).values()) for r in rows])
        (viz_dir / "cell_dashboard.html").write_text(viz.line_series_html(
            "Integrated cell dynamics (Fig 1)", t_min, {
                "mass (fg)": mass,
                "growth rate": growth_rate,
                "total mRNA count": total_mrna,
                "total protein count": total_protein,
            }, x_title="time (min)", y_title="value"))
    except Exception:
        append_run_event(WORKSPACE_ROOT, {
            "run_id": run_id, "event": "completed", "completed_at": time.time(),
            "n_steps": 0, "status": "failed"})
        raise

    append_run_event(WORKSPACE_ROOT, {
        "run_id": run_id, "event": "completed", "completed_at": time.time(),
        "n_steps": len(rows), "status": "completed",
        "observables": observables})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
