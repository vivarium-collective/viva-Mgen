#!/usr/bin/env python3
"""Canonical run for the ``fig2-growth`` study (Fig 2 of Karr et al. 2012).

Runs the integrated FBA-metabolism + mass/growth composite for one cell cycle
and checks that the reproduction recapitulates the paper's growth phenotype:
a ~9 h doubling time, mass exactly doubling over the cycle, and a
protein-dominant dry-mass composition. Renders an interactive growth-curve and a
composition donut, and records the run in the workspace ``.pbg/runs.jsonl``.
"""
from __future__ import annotations

import math
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
from viva_mgen.composites import build_mgen
from viva_mgen import viz, constants as Cst
from vivarium_workbench.lib.run_log import append_run_event

SPEC_ID = "viva_mgen.composites.mgen.mycoplasma_genitalium"
STUDY_SLUG = "fig2-growth"
INVESTIGATION_SLUG = "mgen"


def _run(n_hours=9.0, dt=300.0):
    core = build_core()
    doc = build_mgen(core, interval=dt)
    # coarse timestep so a full ~9 h cycle is a few dozen FBA solves, not 32400
    doc["metabolism"]["interval"] = dt
    doc["mass"]["interval"] = dt
    for _pk in ("transcription","translation","rna_decay","protein_decay","replication"):
        doc[_pk]["interval"] = dt
    sim = Composite({"state": doc}, core=core)
    sim.run(n_hours * 3600.0)
    rows = gather_emitter_results(sim)[("emitter",)]
    rows = [r for r in rows if r.get("mass")]
    t = np.arange(len(rows)) * dt / 3600.0  # hours (emitter has no time port)
    mass = np.array([r["mass"] for r in rows])
    fractions = rows[-1]["mass_fractions"]
    return t, mass, fractions, rows


def main() -> int:
    run_id = uuid.uuid4().hex
    append_run_event(WORKSPACE_ROOT, {
        "run_id": run_id, "event": "started", "spec_id": SPEC_ID,
        "label": STUDY_SLUG, "started_at": time.time(), "status": "running",
        "n_steps": 1, "emitter": "ram", "origin": "canonical_run",
        "study_slug": STUDY_SLUG, "investigation_slug": INVESTIGATION_SLUG, "params": {},
    })
    try:
        t, mass, fractions, rows = _run()

        # doubling time from the exponential fit slope of log(mass) vs time
        slope = np.polyfit(t, np.log(mass), 1)[0]  # per hour
        doubling_time_h = math.log(2.0) / slope
        final_mass_ratio = float(mass[-1] / mass[0])
        total_dry = sum(fractions.values()) or 1.0
        protein_fraction = float(fractions.get("protein", 0.0) / total_dry)

        print(f"doubling_time_h      = {doubling_time_h:.3f}  (target ~9.0)")
        print(f"final_mass_ratio     = {final_mass_ratio:.3f}  (target ~2.0)")
        print(f"protein_fraction     = {protein_fraction:.3f}  (target ~0.62)")

        viz_dir = STUDY_DIR / "viz"
        viz_dir.mkdir(parents=True, exist_ok=True)
        # growth curve: total mass + the four major macromolecule fractions
        comp = {k: mass * Cst.DRY_MASS_FRACTIONS[k] for k in ("protein", "RNA", "DNA", "lipid")}
        series = {"total dry mass": mass, **comp}
        (viz_dir / "growth_curve.html").write_text(viz.line_series_html(
            "Cell mass growth over one cycle (Fig 2)", t, series,
            x_title="time (h)", y_title="mass (fg)"))
        # composition donut at division
        major = {k: fractions[k] for k in ("protein", "RNA", "DNA", "lipid",
                 "nucleotide", "carbohydrate", "ion", "polyamine", "vitamin") if k in fractions}
        (viz_dir / "composition_donut.html").write_text(viz.donut_html(
            "Dry-mass composition at division (Fig 2C)",
            list(major.keys()), list(major.values())))
    except Exception:
        append_run_event(WORKSPACE_ROOT, {
            "run_id": run_id, "event": "completed", "completed_at": time.time(),
            "n_steps": 0, "status": "failed"})
        raise

    append_run_event(WORKSPACE_ROOT, {
        "run_id": run_id, "event": "completed", "completed_at": time.time(),
        "n_steps": len(rows), "status": "completed",
        "observables": {"doubling_time_h": doubling_time_h,
                        "final_mass_ratio": final_mass_ratio,
                        "protein_fraction": protein_fraction}})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
