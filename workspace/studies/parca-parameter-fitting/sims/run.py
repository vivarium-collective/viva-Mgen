#!/usr/bin/env python3
"""Canonical run for the ``parca-parameter-fitting`` study.

Runs the native parameter calculator (ParCa) — the reproduction of Karr 2012
``FitConstants`` — which the seven figure studies depend on. It computes the
per-gene expression/decay/synthesis panel from the real observed knowledge-base
data, fits it under the RNA-mass / DnaA-FtsZ-held / net-supercoiling constraints,
and closes the expression↔metabolism loop. Records the fit summary + feasibility
as report-card observables and renders the parameter figures.
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

from process_bigraph import Composite

from viva_mgen.core import build_core
from viva_mgen.composites.mgen import build_parca
from vivarium_workbench.lib.run_log import append_run_event

SPEC_ID = "viva_mgen.composites.mgen.mycoplasma_parca"
STUDY_SLUG = "parca-parameter-fitting"
INVESTIGATION_SLUG = "mgen"


def _run() -> dict:
    core = build_core()
    comp = Composite({"state": build_parca(core=core)}, core=core)
    comp.run(1)
    return dict(comp.state["parca"])


def main() -> int:
    run_id = uuid.uuid4().hex
    append_run_event(WORKSPACE_ROOT, {
        "run_id": run_id, "event": "started", "spec_id": SPEC_ID,
        "label": STUDY_SLUG, "started_at": time.time(), "status": "running",
        "n_steps": 1, "emitter": "ram", "origin": "canonical_run",
        "study_slug": STUDY_SLUG, "investigation_slug": INVESTIGATION_SLUG, "params": {},
    })
    try:
        s = _run()
        obs = {
            "n_genes": float(s["n_genes"]),
            "mrna_halflife_min_mean": float(s["mrna_halflife_min_mean"]),
            "rrna_halflife_min": float(s["rrna_halflife_min"]),
            "median_mrna_synthesis_rate": float(s["median_mrna_synthesis_rate"]),
            "nmp_au_fraction": float(s["nmp_au_fraction"]),
            "closed_loop_feasible": float(s["closed_loop_feasible"]),
            "nmp_supply_flux": float(s["nmp_supply_flux"]),
            "aa_supply_flux": float(s["aa_supply_flux"]),
            "growth_baseline": float(s["growth_baseline"]),
        }
        for k, v in obs.items():
            print(f"{k:28s} = {v:.6g}")

        import scripts.regen_parca as regen  # noqa: E402  (renders the figures)
        regen.main()
    except Exception:
        append_run_event(WORKSPACE_ROOT, {
            "run_id": run_id, "event": "completed", "completed_at": time.time(),
            "n_steps": 0, "status": "failed"})
        raise

    append_run_event(WORKSPACE_ROOT, {
        "run_id": run_id, "event": "completed", "completed_at": time.time(),
        "n_steps": 1, "status": "completed", "observables": obs})
    return 0


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(WORKSPACE_ROOT))
    raise SystemExit(main())
