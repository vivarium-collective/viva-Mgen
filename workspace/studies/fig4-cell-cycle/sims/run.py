#!/usr/bin/env python3
"""Canonical run for the ``fig4-cell-cycle`` study (Fig 4 of Karr et al. 2012).

Reproduces the *emergent*, genetically-unregulated control of cell-cycle
duration that the whole-cell model uncovered. A population of ~128 single cells
is simulated by stepping the ``ReplicationReproductionProcess`` directly (drawing
a random birth DnaA and dNTP level per cell), and three single-cell correlations
are recovered:

* Fig 4C — more initial DnaA ⇒ shorter replication-*initiation* duration
  (negative ``r`` between ``initial_dnaA`` and ``initiation_duration``).
* Fig 4E — initiation and replication durations are *inversely* related: a longer
  initiation builds a larger dNTP surplus that speeds the subsequent replication,
  buffering total cell-cycle length (negative ``r``).
* Fig 4D — a higher dNTP pool at the start of replication ⇒ shorter replication
  duration (negative ``r``).

A single representative cell is also run through the full ``fig4_cell_cycle``
composite to emit a dynamics trajectory (replicated fraction, dNTP pool, DnaA
complex vs time). Reduced mechanism, but the dynamics are genuine: a
dNTP-buffered three-phase (initiation -> replication -> done) replication cycle.
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
from viva_mgen.processes.replication import ReplicationReproductionProcess
from viva_mgen.composites import build_mgen
from viva_mgen import viz
from vivarium_workbench.lib.run_log import append_run_event

SPEC_ID = "viva_mgen.composites.mgen.mycoplasma_genitalium"
STUDY_SLUG = "fig4-cell-cycle"
INVESTIGATION_SLUG = "mgen"

N_CELLS = 128
DT = 100.0
MAX_STEPS = 500
_DONE = 2.0


def _run_population(core):
    """Step the replication process directly for a population of single cells."""
    init_dur, repl_dur, dntp_start, dnaA0 = [], [], [], []
    for i in range(N_CELLS):
        rng = np.random.default_rng(1000 + i)
        initial_dnaA = float(rng.uniform(0.0, 15.0))
        initial_dntp = float(rng.uniform(0.0, 80000.0))
        p = ReplicationReproductionProcess(
            {"initial_dnaA": initial_dnaA, "initial_dntp": initial_dntp, "seed": int(i)},
            core=core,
        )
        out = None
        for _ in range(MAX_STEPS):
            out = p.update({"dntp_synthesis_scale": 1.0}, DT)
            if out["phase_code"] == _DONE:
                break
        if out is None or out["phase_code"] != _DONE:
            continue  # never finished — skip (should not happen for this range)
        init_dur.append(out["initiation_duration"])
        repl_dur.append(out["replication_duration"])
        dntp_start.append(out["dntp_at_replication_start"])
        dnaA0.append(initial_dnaA)
    return (np.array(init_dur), np.array(repl_dur),
            np.array(dntp_start), np.array(dnaA0))


def _run_trajectory(core, initial_dnaA=8.0, initial_dntp=20000.0, seed=1,
                    interval=100.0, n_seconds=40000.0):
    """One representative cell through the full composite for a dynamics viz."""
    doc = build_mgen(core, initial_dnaA=initial_dnaA, initial_dntp=initial_dntp, seed=seed)
    doc["replication"]["interval"] = interval
    for _pk in ("metabolism","mass","transcription","translation","rna_decay","protein_decay"):
        doc[_pk]["interval"] = n_seconds  # fig4 measures the cell cycle only
    sim = Composite({"state": doc}, core=core)
    sim.run(n_seconds)
    rows = gather_emitter_results(sim)[("emitter",)]
    rows = [r for r in rows if r is not None and "replicated_fraction" in r]
    t = np.arange(len(rows)) * interval / 3600.0  # hours
    frac = np.array([r["replicated_fraction"] for r in rows])
    dntp = np.array([r["dntp_pool"] for r in rows])
    dnaA = np.array([r["dnaA_complex"] for r in rows])
    return t, frac, dntp, dnaA


def main() -> int:
    run_id = uuid.uuid4().hex
    append_run_event(WORKSPACE_ROOT, {
        "run_id": run_id, "event": "started", "spec_id": SPEC_ID,
        "label": STUDY_SLUG, "started_at": time.time(), "status": "running",
        "n_steps": 1, "emitter": "ram", "origin": "canonical_run",
        "study_slug": STUDY_SLUG, "investigation_slug": INVESTIGATION_SLUG, "params": {},
    })
    try:
        core = build_core()
        init_dur, repl_dur, dntp_start, dnaA0 = _run_population(core)
        n_cells_completed = int(len(init_dur))
        if n_cells_completed < 100:
            raise RuntimeError(
                f"only {n_cells_completed} cells completed (<100); raise initial_dnaA range")

        r_init_repl = float(np.corrcoef(init_dur, repl_dur)[0, 1])
        r_dntp_repl = float(np.corrcoef(dntp_start, repl_dur)[0, 1])
        r_dnaA_init = float(np.corrcoef(dnaA0, init_dur)[0, 1])

        print(f"n_cells_completed = {n_cells_completed}")
        print(f"r_init_repl = {r_init_repl:.3f}  (Fig 4E, expect NEGATIVE)")
        print(f"r_dntp_repl = {r_dntp_repl:.3f}  (Fig 4D, expect NEGATIVE)")
        print(f"r_dnaA_init = {r_dnaA_init:.3f}  (Fig 4C, expect NEGATIVE)")

        # single representative cell dynamics trajectory
        t, frac, dntp, dnaA = _run_trajectory(core)

        viz_dir = STUDY_DIR / "viz"
        viz_dir.mkdir(parents=True, exist_ok=True)

        (viz_dir / "init_vs_repl.html").write_text(viz.scatter_html(
            "Inverse initiation↔replication duration control (Fig 4E)",
            init_dur / 3600.0, repl_dur / 3600.0,
            x_title="replication-initiation duration (h)",
            y_title="replication duration (h)"))

        (viz_dir / "dntp_vs_repl.html").write_text(viz.scatter_html(
            "dNTP surplus controls replication duration (Fig 4D)",
            dntp_start, repl_dur / 3600.0,
            x_title="dNTP pool at replication start (molecules)",
            y_title="replication duration (h)"))

        dntp_norm = dntp / (dntp.max() or 1.0)
        dnaA_norm = dnaA / (dnaA.max() or 1.0)
        (viz_dir / "cell_cycle_trajectory.html").write_text(viz.line_series_html(
            "Single-cell cell-cycle trajectory (Fig 4B)", t,
            {"replicated fraction": frac,
             "dNTP pool (normalized)": dntp_norm,
             "DnaA complex (normalized)": dnaA_norm},
            x_title="time (h)", y_title="fraction / normalized level"))
    except Exception:
        append_run_event(WORKSPACE_ROOT, {
            "run_id": run_id, "event": "completed", "completed_at": time.time(),
            "n_steps": 0, "status": "failed"})
        raise

    observables = {
        "r_init_repl": r_init_repl,
        "r_dntp_repl": r_dntp_repl,
        "r_dnaA_init": r_dnaA_init,
        "n_cells_completed": n_cells_completed,
    }
    append_run_event(WORKSPACE_ROOT, {
        "run_id": run_id, "event": "completed", "completed_at": time.time(),
        "n_steps": n_cells_completed, "status": "completed",
        "observables": observables})
    print("observables:", observables)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
