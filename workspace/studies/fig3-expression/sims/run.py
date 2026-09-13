#!/usr/bin/env python3
"""Canonical run for the ``fig3-expression`` study (Fig 3 of Karr et al. 2012).

Reproduces the chromosome DNA–protein-interaction phenomena of Figure 3 and
records paper-grounded observables for the study's report cards:

* chromosome-exploration kinetics — % of the genome protein-bound at 6 min and
  20 min, and the time for RNA polymerase to bind 90 % of the chromosome (Fig 3B);
* RNA-expression kinetics — t50, the time for 50 % of genes to be expressed (Fig 3C);
* protein–DNA collisions over a cell cycle — total count, the fraction caused by
  RNA polymerase, the fraction displacing SMC, and the collisions-vs-density
  correlation (Fig 3E/3F).

Chromosome dynamics come from :class:`ChromosomeDynamicsReproductionProcess`
(read, not modified); RNA-expression timing from the transcription process on the
representative gene panel.

REDUCED FIDELITY: binned genome, representative protein counts, a 6-gene
expression panel — so absolute collision counts and t50 are representative, not
KB-fitted. The report cards are set to the paper's values so the divergences are
explicit and can be closed later via the model mechanisms.
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

from viva_mgen.core import build_core
from viva_mgen.processes.chromosome import ChromosomeDynamicsReproductionProcess
from viva_mgen.processes.transcription import TranscriptionReproductionProcess
from viva_mgen.expression_defaults import DEFAULT_GENES
from vivarium_workbench.lib.run_log import append_run_event

SPEC_ID = "viva_mgen.composites.mgen.mycoplasma_genitalium"
STUDY_SLUG = "fig3-expression"
INVESTIGATION_SLUG = "mgen"


def _chromosome_metrics(core, dt=30.0, cycle_s=9 * 3600.0):
    """Run the chromosome process over one cell cycle; return Fig-3 observables."""
    ch = ChromosomeDynamicsReproductionProcess(config={"seed": 0}, core=core)
    t_min, explored, rnap_expl, dens, ncoll = [], [], [], [], []
    out = None
    n = int(cycle_s / dt)
    for k in range(n):
        out = ch.update({"rna_polymerase": 120.0, "replication_active": 1.0}, dt)
        t_min.append(k * dt / 60.0)
        explored.append(out["fraction_explored"] * 100.0)
        rnap_expl.append(out["percent_rnap_explored"] * 100.0)
        dens.append(out["dna_binding_density"])
        ncoll.append(out["n_collisions"])
    t_min = np.array(t_min); explored = np.array(explored)
    rnap_expl = np.array(rnap_expl); dens = np.array(dens); ncoll = np.array(ncoll)

    def _at(tmin):
        return float(np.interp(tmin, t_min, explored))

    def _first_time_ge(arr, thr):
        idx = np.where(arr >= thr)[0]
        return float(t_min[idx[0]]) if len(idx) else float("inf")

    collisions = out["collisions"] or {}
    total = float(sum(collisions.values())) or 1.0
    by_rnap = sum(v for k, v in collisions.items() if k.split("||")[0] == "RNA Pol")
    disp_smc = sum(v for k, v in collisions.items() if k.split("||")[-1] == "SMC")
    coll_delta = np.diff(np.concatenate([[0.0], ncoll]))
    # correlate per-step new collisions with the instantaneous binding density
    r_cd = float(np.corrcoef(dens, coll_delta)[0, 1]) if np.std(dens) and np.std(coll_delta) else 0.0

    return {
        "pct_explored_at_6min": _at(6.0),
        "pct_explored_at_20min": _at(20.0),
        "rnap_90pct_time_min": _first_time_ge(rnap_expl, 90.0),
        "n_collisions_per_cycle": float(ncoll[-1]),
        "frac_collisions_by_rnap": float(by_rnap / total),
        "frac_collisions_displacing_smc": float(disp_smc / total),
        "collisions_density_pearson_r": r_cd,
    }


def _rna_expression_t50(core, dt=1.0, max_min=150.0):
    """Time for 50 % (and track 90 %) of the gene panel to be expressed (Fig 3C)."""
    txn = TranscriptionReproductionProcess(config={"seed": 3}, core=core)
    seen = set()
    n_genes = len(DEFAULT_GENES)
    t50 = t90 = float("inf")
    for k in range(int(max_min * 60 / dt)):
        d = txn.update({"ntp": 1e12, "rna_pol": 100.0}, dt)
        for g, nsyn in (d.get("rna_counts") or {}).items():
            if nsyn > 0:
                seen.add(g)
        frac = len(seen) / n_genes
        tmin = k * dt / 60.0
        if frac >= 0.5 and t50 == float("inf"):
            t50 = tmin
        if frac >= 0.9 and t90 == float("inf"):
            t90 = tmin
            break
    return {"rna_t50_min": t50, "rna_t90_min": t90}


def main() -> int:
    run_id = uuid.uuid4().hex
    append_run_event(WORKSPACE_ROOT, {
        "run_id": run_id, "event": "started", "spec_id": SPEC_ID,
        "label": STUDY_SLUG, "started_at": time.time(), "status": "running",
        "n_steps": 1, "emitter": "ram", "origin": "canonical_run",
        "study_slug": STUDY_SLUG, "investigation_slug": INVESTIGATION_SLUG,
        "params": {"seed": 0},
    })
    try:
        core = build_core()
        obs = {}
        obs.update(_chromosome_metrics(core))
        obs.update(_rna_expression_t50(core))
        for k, v in obs.items():
            print(f"{k:34s} = {v}")
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
    raise SystemExit(main())
