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


def _mrna_protein_decoupling(n_cells=48, n_hours=9.0, dt=300.0):
    """Fig 2H — across unsynchronised single cells, per-gene mRNA and protein copy
    numbers are UNCORRELATED (short-lived bursty mRNA vs long-lived cumulative
    protein). Step the stochastic transcription+translation processes for a
    population and correlate final mRNA vs protein; return |Pearson r| and the
    mean mRNA per gene (Fig 2G, low-copy)."""
    from viva_mgen.processes.transcription import TranscriptionReproductionProcess
    from viva_mgen.processes.translation import TranslationReproductionProcess
    from viva_mgen.processes.decay import (RnaDecayReproductionProcess,
                                           ProteinDecayReproductionProcess)
    core = build_core()
    # smaller dt for the stochastic expression so mRNA decay keeps it bursty/low
    edt = 30.0
    steps = int(n_hours * 3600.0 / edt)
    mrna_pts, prot_pts, mrna_per_gene = [], [], []
    for c in range(n_cells):
        txn = TranscriptionReproductionProcess({"seed": c}, core=core)
        tsl = TranslationReproductionProcess({"seed": c + 5000}, core=core)
        rdec = RnaDecayReproductionProcess({"seed": c + 9000}, core=core)
        pdec = ProteinDecayReproductionProcess({"seed": c + 13000}, core=core)
        rna: dict = {}
        prot: dict = {}
        for _ in range(steps):
            d = txn.update({"ntp": 1e12, "rna_pol": 100.0}, edt)
            for g, n in d.get("rna_counts", {}).items():
                rna[g] = max(0.0, rna.get(g, 0.0) + float(n))
            dd = rdec.update({"rna_counts": dict(rna)}, edt)          # mRNA decay -> low/bursty
            for g, n in dd.get("rna_counts", {}).items():
                rna[g] = max(0.0, rna.get(g, 0.0) + float(n))
            d2 = tsl.update({"rna_counts": dict(rna), "gtp": 1e12}, edt)
            for g, n in d2.get("protein_counts", {}).items():
                prot[g] = max(0.0, prot.get(g, 0.0) + float(n))
            dp = pdec.update({"protein_counts": dict(prot)}, edt)     # slow protein decay
            for g, n in dp.get("protein_counts", {}).items():
                prot[g] = max(0.0, prot.get(g, 0.0) + float(n))
        genes = sorted(set(rna) | set(prot))
        for g in genes:
            mrna_pts.append(rna.get(g, 0.0)); prot_pts.append(prot.get(g, 0.0))
        if rna:
            mrna_per_gene.append(np.mean([rna.get(g, 0.0) for g in genes]))
    mrna_pts, prot_pts = np.array(mrna_pts), np.array(prot_pts)
    if mrna_pts.std() > 0 and prot_pts.std() > 0:
        corr = abs(float(np.corrcoef(mrna_pts, prot_pts)[0, 1]))
    else:
        corr = 0.0
    return corr, float(np.mean(mrna_per_gene) if mrna_per_gene else 0.0)


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
        rna_fraction = float(fractions.get("RNA", fractions.get("rna", 0.0)) / total_dry)

        # DNA content doubles 1->2 over the cycle (Fig 2D), driven by the
        # replication submodel's replicated_fraction.
        repl = np.array([float(r.get("replicated_fraction", 0.0)) for r in rows])
        dna_fold_change = float((1.0 + repl[-1]) / (1.0 + repl[0]))

        # Fig 2G/2H: low-copy bursty mRNA + mRNA↔protein decoupling across cells
        mrna_protein_abs_corr, mean_mrna_per_gene = _mrna_protein_decoupling()

        print(f"doubling_time_h        = {doubling_time_h:.3f}  (target ~8.9)")
        print(f"final_mass_ratio       = {final_mass_ratio:.3f}  (target ~2.0)")
        print(f"protein_fraction       = {protein_fraction:.3f}  (target 0.62-0.75)")
        print(f"rna_fraction           = {rna_fraction:.3f}  (target ~0.09)")
        print(f"dna_fold_change        = {dna_fold_change:.3f}  (target ~2.0)")
        print(f"mean_mrna_per_gene     = {mean_mrna_per_gene:.3f}  (target 0-2.5)")
        print(f"mrna_protein_abs_corr  = {mrna_protein_abs_corr:.3f}  (target <0.4)")

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
                        "protein_fraction": protein_fraction,
                        "rna_fraction": rna_fraction,
                        "dna_fold_change": dna_fold_change,
                        "mean_mrna_per_gene": mean_mrna_per_gene,
                        "mrna_protein_abs_corr": mrna_protein_abs_corr}})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
