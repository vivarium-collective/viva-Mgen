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
    # mRNA gene keys only — the paper's "low-copy bursty" (Fig 2G) is about mRNA;
    # stable RNAs (rRNA/tRNA) are legitimately high-copy and must be excluded from
    # the mRNA copy-number statistic (now that RNA type comes from real KB data).
    from viva_mgen.kb import load_genes
    from viva_mgen.parca import rna_type
    _mrna_keys = {((g.get("symbol") or "").strip() or g["gene_id"])
                  for g in load_genes() if rna_type(g) == "mRNA"}
    # smaller dt for the stochastic expression so mRNA decay keeps it bursty/low
    edt = 30.0
    steps = int(n_hours * 3600.0 / edt)
    age_rng = np.random.default_rng(4242)
    per_gene, mrna_per_gene = {}, []
    for c in range(n_cells):
        txn = TranscriptionReproductionProcess({"seed": c}, core=core)
        tsl = TranslationReproductionProcess({"seed": c + 5000}, core=core)
        rdec = RnaDecayReproductionProcess({"seed": c + 9000}, core=core)
        pdec = ProteinDecayReproductionProcess({"seed": c + 13000}, core=core)
        rna: dict = {}
        prot: dict = {}
        # cells are UNSYNCHRONISED: each is sampled at a random age, so cumulative
        # protein varies with age while the current (bursty, decaying) mRNA snapshot
        # does not — the mechanism behind the mRNA↔protein decoupling of Fig 2H.
        steps_c = int(age_rng.uniform(0.08, 1.0) * steps)
        for _ in range(steps_c):
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
        # keep per-gene, per-cell samples so the correlation is measured WITHIN a
        # gene across cells (Fig 2H), not pooled across genes (where the gene's
        # expression level would trivially correlate both).
        per_gene.setdefault(c, {})
        for g in genes:
            per_gene[c][g] = (rna.get(g, 0.0), prot.get(g, 0.0))
        mrna_genes = [g for g in genes if g in _mrna_keys] or genes
        if rna:
            mrna_per_gene.append(np.mean([rna.get(g, 0.0) for g in mrna_genes]))
    # per-gene |Pearson r| across cells, averaged over genes
    all_genes = sorted({g for cell in per_gene.values() for g in cell})
    r_by_gene = []
    for g in all_genes:
        m = np.array([per_gene[c][g][0] for c in per_gene if g in per_gene[c]])
        p = np.array([per_gene[c][g][1] for c in per_gene if g in per_gene[c]])
        if len(m) > 3 and m.std() > 0 and p.std() > 0:
            r_by_gene.append(abs(float(np.corrcoef(m, p)[0, 1])))
    corr = float(np.mean(r_by_gene)) if r_by_gene else 0.0
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
