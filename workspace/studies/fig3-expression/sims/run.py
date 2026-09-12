#!/usr/bin/env python3
"""Canonical run for the ``fig3-expression`` study (Fig 3 / 2G-2H of Karr et al. 2012).

Runs the stochastic transcription + translation + RNA/protein decay processes on
the representative gene panel through shared bigraph stores for one hour, and
checks that the reproduction recapitulates the *qualitative* single-cell
expression phenomena of Fig 2G/2H: bursty, low-copy mRNA and protein that
accumulates to much higher, more stable copy numbers (the mRNA↔protein
decoupling). Renders an interactive expression time-series and a per-gene
mRNA-vs-protein comparison, and records the run in ``.pbg/runs.jsonl``.

REDUCED FIDELITY: a 6-gene representative panel, not the full ~480-gene set. It
reproduces the qualitative decoupling, NOT the paper's quantitative t50 = 18 min
mRNA half-life-of-expression value (that requires the full gene set) — see the
study conclusion for this divergence.
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
from viva_mgen.composites import build_mgen
from viva_mgen import viz
from viva_mgen.expression_defaults import DEFAULT_GENES
from vivarium_workbench.lib.run_log import append_run_event

SPEC_ID = "viva_mgen.composites.mgen.mycoplasma_genitalium"
STUDY_SLUG = "fig3-expression"
INVESTIGATION_SLUG = "mgen"


def _run(n_seconds=3600.0, dt=1.0, seed=0):
    core = build_core()
    doc = build_mgen(core, seed=seed)
    for _pk in ("metabolism","mass","replication"):
        doc[_pk]["interval"] = n_seconds  # run once; fig3 measures expression only
    # The map[float] store apply accumulates deltas only into keys that already
    # exist; a bare {} drops every synthesized species. Seed the panel gene keys
    # at 0.0 so the stochastic synthesis/decay deltas land and accumulate.
    doc["stores"]["rna_counts"] = {g: 0.0 for g in DEFAULT_GENES}
    doc["stores"]["protein_counts"] = {g: 0.0 for g in DEFAULT_GENES}
    sim = Composite({"state": doc}, core=core)
    sim.run(n_seconds)
    rows = gather_emitter_results(sim)[("emitter",)]
    rows = [r for r in rows if r is not None]
    return rows, dt


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
        rows, dt = _run()

        n_genes = len(DEFAULT_GENES)
        # per-row totals across the gene panel
        total_mrna = np.array([sum((r.get("rna_counts") or {}).values()) for r in rows], dtype=float)
        total_protein = np.array([sum((r.get("protein_counts") or {}).values()) for r in rows], dtype=float)
        t_min = np.arange(len(rows)) * dt / 60.0  # minutes

        final_rna = rows[-1].get("rna_counts") or {}
        final_protein = rows[-1].get("protein_counts") or {}

        # observables
        mean_mrna_per_gene = float(total_mrna.mean() / n_genes)
        final_total_mrna = float(total_mrna[-1])
        final_total_protein = float(total_protein[-1])
        protein_to_mrna_ratio = float(final_total_protein / max(1.0, final_total_mrna))
        fraction_genes_expressed = float(
            sum(1 for g in DEFAULT_GENES if float(final_protein.get(g, 0.0)) > 0.0) / n_genes)

        print(f"mean_mrna_per_gene       = {mean_mrna_per_gene:.3f}  (low/bursty, expect ~0-2)")
        print(f"final_total_mrna         = {final_total_mrna:.1f}")
        print(f"final_total_protein      = {final_total_protein:.1f}  (>> mRNA)")
        print(f"protein_to_mrna_ratio    = {protein_to_mrna_ratio:.1f}  (>> 1, Fig 2H decoupling)")
        print(f"fraction_genes_expressed = {fraction_genes_expressed:.3f}  (expect ~1.0)")

        viz_dir = STUDY_DIR / "viz"
        viz_dir.mkdir(parents=True, exist_ok=True)

        # 1) expression time series — total mRNA vs total protein over the hour
        (viz_dir / "expression_timeseries.html").write_text(viz.line_series_html(
            "Gene expression dynamics: bursty mRNA vs accumulating protein (Fig 2G)",
            t_min, {"total mRNA": total_mrna, "total protein": total_protein},
            x_title="time (min)", y_title="copy number"))

        # 2) per-gene final protein vs mRNA — grouped bars showing protein >> mRNA
        genes = list(DEFAULT_GENES)
        (viz_dir / "protein_vs_mrna.html").write_text(viz.grouped_bar_html(
            "mRNA vs protein copy number by gene (Fig 2H)",
            genes,
            {"mRNA": [float(final_rna.get(g, 0.0)) for g in genes],
             "protein": [float(final_protein.get(g, 0.0)) for g in genes]},
            x_title="gene", y_title="final copy number"))

        observables = {
            "mean_mrna_per_gene": mean_mrna_per_gene,
            "final_total_mrna": final_total_mrna,
            "final_total_protein": final_total_protein,
            "protein_to_mrna_ratio": protein_to_mrna_ratio,
            "fraction_genes_expressed": fraction_genes_expressed,
        }
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
