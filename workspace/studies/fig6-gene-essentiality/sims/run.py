#!/usr/bin/env python3
"""Canonical run for the ``fig6-gene-essentiality`` study (Fig 6A of Karr et al. 2012).

Reproduces Figure 6A: the single-gene-disruption essentiality phenotype, model
vs experiment. For every gene in the iPS189 metabolic reconstruction that also
carries a reference essentiality call, an in-silico single-gene knockout is run
via FBA and the gene is classified essential (growth collapses) or non-essential.
Model predictions are cross-tabulated against the reference calls into a 2x2
confusion matrix and scored (accuracy, sensitivity, specificity).

FIDELITY: HIGH over the ~125 metabolic genes iPS189 covers — the paper notes
metabolic-gene disruptions are the most debilitating. Essentiality is assessed
over this metabolic subset, not all 525 M. genitalium genes. A gene is called
model-essential when its single knockout drops FBA growth_fraction below 0.05.
Reference calls are the ``essential_ref`` column of ``datasets/genes.csv``,
extracted from the Karr 2012 repo's allDeletionSimulations.json.
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

from process_bigraph import Composite, gather_emitter_results

from viva_mgen.core import build_core
from viva_mgen import kb, viz
from viva_mgen.composites.genetics import fig6_gene_essentiality
from viva_mgen.processes.metabolism import MetabolismFbaReproductionProcess
from vivarium_workbench.lib.run_log import append_run_event

SPEC_ID = "viva_mgen.composites.genetics.fig6_gene_essentiality"
STUDY_SLUG = "fig6-gene-essentiality"
INVESTIGATION_SLUG = "mgen"

# growth_fraction below this counts as an essential-gene (growth-collapse) call
ESSENTIAL_THRESHOLD = 0.05


def _scan():
    """Single-gene-knockout FBA scan over metabolic genes with a reference call."""
    core = build_core()
    ref = kb.gene_essentiality_reference()  # normalized id -> bool essential
    # (normalized_id, model_gene_id) for metabolic genes that have a reference call
    genes = [(kb.normalize_gene_id(mid), mid) for mid in kb.metabolic_gene_ids()]
    genes = [(nid, mid) for nid, mid in genes if nid in ref]

    tp = fn = fp = tn = 0
    for nid, mid in genes:
        proc = MetabolismFbaReproductionProcess(config={"disrupted_genes": [mid]}, core=core)
        out = proc.update({"nutrient_scale": 1.0}, 1.0)
        model_essential = float(out["growth_fraction"]) < ESSENTIAL_THRESHOLD
        ref_essential = bool(ref[nid])
        if ref_essential and model_essential:
            tp += 1
        elif ref_essential and not model_essential:
            fn += 1
        elif not ref_essential and model_essential:
            fp += 1
        else:
            tn += 1
    return tp, fn, fp, tn


def _baseline_composite_run(gene="MG_023"):
    """Run the fig6_gene_essentiality composite once for a representative essential
    gene to record a canonical baseline and confirm growth collapses."""
    core = build_core()
    doc = fig6_gene_essentiality(core, disrupted_gene=gene)
    sim = Composite({"state": doc}, core=core)
    sim.run(1.0)
    rows = gather_emitter_results(sim)[("emitter",)]
    gf = [r["growth_fraction"] for r in rows if "growth_fraction" in r]
    return float(gf[-1]) if gf else 0.0


def main() -> int:
    run_id = uuid.uuid4().hex
    append_run_event(WORKSPACE_ROOT, {
        "run_id": run_id, "event": "started", "spec_id": SPEC_ID,
        "label": STUDY_SLUG, "started_at": time.time(), "status": "running",
        "n_steps": 1, "emitter": "ram", "origin": "canonical_run",
        "study_slug": STUDY_SLUG, "investigation_slug": INVESTIGATION_SLUG,
        "params": {"disrupted_gene": "MG_023"},
    })
    try:
        tp, fn, fp, tn = _scan()
        total = tp + fn + fp + tn
        accuracy = (tp + tn) / total if total else 0.0
        n_true_essential = tp + fn
        n_true_nonessential = fp + tn
        sensitivity = tp / (tp + fn) if (tp + fn) else 0.0
        specificity = tn / (tn + fp) if (tn + fp) else 0.0
        matrix = [[tp, fn], [fp, tn]]

        # canonical composite baseline: a representative essential gene collapses growth
        baseline_gf = _baseline_composite_run("MG_023")

        observables = {
            "essentiality_accuracy": accuracy,
            "n_genes_tested": total,
            "n_true_essential": n_true_essential,
            "n_true_nonessential": n_true_nonessential,
            "sensitivity": sensitivity,
            "specificity": specificity,
            "baseline_growth_fraction": baseline_gf,
        }

        print(f"n_genes_tested        = {total}  (metabolic genes with a reference call)")
        print(f"n_true_essential      = {n_true_essential}")
        print(f"n_true_nonessential   = {n_true_nonessential}")
        print(f"essentiality_accuracy = {accuracy:.3f}  (paper overall ~0.79; iPS189 metabolic ~0.87)")
        print(f"sensitivity (TPR)     = {sensitivity:.3f}")
        print(f"specificity (TNR)     = {specificity:.3f}")
        print(f"baseline MG_023 KO growth_fraction = {baseline_gf:.4f}  (essential; should collapse)")
        print("confusion matrix [[TP, FN], [FP, TN]] "
              "(rows=ref [ess, non-ess], cols=model [ess, non-ess]):")
        print(f"  ref essential     : model-ess {tp:4d}   model-noness {fn:4d}")
        print(f"  ref non-essential : model-ess {fp:4d}   model-noness {tn:4d}")

        viz_dir = STUDY_DIR / "viz"
        viz_dir.mkdir(parents=True, exist_ok=True)
        (viz_dir / "confusion_matrix.html").write_text(viz.confusion_matrix_html(
            title="Gene-essentiality prediction: model vs reference (Fig 6A)",
            matrix=matrix))
        (viz_dir / "accuracy_summary.html").write_text(viz.grouped_bar_html(
            title="Essentiality prediction performance (Fig 6A)",
            categories=["accuracy", "sensitivity", "specificity"],
            groups={"viva-mGen (% )": [accuracy * 100.0,
                                       sensitivity * 100.0,
                                       specificity * 100.0]},
            x_title="metric", y_title="percent"))
    except Exception:
        append_run_event(WORKSPACE_ROOT, {
            "run_id": run_id, "event": "completed", "completed_at": time.time(),
            "n_steps": 0, "status": "failed"})
        raise

    append_run_event(WORKSPACE_ROOT, {
        "run_id": run_id, "event": "completed", "completed_at": time.time(),
        "n_steps": total, "status": "completed", "observables": observables})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
