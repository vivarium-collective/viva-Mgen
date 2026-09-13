#!/usr/bin/env python3
"""Canonical run for the ``fig5-energy`` study (Fig 5 of Karr et al. 2012).

Reproduces Figure 5 — the global distribution of cellular energy — at REDUCED
fidelity. Two readouts:

* Fig 5A (synthesis rates): FBA metabolism produces ATP and GTP; the LP flux for
  ATP synthase (ATPS4r) far exceeds the GTP-producing (NDPK/GK) flux, so
  ATP > GTP > 0, echoing Fig 5A's ordering (ATP/GTP dominate the carriers).

* Fig 5D (energy budget by process): the paper's budget is dominated by
  translation (~29%), then tRNA-aminoacylation (~15%), then transcription
  (~7%), with a large unaccounted remainder (~44%). This reduced model tracks
  only the two expression processes actually wired in the composite —
  transcription (NTP) and translation (GTP) — and computes their shares of the
  MODELED expression energy by directly stepping the two processes for one hour
  and summing the magnitude of their returned NTP/GTP consumption. We reproduce
  the ORDERING (translation >> transcription) and the rough structure, not the
  exact 29/15/7/44 split: NAD/FAD/NADP carriers and tRNA-aminoacylation are not
  modeled here.
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
from viva_mgen.processes.transcription import TranscriptionReproductionProcess
from viva_mgen.processes.translation import TranslationReproductionProcess
from viva_mgen import viz
from vivarium_workbench.lib.run_log import append_run_event

SPEC_ID = "viva_mgen.composites.mgen.mycoplasma_genitalium"
STUDY_SLUG = "fig5-energy"
INVESTIGATION_SLUG = "mgen"

# Fig 5D representative "unaccounted" share (NAD/FAD/NADP + tRNA-aminoacylation +
# everything else) — used ONLY to echo the paper's donut structure, clearly
# labelled representative; it is not a measured quantity of this reduced model.
FIG5D_UNACCOUNTED = 0.44


def _run_metabolism(core, n_seconds=3600.0, dt=60.0):
    """Integrated cell for one hour; time-averaged ATP/GTP production flux."""
    doc = build_mgen(core, interval=dt)
    doc["metabolism"]["interval"] = dt
    doc["mass"]["interval"] = dt
    for _pk in ("replication","transcription","translation","rna_decay","protein_decay"):
        doc[_pk]["interval"] = dt
    # expression stays at 1.0 s (its native stochastic step)
    sim = Composite({"state": doc}, core=core)
    sim.run(n_seconds)
    rows = gather_emitter_results(sim)[("emitter",)]
    atp = np.array([float(r.get("atp_production", 0.0)) for r in rows])
    gtp = np.array([float(r.get("gtp_production", 0.0)) for r in rows])
    # drop the pre-solve initial emit (production still 0.0)
    keep = atp > 0.0
    atp, gtp = atp[keep], gtp[keep]
    return float(atp.mean()), float(gtp.mean()), rows


def _run_expression_accounting(core, n_seconds=3600, dt=1.0):
    """Step transcription + translation directly for one hour and sum the
    magnitude of their NTP/GTP consumption. Transcription feeds its synthesized
    mRNA forward into an accumulating pool that translation reads each step
    (decay omitted for the energy tally, per the reduced accounting). Pools are
    kept large each step so consumption is never supply-capped."""
    txn = TranscriptionReproductionProcess({"seed": 0}, core=core)
    tsl = TranslationReproductionProcess({"seed": 1}, core=core)
    rna_counts: dict[str, float] = {}
    ntp_used_total = 0.0
    gtp_used_total = 0.0
    for _ in range(int(n_seconds)):
        tx = txn.update({"ntp": 1e12, "rna_pol": 100.0}, dt)
        for gene, n in tx.get("rna_counts", {}).items():
            rna_counts[gene] = rna_counts.get(gene, 0.0) + float(n)
        ntp_used_total += abs(float(tx.get("ntp", 0.0)))

        ts = tsl.update({"rna_counts": dict(rna_counts), "gtp": 1e12}, dt)
        gtp_used_total += abs(float(ts.get("gtp", 0.0)))
    return ntp_used_total, gtp_used_total


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
        atp_production_rate, gtp_production_rate, rows = _run_metabolism(core)
        transcription_ntp_used, translation_gtp_used = _run_expression_accounting(core)

        assert atp_production_rate > gtp_production_rate > 0.0, (
            f"expected atp>gtp>0, got atp={atp_production_rate}, gtp={gtp_production_rate}")

        transcription_energy = transcription_ntp_used
        translation_energy = translation_gtp_used
        modeled_total = translation_energy + transcription_energy or 1.0
        translation_share = translation_energy / modeled_total
        transcription_share = 1.0 - translation_share
        atp_to_gtp_ratio = atp_production_rate / max(1e-9, gtp_production_rate)

        observables = {
            "atp_production_rate": atp_production_rate,
            "gtp_production_rate": gtp_production_rate,
            "atp_to_gtp_ratio": atp_to_gtp_ratio,
            "transcription_energy": transcription_energy,
            "translation_energy": translation_energy,
            "translation_share": translation_share,
            "transcription_share": transcription_share,
        }

        print(f"atp_production_rate  = {atp_production_rate:.4f}")
        print(f"gtp_production_rate  = {gtp_production_rate:.4f}")
        print(f"atp_to_gtp_ratio     = {atp_to_gtp_ratio:.4f}  (target >1)")
        print(f"transcription_energy = {transcription_energy:.0f} NTP")
        print(f"translation_energy   = {translation_energy:.0f} GTP")
        print(f"translation_share    = {translation_share:.4f}  (of modeled expression energy)")
        print(f"transcription_share  = {transcription_share:.4f}")

        viz_dir = STUDY_DIR / "viz"
        viz_dir.mkdir(parents=True, exist_ok=True)

        # Fig 5D — energy allocation. Modeled shares (translation, transcription)
        # scaled to fill the accounted 56%, plus a representative 44% unaccounted
        # slice to echo the paper's donut structure.
        accounted = 1.0 - FIG5D_UNACCOUNTED
        (viz_dir / "energy_allocation.html").write_text(viz.donut_html(
            "Energy allocation across processes (Fig 5D)",
            ["translation", "transcription",
             "other/unaccounted (representative, ~44% per Fig 5D)"],
            [translation_share * accounted, transcription_share * accounted,
             FIG5D_UNACCOUNTED]))

        # Fig 5A — ATP vs GTP synthesis rate.
        (viz_dir / "synthesis_rates.html").write_text(viz.grouped_bar_html(
            "ATP and GTP synthesis (Fig 5A)",
            ["ATP", "GTP"],
            {"synthesis rate (FBA flux)": [atp_production_rate, gtp_production_rate]},
            x_title="energy carrier", y_title="production flux (iPS189 units)"))
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
