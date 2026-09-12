#!/usr/bin/env python3
"""Canonical run for the ``fig7-kinetic-parameters`` study (Fig 7E-G of Karr et al. 2012).

Reproduces the quantitative gene-disruption characterization of Figure 7E-G:
predicted growth rate depends on an enzyme's kinetic parameter (kcat/Vmax),
producing a SIGMOIDAL growth-vs-kcat curve that rises from near-zero at low
enzyme activity and saturates at the wild-type rate once the enzyme is no longer
rate-limiting. In the paper this dependence is what lets tuning kcat reconcile
model and experiment for lpdA, deoD and thyA.

FIDELITY: REDUCED / representative. kcat is modelled as a multiplicative scale
on a target reaction's FBA flux bound (a Vmax proxy) via
``MetabolismFbaReproductionProcess(config={"reaction_bound_scale": {rid: s}})``,
NOT explicit enzyme kinetics. The specific lpdA/deoD/thyA genetics are not
individually modelled; instead we identify, programmatically, a genuinely
growth-limiting reaction in the iPS189 network and sweep its flux bound. The
selected target is the leucine-uptake reaction ``EX_leu_DASH_L_e`` (an
amino-acid transporter Vmax proxy): M. genitalium is auxotrophic for most amino
acids, so amino-acid provisioning is a real growth-limiting kinetic step — the
same qualitative role the paper's biosynthetic enzymes play.
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
from viva_mgen import kb, viz
from viva_mgen.composites.genetics import fig7_kinetic_parameters
from viva_mgen.processes.metabolism import MetabolismFbaReproductionProcess
from vivarium_workbench.lib.run_log import append_run_event

SPEC_ID = "viva_mgen.composites.genetics.fig7_kinetic_parameters"
STUDY_SLUG = "fig7-kinetic-parameters"
INVESTIGATION_SLUG = "mgen"

# kcat proxy sweep: 0.01 .. 10 x wild-type flux bound (16 log-spaced points)
SWEEP_SCALES = np.logspace(-2, 1, 16)

# Preferred primary target (an amino-acid transporter Vmax proxy) — confirmed
# growth-limiting-and-saturating by the programmatic scan below; falls back to
# the widest-range limiting reaction if it is ever not found.
PREFERRED_TARGET = "EX_leu_DASH_L_e"

# central-carbon / ATP / nucleotide shortlist the scan seeds itself with, on top
# of every reaction that carries flux at the wild-type optimum
CANDIDATE_SHORTLIST = ["ATPS4r", "PGK", "PYK", "ENO", "PGI", "PFK", "FBA", "GAPD",
                       "TPI", "PGM", "HEX1", "NDPK1", "NDPK2", "GK1", "PRPPS"]


def _growth_fraction(core, rid, scale):
    """Growth fraction with the target reaction's flux bound scaled by ``scale``."""
    proc = MetabolismFbaReproductionProcess(
        config={"reaction_bound_scale": {rid: float(scale)}}, core=core)
    return float(proc.update({"nutrient_scale": 1.0}, 1.0)["growth_fraction"])


def _find_targets(core):
    """Programmatically find growth-LIMITING reactions.

    A reaction is limiting-and-saturating if throttling its bound to 0.05x drops
    growth well below the wild-type (bound is binding) AND growth has already
    recovered to the wild-type plateau by 0.2x (so the sweep shows a real
    rise-then-plateau sigmoid, not an unbounded ramp). Returns
    ``(primary, secondaries, ranked)`` where ranked is [(rid, dynamic_range), ...].
    """
    model = kb.load_metabolic_model()
    sol = model.optimize()
    active = [r.id for r in model.reactions
              if abs(sol.fluxes.get(r.id, 0.0)) > 1e-6]
    candidates, seen = [], set()
    for rid in CANDIDATE_SHORTLIST + active:
        if rid in seen or rid not in model.reactions:
            continue
        seen.add(rid)
        candidates.append(rid)

    ranked = []
    for rid in candidates:
        g05 = _growth_fraction(core, rid, 0.05)
        g20 = _growth_fraction(core, rid, 0.2)
        g10 = _growth_fraction(core, rid, 1.0)
        is_limiting = g05 < g10 - 0.02              # bound binds when throttled
        saturates = g10 > 0 and g20 >= 0.97 * g10   # plateau reached before top
        if is_limiting and saturates:
            ranked.append((rid, g10 - g05))
    ranked.sort(key=lambda x: x[1], reverse=True)

    ids = [rid for rid, _ in ranked]
    primary = PREFERRED_TARGET if PREFERRED_TARGET in ids else (ids[0] if ids else PREFERRED_TARGET)
    secondaries = [rid for rid in ids if rid != primary][:2]
    return primary, secondaries, ranked


def _sweep(core, rid):
    return np.array([_growth_fraction(core, rid, s) for s in SWEEP_SCALES])


def main() -> int:
    run_id = uuid.uuid4().hex
    core = build_core()

    # discover the target(s) first so the started-event params are honest
    primary, secondaries, ranked = _find_targets(core)

    append_run_event(WORKSPACE_ROOT, {
        "run_id": run_id, "event": "started", "spec_id": SPEC_ID,
        "label": STUDY_SLUG, "started_at": time.time(), "status": "running",
        "n_steps": len(SWEEP_SCALES), "emitter": "ram", "origin": "canonical_run",
        "study_slug": STUDY_SLUG, "investigation_slug": INVESTIGATION_SLUG,
        "params": {"reaction_id": primary, "bound_scale": 1.0},
    })
    try:
        targets = [primary] + secondaries
        sweeps = {rid: _sweep(core, rid) for rid in targets}
        g = sweeps[primary]

        # sanity: run the declared composite once at the wild-type bound
        doc = fig7_kinetic_parameters(core, reaction_id=primary, bound_scale=1.0)
        sim = Composite({"state": doc}, core=core)
        sim.run(1.0)
        crows = gather_emitter_results(sim)[("emitter",)]
        cgf = [r["growth_fraction"] for r in crows if "growth_fraction" in r]
        baseline_gf = float(cgf[-1]) if cgf else 0.0

        # observables (computed on the primary target's sweep)
        eps = 1e-3
        monotonic = all(g[i + 1] >= g[i] - eps for i in range(len(g) - 1))
        saturates = g[-1] >= 0.97 * float(g.max())
        dynamic_range = float(g.max() - g.min())

        observables = {
            "target_reaction": primary,
            "growth_is_monotonic_in_kcat": 1.0 if monotonic else 0.0,
            "growth_saturates": 1.0 if saturates else 0.0,
            "dynamic_range": dynamic_range,
            "growth_at_min_kcat": float(g.min()),
            "growth_at_max_kcat": float(g[-1]),
            "baseline_growth_fraction": baseline_gf,
        }

        # ---- print sweep table + observables --------------------------------
        print(f"target_reaction (primary)     = {primary}")
        print(f"secondary series              = {', '.join(secondaries) or '(none)'}")
        print(f"n limiting+saturating found   = {len(ranked)}")
        print("\nkcat-proxy sweep (flux-bound scale -> growth_fraction):")
        header = "  scale   " + "".join(f"{rid[:12]:>13s}" for rid in targets)
        print(header)
        for i, s in enumerate(SWEEP_SCALES):
            row = f"  {s:7.3f}  " + "".join(f"{sweeps[rid][i]:13.3f}" for rid in targets)
            print(row)
        print("\nobservables:")
        for k, v in observables.items():
            print(f"  {k:30s} = {v}")

        # ---- viz: kcat -> growth sigmoid (Fig 7E) ---------------------------
        viz_dir = STUDY_DIR / "viz"
        viz_dir.mkdir(parents=True, exist_ok=True)
        series = {rid: sweeps[rid].tolist() for rid in targets}
        html = viz.line_series_html(
            "Growth rate vs kinetic parameter (kcat proxy) (Fig 7E)",
            SWEEP_SCALES.tolist(), series,
            x_title="flux-bound scale (kcat / Vmax proxy, x wild-type)",
            y_title="growth fraction", modes="lines+markers")
        # render the kcat axis logarithmically (a kcat sweep spans decades) without
        # touching the shared viz helper — inject a single layout tweak pre-plot
        html = html.replace(
            "Plotly.newPlot(",
            "spec.layout.xaxis=Object.assign(spec.layout.xaxis||{},{type:'log'});\n"
            "Plotly.newPlot(", 1)
        (viz_dir / "kcat_growth_curve.html").write_text(html)
    except Exception:
        append_run_event(WORKSPACE_ROOT, {
            "run_id": run_id, "event": "completed", "completed_at": time.time(),
            "n_steps": 0, "status": "failed"})
        raise

    append_run_event(WORKSPACE_ROOT, {
        "run_id": run_id, "event": "completed", "completed_at": time.time(),
        "n_steps": len(SWEEP_SCALES), "status": "completed",
        "observables": observables})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
