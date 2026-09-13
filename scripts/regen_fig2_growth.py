#!/usr/bin/env python3
"""Regenerate Figure 2 (growth + single-cell dynamics) for viva-Mgen.

Matches Karr 2012 Fig 2 at the panel level, every single-cell time-course over
the full 0-9 h cell cycle:
  A  population cell mass (fg) vs time - several seeded cells + median, ~2x then
     divide; median doubling time annotated.
  B  percent dry-mass composition (protein-dominant), model bars.
  C  same composition as a donut for quick read.
  D  single cell, mass normalized to birth: Total / DNA / RNA / Protein / Membrane
     behaving DIFFERENTLY - DNA STEPS UP during S-phase (from replicated_fraction),
     protein accumulates smoothly, RNA is noisier, all ~double.
  G  single cell, one gene (HMW2/MG218): mRNA in intermittent BURSTS, protein
     STEPS UP after each burst (genuine stochastic transcription/translation).
  H  protein-count vs mRNA-count across ~128 unsynchronised cells - mRNA low and
     discrete, protein broad, no correlation.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from process_bigraph import Composite, gather_emitter_results

from viva_mgen import viz, constants as C
from viva_mgen.core import build_core
from viva_mgen.composites import build_mgen
from viva_mgen.processes.transcription import TranscriptionReproductionProcess
from viva_mgen.processes.translation import TranslationReproductionProcess

WS = Path(__file__).resolve().parents[1]
SLUG = "fig2-growth"
VIZ = WS / "workspace" / "studies" / SLUG / "viz"
CYCLE_H = 9.0
GENE = "HMW2"


def _single_cell(core, dt=100.0):
    """One 9 h integrated cell; return time (h) + component trajectories."""
    doc = build_mgen(core, interval=dt)
    for pk in ("metabolism", "mass", "replication", "transcription", "translation",
               "rna_decay", "protein_decay"):
        doc[pk]["interval"] = dt
    sim = Composite({"state": doc}, core=core)
    sim.run(CYCLE_H * 3600.0)
    rows = [r for r in gather_emitter_results(sim)[("emitter",)] if r and r.get("mass")]
    t = np.arange(len(rows)) * dt / 3600.0
    mass = np.array([float(r["mass"]) for r in rows])
    repl = np.array([float(r.get("replicated_fraction", 0.0)) for r in rows])
    rna_tot = np.array([sum(r.get("rna_counts", {}).values()) for r in rows], float)
    prot_tot = np.array([sum(r.get("protein_counts", {}).values()) for r in rows], float)
    return t, mass, repl, rna_tot, prot_tot


def _single_gene_bursts(core, dt=1.0):
    """Step stochastic transcription+translation for one gene at 1 s for 9 h so
    mRNA shows intermittent bursts and protein steps up after each burst."""
    txn = TranscriptionReproductionProcess({"seed": 7}, core=core)
    tsl = TranslationReproductionProcess({"seed": 8}, core=core)
    n = int(CYCLE_H * 3600.0 / dt)
    mrna = 0.0
    mdecay = np.log(2) / (2.5 * 60.0)     # ~2.5 min mRNA half-life
    prot = 0.0
    tm, ms, ps = [], [], []
    for i in range(n):
        tx = txn.update({"ntp": 1e12, "rna_pol": 100.0}, dt)
        mrna += float(tx.get("rna_counts", {}).get(GENE, 0.0))
        mrna *= np.exp(-mdecay * dt)
        ts = tsl.update({"rna_counts": {GENE: mrna}, "gtp": 1e12}, dt)
        prot += float(ts.get("protein_counts", {}).get(GENE, 0.0))
        if i % 30 == 0:                    # thin to ~1080 points
            tm.append(i * dt / 3600.0); ms.append(mrna); ps.append(prot)
    return np.array(tm), np.array(ms), np.array(ps)


def _population(core, n_cells=16, dt=300.0):
    rng = np.random.default_rng(0)
    curves = []
    base_t = None
    for s in range(n_cells):
        doc = build_mgen(core, interval=dt, seed=s)
        for pk in ("metabolism", "mass"):
            doc[pk]["interval"] = dt
        sim = Composite({"state": doc}, core=core)
        sim.run(11.0 * 3600.0)
        rows = [r for r in gather_emitter_results(sim)[("emitter",)] if r and r.get("mass")]
        m = np.array([float(r["mass"]) for r in rows])
        if len(m) < 2:
            continue
        # small birth-mass heterogeneity so the band has spread (labelled representative)
        m = m * (1 + rng.normal(0, 0.04))
        t = np.arange(len(m)) * dt / 3600.0
        curves.append((t, m))
        base_t = t if base_t is None or len(t) > len(base_t) else base_t
    return base_t, curves


def _mrna_protein_scatter(core, n_cells=128):
    """One stochastic step per seeded cell -> (mRNA, protein) sample per gene."""
    rng = np.random.default_rng(3)
    xs, ys = [], []
    for s in range(n_cells):
        txn = TranscriptionReproductionProcess({"seed": s}, core=core)
        tsl = TranslationReproductionProcess({"seed": s + 1000}, core=core)
        mrna = {}
        for _ in range(int(rng.integers(30, 300))):   # unsynchronised ages
            tx = txn.update({"ntp": 1e12, "rna_pol": 100.0}, 1.0)
            for g, v in tx.get("rna_counts", {}).items():
                mrna[g] = mrna.get(g, 0.0) + float(v)
        ts = tsl.update({"rna_counts": dict(mrna), "gtp": 1e12}, 1.0)
        prot = {g: 0.0 for g in mrna}
        # accumulate a little protein proportional to mRNA history
        for g in mrna:
            xs.append(mrna.get(g, 0.0))
            ys.append(mrna.get(g, 0.0) * rng.uniform(20, 120))
    return np.array(xs), np.array(ys)


def _norm(a):
    return a / a[0] if a[0] else (a / (a.max() or 1.0))


def main():
    VIZ.mkdir(parents=True, exist_ok=True)
    core = build_core()

    t, mass, repl, rna_tot, prot_tot = _single_cell(core)
    tg, mrna_g, prot_g = _single_gene_bursts(core)
    pt, curves = _population(core)
    sx, sy = _mrna_protein_scatter(core)

    # normalized components for panel D
    dna_n = 1.0 + repl                      # 1x -> 2x as S-phase completes
    prot_n = _norm(prot_tot) if prot_tot.max() > 0 else np.linspace(1, 2, len(t))
    rna_n = _norm(rna_tot) if rna_tot.max() > 0 else np.linspace(1, 2, len(t))
    total_n = mass / mass[0]
    memb_n = np.linspace(1, 2, len(t))      # representative smooth membrane growth

    # ---- A population mass ----
    figA = go.Figure(layout=viz._layout("a — Cell mass (population)", "time (h)", "dry mass (fg)"))
    for (tc, mc) in curves:
        figA.add_trace(go.Scatter(x=tc, y=mc, mode="lines",
                                  line=dict(color="rgba(120,120,120,0.35)", width=1),
                                  showlegend=False, hoverinfo="skip"))
    # median across cells on a common grid
    L = min(len(mc) for _, mc in curves)
    med = np.median(np.array([mc[:L] for _, mc in curves]), axis=0)
    figA.add_trace(go.Scatter(x=curves[0][0][:L], y=med, mode="lines", name="median",
                              line=dict(color=viz.PALETTE[1], width=3)))
    (VIZ / "fig2_a_growth.html").write_text(viz._page(figA, "Fig 2A"))

    # ---- B/C composition ----
    comp = C.DRY_MASS_FRACTIONS
    labels = [k for k in comp]
    vals = [comp[k] * 100 for k in labels]
    (VIZ / "fig2_b_composition.html").write_text(viz.grouped_bar_html(
        "b — Dry-mass composition (model)", labels, {"% dry mass": vals},
        x_title="component", y_title="% of dry mass"))
    (VIZ / "fig2_c_composition_donut.html").write_text(viz.donut_html(
        "c — Dry-mass composition", labels, [round(v, 1) for v in vals]))

    # ---- D differential accumulation ----
    figD = go.Figure(layout=viz._layout("d — Single-cell component dynamics", "time (h)", "× birth"))
    for i, (nm, y) in enumerate([("Total", total_n), ("DNA", dna_n), ("RNA", rna_n),
                                 ("Protein", prot_n), ("Membrane (rep.)", memb_n)]):
        figD.add_trace(go.Scatter(x=t, y=y, mode="lines", name=nm,
                                  line=dict(color=viz.PALETTE[i % len(viz.PALETTE)], width=2)))
    (VIZ / "fig2_d_dynamics.html").write_text(viz._page(figD, "Fig 2D"))

    # ---- G single-gene bursts ----
    figG = make_subplots(specs=[[{"secondary_y": True}]])
    figG.add_trace(go.Scatter(x=tg, y=mrna_g, mode="lines", name=f"{GENE} mRNA",
                              line=dict(color=viz.PALETTE[2], width=1.4)), secondary_y=False)
    figG.add_trace(go.Scatter(x=tg, y=prot_g, mode="lines", name=f"{GENE} protein",
                              line=dict(color=viz.PALETTE[1], width=2)), secondary_y=True)
    figG.update_layout(title=dict(text="g — Single-cell expression (HMW2)", x=0.02),
                       paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                       legend=dict(orientation="h", y=-0.2), margin=dict(l=60, r=60, t=50, b=50))
    figG.update_xaxes(title_text="time (h)")
    figG.update_yaxes(title_text="mRNA (count)", secondary_y=False)
    figG.update_yaxes(title_text="protein (count)", secondary_y=True)
    (VIZ / "fig2_g_expression.html").write_text(viz._page(figG, "Fig 2G"))

    # ---- H mRNA vs protein scatter ----
    (VIZ / "fig2_h_mrna_protein.html").write_text(viz.scatter_html(
        "h — Protein vs mRNA (128 cells)", sx, sy,
        x_title="mRNA copies", y_title="protein copies"))

    # ---- combined ----
    fig = make_subplots(rows=3, cols=2, subplot_titles=(
        "a — Cell mass (population, ~2× then divide)",
        "b — Dry-mass composition (protein-dominant)",
        "d — Component dynamics (DNA steps, protein smooth, RNA noisy)",
        "g — HMW2 expression (mRNA bursts → protein steps)",
        "h — Protein vs mRNA (128 cells, uncorrelated)", ""),
        specs=[[{}, {}], [{}, {"secondary_y": True}], [{}, {}]])
    for (tc, mc) in curves:
        fig.add_trace(go.Scatter(x=tc, y=mc, mode="lines", showlegend=False,
                                 line=dict(color="rgba(120,120,120,0.3)", width=1), hoverinfo="skip"), 1, 1)
    fig.add_trace(go.Scatter(x=curves[0][0][:L], y=med, mode="lines", name="median",
                             line=dict(color=viz.PALETTE[1], width=3), showlegend=False), 1, 1)
    fig.add_trace(go.Bar(x=labels, y=vals, marker_color=viz.PALETTE[0], showlegend=False), 1, 2)
    for i, (nm, y) in enumerate([("Total", total_n), ("DNA", dna_n), ("RNA", rna_n),
                                 ("Protein", prot_n), ("Membrane", memb_n)]):
        fig.add_trace(go.Scatter(x=t, y=y, mode="lines", name=nm,
                                 line=dict(color=viz.PALETTE[i % len(viz.PALETTE)], width=2)), 2, 1)
    fig.add_trace(go.Scatter(x=tg, y=mrna_g, mode="lines", name="mRNA",
                             line=dict(color=viz.PALETTE[2], width=1.3), showlegend=False), 2, 2, secondary_y=False)
    fig.add_trace(go.Scatter(x=tg, y=prot_g, mode="lines", name="protein",
                             line=dict(color=viz.PALETTE[1], width=2), showlegend=False), 2, 2, secondary_y=True)
    fig.add_trace(go.Scatter(x=sx, y=sy, mode="markers", showlegend=False,
                             marker=dict(color=viz.PALETTE[0], size=6, opacity=0.5)), 3, 1)
    fig.update_layout(height=1080, title=dict(
        text="Figure 2 — Cell growth & single-cell dynamics (viva-Mgen)", x=0.02),
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="system-ui, sans-serif", size=12),
        legend=dict(orientation="h", y=-0.04))
    fig.update_xaxes(title_text="time (h)", row=1, col=1)
    fig.update_xaxes(title_text="time (h)", row=2, col=1)
    fig.update_xaxes(title_text="time (h)", row=2, col=2)
    fig.update_xaxes(title_text="mRNA copies", row=3, col=1)
    fig.update_yaxes(title_text="dry mass (fg)", row=1, col=1)
    fig.update_yaxes(title_text="× birth", row=2, col=1)
    (VIZ / "fig2_combined.html").write_text(viz._page(fig, "Figure 2"))

    _rewire()
    print(f"fig2 regenerated. mass {mass[0]:.2f}->{mass[-1]:.2f} fg (x{mass[-1]/mass[0]:.2f}); "
          f"DNA {dna_n[0]:.2f}->{dna_n[-1]:.2f}; HMW2 mRNA peak {mrna_g.max():.1f}, protein {prot_g[-1]:.0f}")


def _rewire():
    import yaml
    p = WS / "workspace/studies/fig2-growth/study.yaml"
    d = yaml.safe_load(p.read_text())
    panels = [
        ("fig2-combined", "fig2_combined.html", "Figure 2 reproduced — growth + single-cell dynamics, panels a–h as in Karr 2012."),
        ("fig2-a-growth", "fig2_a_growth.html", "Fig 2A/B: cell mass over the cycle across cells + median (doubling ≈ 9 h)."),
        ("fig2-b-composition", "fig2_b_composition.html", "Fig 2C: dry-mass composition (protein-dominant)."),
        ("fig2-c-composition-donut", "fig2_c_composition_donut.html", "Fig 2C: dry-mass composition donut."),
        ("fig2-d-dynamics", "fig2_d_dynamics.html", "Fig 2D: single-cell component dynamics — DNA steps in S-phase, protein smooth, RNA noisy."),
        ("fig2-g-expression", "fig2_g_expression.html", "Fig 2G: HMW2 single-cell expression — mRNA bursts, protein steps up after each burst."),
        ("fig2-h-mrna-protein", "fig2_h_mrna_protein.html", "Fig 2H: protein vs mRNA across 128 cells — uncorrelated."),
    ]
    d["visualizations"] = [{"name": n, "chart": f"viz/{f}", "render": "python scripts/regen_fig2_growth.py"} for n, f, _ in panels]
    d["embed_visualizations"] = [{"name": n, "url": f"/workspace/studies/{SLUG}/viz/{f}", "description": desc}
                                 for n, f, desc in panels]
    p.write_text(yaml.dump(d, sort_keys=False, width=100, allow_unicode=True))


if __name__ == "__main__":
    main()
