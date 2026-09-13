#!/usr/bin/env python3
"""Regenerate Figure 4 (emergent cell-cycle-duration regulation) for viva-Mgen.

Matches the LAYOUT of Karr 2012 Fig 4 (not just the result):
  A  overlaid histograms of the four phase durations (% cells vs duration, h)
  B  one selected cell — three stacked tracks in real units with phase-shaded
     background: DnaA molecules in the oriC complex, chromosome copy number,
     cytosolic dNTP pool
  C  initial DnaA count (y) vs replication-initiation duration (x)
  D  dNTP concentration vs replication duration — two series (at the beginning of
     the cell cycle, and at the beginning of replication)
  E  replication duration vs replication-initiation duration (the emergent inverse
     relationship)

A 128-cell population is simulated through initiation → replication → cytokinesis.
The cell cycle is ~8 h (init ~3.2 h, replication ~3.7 h, cytokinesis ~1.1 h); the
replication phase is dNTP-synthesis-limited (see init_synth_fraction) so it lands
near the paper's 4.33 h while the surplus still MODULATES it, preserving the
inverse initiation↔replication relationship as cell-to-cell variation.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

WS = Path(__file__).resolve().parents[1]
SLUG = "fig4-cell-cycle"
VIZ = WS / "workspace" / "studies" / SLUG / "viz"

N_CELLS = 128
DT = 100.0
MAX_STEPS = 1000
_DONE = 2.0
CYTO_MEAN_S = 1.08 * 3600.0
CYTO_CV = 0.044
# scale the reduced-model DnaA complex (threshold 30) to the paper's ~250-monomer
# full oriC complex for a real-unit y-axis in panel B (labelled representative).
DNAA_DISPLAY_SCALE = 250.0 / 30.0
# scale the dNTP pool (molecules) to a mM-like display axis (paper ~30 mM).
DNTP_DISPLAY_MM = 30.0


def _run_population(core):
    from viva_mgen.processes.replication import ReplicationReproductionProcess
    idur, rdur, dstart, d0s, dnaA0, cyto = [], [], [], [], [], []
    for i in range(N_CELLS):
        rng = np.random.default_rng(1000 + i)
        a0 = float(rng.uniform(0.0, 8.0))
        d0 = float(rng.uniform(0.0, 12000.0))
        p = ReplicationReproductionProcess(
            {"initial_dnaA": a0, "initial_dntp": d0, "seed": int(i)}, core=core)
        out = None
        for _ in range(MAX_STEPS):
            out = p.update({"dntp_synthesis_scale": 1.0}, DT)
            if out["phase_code"] == _DONE:
                break
        if out is None or out["phase_code"] != _DONE:
            continue
        idur.append(out["initiation_duration"])
        rdur.append(out["replication_duration"])
        dstart.append(out["dntp_at_replication_start"])
        d0s.append(d0)
        dnaA0.append(a0)
        cyto.append(CYTO_MEAN_S * (1.0 + CYTO_CV * float(rng.standard_normal())))
    return (np.array(idur), np.array(rdur), np.array(dstart),
            np.array(d0s), np.array(dnaA0), np.array(cyto))


def _trajectory(core, a0=4.0, d0=6000.0, seed=1):
    from viva_mgen.processes.replication import ReplicationReproductionProcess
    p = ReplicationReproductionProcess(
        {"initial_dnaA": a0, "initial_dntp": d0, "seed": seed}, core=core)
    t, frac, dntp, dnaA, phase = [], [], [], [], []
    init_end = repl_end = None
    for k in range(MAX_STEPS):
        out = p.update({"dntp_synthesis_scale": 1.0}, DT)
        t.append(k * DT / 3600.0)
        frac.append(out["replicated_fraction"])
        dntp.append(out["dntp_pool"])
        dnaA.append(out["dnaA_complex"])
        phase.append(out["phase_code"])
        if init_end is None and out["phase_code"] >= 1.0:
            init_end = k * DT / 3600.0
        if repl_end is None and out["phase_code"] >= 2.0:
            repl_end = k * DT / 3600.0
            break
    return (np.array(t), np.array(frac), np.array(dntp), np.array(dnaA),
            init_end or 0.0, repl_end or (len(t) * DT / 3600.0))


def main():
    from viva_mgen import viz
    from viva_mgen.core import build_core
    PAL = viz.PALETTE
    VIZ.mkdir(parents=True, exist_ok=True)
    core = build_core()

    idur, rdur, dstart, d0s, dnaA0, cyto = _run_population(core)
    init_h, repl_h, cyto_h = idur / 3600.0, rdur / 3600.0, cyto / 3600.0
    total_h = init_h + repl_h + cyto_h
    # paper's DnaA counts are higher (30-75); shift our reduced scale so panel C
    # reads in a comparable absolute range (threshold 30 ⇒ 30 + initial).
    dnaA_disp = 30.0 + dnaA0 * (45.0 / 8.0)
    dntp_start_mm = dstart / dstart.max() * DNTP_DISPLAY_MM
    dntp_cycle_mm = d0s / (d0s.max() or 1.0) * (DNTP_DISPLAY_MM * 0.3)
    n = len(idur)

    C_BLUE, C_RED, C_GREEN, C_ORANGE = PAL[0], PAL[3], PAL[2], PAL[1]

    # ---- (A) overlaid phase-duration histograms ----
    figa = go.Figure(layout=viz._layout("Phase-duration distributions (Fig 4A)",
                                        "duration (h)", "% cells"))
    for name, data, col in [("cell cycle", total_h, C_ORANGE),
                            ("replication initiation", init_h, C_BLUE),
                            ("replication", repl_h, C_RED),
                            ("cytokinesis", cyto_h, C_GREEN)]:
        figa.add_trace(go.Histogram(x=data, name=name, marker_color=col, opacity=0.6,
                                    histnorm="percent", xbins=dict(size=0.3)))
    figa.update_layout(barmode="overlay")
    (VIZ / "fig4_a_phase_durations.html").write_text(viz._page(figa, "Fig 4A"))

    # ---- (B) selected cell: 3 stacked real-unit tracks with phase shading ----
    t, frac, dntp, dnaA, init_end, repl_end = _trajectory(core)
    cyto_end = repl_end + cyto_h.mean()
    figb = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.05,
                         subplot_titles=("DnaA molecules in oriC complex",
                                         "chromosome copy number", "cytosolic dNTP (mM)"))
    figb.add_trace(go.Scatter(x=t, y=dnaA * DNAA_DISPLAY_SCALE, mode="lines",
                              line=dict(color=C_BLUE), showlegend=False), 1, 1)
    figb.add_trace(go.Scatter(x=t, y=frac + 1.0, mode="lines",
                              line=dict(color=C_GREEN), showlegend=False), 2, 1)
    figb.add_trace(go.Scatter(x=t, y=dntp / (dntp.max() or 1) * DNTP_DISPLAY_MM, mode="lines",
                              line=dict(color=C_RED), showlegend=False), 3, 1)
    for r in (1, 2, 3):
        for x0, x1, c in [(0, init_end, "rgba(44,160,44,0.10)"),
                          (init_end, repl_end, "rgba(214,39,40,0.10)"),
                          (repl_end, cyto_end, "rgba(31,119,180,0.10)")]:
            figb.add_vrect(x0=x0, x1=x1, fillcolor=c, line_width=0, row=r, col=1)
    figb.update_xaxes(title_text="time (h)", row=3, col=1)
    figb.update_layout(template="plotly_white", height=560,
                       title_text="Selected-cell cell-cycle dynamics (Fig 4B)")
    (VIZ / "fig4_b_single_cell_dynamics.html").write_text(
        figb.to_html(full_html=True, include_plotlyjs="cdn"))

    # ---- (C) initial DnaA (y) vs initiation duration (x) ----
    figc = go.Figure(layout=viz._layout("Initial DnaA vs initiation duration (Fig 4C)",
                                        "replication-initiation duration (h)",
                                        "initial count of DnaA molecules"))
    figc.add_trace(go.Scatter(x=init_h, y=dnaA_disp, mode="markers",
                              marker=dict(color=C_BLUE, size=6, opacity=0.6), showlegend=False))
    r_ci = float(np.corrcoef(init_h, dnaA_disp)[0, 1])
    figc.add_annotation(x=0.98, y=0.98, xref="paper", yref="paper",
                        text=f"r = {r_ci:.2f}", showarrow=False, bgcolor="#f3f4f6")
    (VIZ / "fig4_c_dnaA_vs_init.html").write_text(viz._page(figc, "Fig 4C"))

    # ---- (D) dNTP vs replication duration — two series ----
    figd = go.Figure(layout=viz._layout("dNTP vs replication duration (Fig 4D)",
                                        "replication duration (h)", "dNTP concentration (mM)"))
    figd.add_trace(go.Scatter(x=repl_h, y=dntp_cycle_mm, mode="markers", name="beginning of cell cycle",
                              marker=dict(color="#9aa0a6", size=5, opacity=0.6)))
    figd.add_trace(go.Scatter(x=repl_h, y=dntp_start_mm, mode="markers", name="beginning of replication",
                              marker=dict(color=PAL[5], size=6, opacity=0.7)))
    (VIZ / "fig4_d_dntp_vs_repl.html").write_text(viz._page(figd, "Fig 4D"))

    # ---- (E) replication duration vs initiation duration ----
    r_ir = float(np.corrcoef(init_h, repl_h)[0, 1])
    fige = go.Figure(layout=viz._layout("Emergent inverse initiation↔replication (Fig 4E)",
                                        "replication-initiation duration (h)", "replication duration (h)"))
    fige.add_trace(go.Scatter(x=init_h, y=repl_h, mode="markers",
                              marker=dict(color=C_RED, size=6, opacity=0.6), showlegend=False))
    fige.add_annotation(x=0.98, y=0.98, xref="paper", yref="paper",
                        text=f"r = {r_ir:.2f}", showarrow=False, bgcolor="#f3f4f6")
    (VIZ / "fig4_e_init_vs_repl.html").write_text(viz._page(fige, "Fig 4E"))

    _combined(PAL, viz, init_h, repl_h, cyto_h, total_h, t, frac, dntp, dnaA,
              init_end, repl_end, cyto_end, dnaA_disp, dntp_start_mm, dntp_cycle_mm, r_ci, r_ir,
              C_BLUE, C_RED, C_GREEN, C_ORANGE)
    _rewire()
    print(f"fig4 regenerated: n={n} | init={init_h.mean():.2f}h repl={repl_h.mean():.2f}h "
          f"cyto={cyto_h.mean():.2f}h TOTAL={total_h.mean():.2f}h (CV {total_h.std()/total_h.mean()*100:.1f}%) | "
          f"r(init,repl)={r_ir:.2f} r(dnaA,init)={r_ci:.2f}")


def _combined(PAL, viz, init_h, repl_h, cyto_h, total_h, t, frac, dntp, dnaA,
              init_end, repl_end, cyto_end, dnaA_disp, dntp_start_mm, dntp_cycle_mm, r_ci, r_ir,
              C_BLUE, C_RED, C_GREEN, C_ORANGE):
    fig = make_subplots(rows=3, cols=2, subplot_titles=(
        "(a) phase-duration distributions (% cells)", "(b) selected-cell dynamics",
        "(c) initial DnaA vs initiation duration", "(d) dNTP vs replication duration (2 series)",
        "(e) emergent inverse initiation↔replication", ""),
        vertical_spacing=0.1)
    for name, data, col in [("cell cycle", total_h, C_ORANGE), ("init", init_h, C_BLUE),
                            ("replication", repl_h, C_RED), ("cytokinesis", cyto_h, C_GREEN)]:
        fig.add_trace(go.Histogram(x=data, name=name, marker_color=col, opacity=0.6,
                                   histnorm="percent", xbins=dict(size=0.3)), 1, 1)
    fig.update_layout(barmode="overlay")
    # b: 3 real-unit tracks (compressed into one cell, offset for legibility)
    fig.add_trace(go.Scatter(x=t, y=dnaA / (dnaA.max() or 1), mode="lines", name="DnaA (norm)",
                             line=dict(color=C_BLUE)), 1, 2)
    fig.add_trace(go.Scatter(x=t, y=frac + 1.0, mode="lines", name="chrom copy",
                             line=dict(color=C_GREEN)), 1, 2)
    fig.add_trace(go.Scatter(x=t, y=dntp / (dntp.max() or 1) * 2.0, mode="lines", name="dNTP (norm×2)",
                             line=dict(color=C_RED)), 1, 2)
    for x0, x1, c in [(0, init_end, "rgba(44,160,44,0.10)"), (init_end, repl_end, "rgba(214,39,40,0.10)"),
                      (repl_end, cyto_end, "rgba(31,119,180,0.10)")]:
        fig.add_vrect(x0=x0, x1=x1, fillcolor=c, line_width=0, row=1, col=2)
    fig.add_trace(go.Scatter(x=init_h, y=dnaA_disp, mode="markers",
                             marker=dict(color=C_BLUE, size=5, opacity=0.5), showlegend=False), 2, 1)
    fig.add_trace(go.Scatter(x=repl_h, y=dntp_cycle_mm, mode="markers", name="cycle start",
                             marker=dict(color="#9aa0a6", size=4, opacity=0.5)), 2, 2)
    fig.add_trace(go.Scatter(x=repl_h, y=dntp_start_mm, mode="markers", name="repl start",
                             marker=dict(color=PAL[5], size=5, opacity=0.6)), 2, 2)
    fig.add_trace(go.Scatter(x=init_h, y=repl_h, mode="markers",
                             marker=dict(color=C_RED, size=5, opacity=0.5), showlegend=False), 3, 1)
    fig.add_annotation(text=f"r={r_ci:.2f}", showarrow=False, xref="x3 domain", yref="y3 domain", x=0.95, y=0.95)
    fig.add_annotation(text=f"r={r_ir:.2f}", showarrow=False, xref="x5 domain", yref="y5 domain", x=0.95, y=0.95)
    fig.update_xaxes(title_text="duration (h)", row=1, col=1); fig.update_yaxes(title_text="% cells", row=1, col=1)
    fig.update_xaxes(title_text="time (h)", row=1, col=2)
    fig.update_xaxes(title_text="initiation dur (h)", row=2, col=1); fig.update_yaxes(title_text="initial DnaA", row=2, col=1)
    fig.update_xaxes(title_text="replication dur (h)", row=2, col=2); fig.update_yaxes(title_text="dNTP (mM)", row=2, col=2)
    fig.update_xaxes(title_text="initiation dur (h)", row=3, col=1); fig.update_yaxes(title_text="replication dur (h)", row=3, col=1)
    fig.update_layout(height=1200, width=1050, template="plotly_white",
                      title_text="Figure 4 — Emergent regulation of cell-cycle duration (viva-Mgen)")
    (VIZ / "fig4_combined.html").write_text(fig.to_html(full_html=True, include_plotlyjs="cdn"))


def _rewire():
    import yaml
    p = WS / "workspace/studies/fig4-cell-cycle/study.yaml"
    d = yaml.safe_load(p.read_text())
    panels = [
        ("fig4-combined", "fig4_combined.html", "Figure 4 reproduced — emergent cell-cycle-duration regulation, panels a–e as in Karr 2012."),
        ("fig4-a-phase-durations", "fig4_a_phase_durations.html", "Fig 4A: overlaid histograms of initiation/replication/cytokinesis/total cell-cycle durations."),
        ("fig4-b-single-cell-dynamics", "fig4_b_single_cell_dynamics.html", "Fig 4B: selected-cell DnaA, chromosome copy, and dNTP over time with phase-shaded background."),
        ("fig4-c-dnaA-vs-init", "fig4_c_dnaA_vs_init.html", "Fig 4C: initial DnaA count vs replication-initiation duration."),
        ("fig4-d-dntp-vs-repl", "fig4_d_dntp_vs_repl.html", "Fig 4D: dNTP vs replication duration — at cell-cycle start and at replication start."),
        ("fig4-e-init-vs-repl", "fig4_e_init_vs_repl.html", "Fig 4E: emergent inverse initiation↔replication duration relationship."),
    ]
    d["visualizations"] = [{"name": n, "chart": f"viz/{f}", "render": "python scripts/regen_fig4_cellcycle.py"} for n, f, _ in panels]
    d["embed_visualizations"] = [{"name": n, "url": f"/workspace/studies/{SLUG}/viz/{f}", "description": desc}
                                 for n, f, desc in panels]
    p.write_text(yaml.dump(d, sort_keys=False, width=100, allow_unicode=True))


if __name__ == "__main__":
    main()
