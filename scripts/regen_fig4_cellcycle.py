#!/usr/bin/env python3
"""Regenerate Figure 4 (emergent cell-cycle-duration regulation) for viva-Mgen.

Matches Karr 2012 Fig 4: a 128-cell population of single cells is simulated
through replication initiation → replication → cytokinesis. The cell cycle is
now ~9 h (initiation ~3.2 h, replication ~3.8 h, cytokinesis ~1.1 h) — the
replication phase is dNTP-synthesis-limited (see ReplicationReproductionProcess.
init_synth_fraction) rather than polymerase-limited, so it lands near the paper's
4.33 h instead of ~1 h, while the emergent inverse initiation↔replication
relationship (the paper's headline result) is preserved as cell-to-cell
variation.

Writes fig4_a..e + fig4_combined into the study viz/ dir and rewires the study's
visualizations/embed_visualizations via a YAML round-trip (never a regex edit).
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
CYTO_MEAN_S = 1.08 * 3600.0   # paper's fixed cytokinesis duration
CYTO_CV = 0.044               # paper's cytokinesis CV (4.4%)


def _run_population(core):
    from viva_mgen.processes.replication import ReplicationReproductionProcess
    init_dur, repl_dur, dntp_start, dnaA0, cyto = [], [], [], [], []
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
        init_dur.append(out["initiation_duration"])
        repl_dur.append(out["replication_duration"])
        dntp_start.append(out["dntp_at_replication_start"])
        dnaA0.append(a0)
        cyto.append(CYTO_MEAN_S * (1.0 + CYTO_CV * float(rng.standard_normal())))
    return (np.array(init_dur), np.array(repl_dur), np.array(dntp_start),
            np.array(dnaA0), np.array(cyto))


def _trajectory(core, a0=4.0, d0=6000.0, seed=1):
    from viva_mgen.processes.replication import ReplicationReproductionProcess
    p = ReplicationReproductionProcess(
        {"initial_dnaA": a0, "initial_dntp": d0, "seed": seed}, core=core)
    t, frac, dntp, dnaA, phase = [], [], [], [], []
    for k in range(MAX_STEPS):
        out = p.update({"dntp_synthesis_scale": 1.0}, DT)
        t.append(k * DT / 3600.0)
        frac.append(out["replicated_fraction"])
        dntp.append(out["dntp_pool"])
        dnaA.append(out["dnaA_complex"])
        phase.append(out["phase_code"])
        if out["phase_code"] == _DONE and k > 5:
            # run a little past completion for context
            if k * DT / 3600.0 > (out["initiation_duration"] + out["replication_duration"]) / 3600.0 + 0.5:
                break
    return (np.array(t), np.array(frac), np.array(dntp), np.array(dnaA), np.array(phase))


def _fit(x, y):
    m, b = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    return xs, m * xs + b, float(np.corrcoef(x, y)[0, 1])


def main():
    from viva_mgen import viz
    from viva_mgen.core import build_core
    PAL = viz.PALETTE
    VIZ.mkdir(parents=True, exist_ok=True)
    core = build_core()

    idur, rdur, dstart, dnaA0, cyto = _run_population(core)
    init_h, repl_h, cyto_h = idur / 3600.0, rdur / 3600.0, cyto / 3600.0
    total_h = init_h + repl_h + cyto_h
    n = len(idur)

    # (a) phase-duration distributions (violins)
    figa = go.Figure(layout=viz._layout("Phase-duration distributions (Fig 4A)",
                                         "", "duration (h)"))
    for name, data, col in [("initiation", init_h, PAL[0]), ("replication", repl_h, PAL[3]),
                            ("cytokinesis", cyto_h, PAL[2]), ("cell cycle", total_h, PAL[1])]:
        figa.add_trace(go.Violin(y=data, name=f"{name}\nCV={data.std()/data.mean()*100:.0f}%",
                                 line_color=col, meanline_visible=True, points=False))
    figa.update_layout(showlegend=False)
    (VIZ / "fig4_a_phase_durations.html").write_text(viz._page(figa, "Fig 4A"))

    # (b) single-cell trajectory with phase shading
    t, frac, dntp, dnaA, phase = _trajectory(core)
    figb = go.Figure(layout=viz._layout("Single-cell cell-cycle dynamics (Fig 4B)",
                                        "time (h)", "level"))
    figb.add_trace(go.Scatter(x=t, y=dnaA / (dnaA.max() or 1), mode="lines",
                              name="DnaA complex", line=dict(color=PAL[0])))
    figb.add_trace(go.Scatter(x=t, y=frac + 1.0, mode="lines",
                              name="chromosome copy (1→2)", line=dict(color=PAL[2])))
    figb.add_trace(go.Scatter(x=t, y=dntp / (dntp.max() or 1), mode="lines",
                              name="dNTP pool (norm)", line=dict(color=PAL[3])))
    # shade phases
    init_end = idur.mean() / 3600.0
    repl_end = init_end + rdur.mean() / 3600.0
    for x0, x1, c in [(0, init_end, "rgba(31,119,180,0.07)"),
                      (init_end, repl_end, "rgba(214,39,40,0.07)"),
                      (repl_end, repl_end + cyto_h.mean(), "rgba(44,160,44,0.07)")]:
        figb.add_vrect(x0=x0, x1=x1, fillcolor=c, line_width=0)
    (VIZ / "fig4_b_single_cell_dynamics.html").write_text(viz._page(figb, "Fig 4B"))

    # (c) initial DnaA vs initiation duration
    xs, ys, r_ci = _fit(dnaA0, init_h)
    figc = go.Figure(layout=viz._layout("Initial DnaA vs initiation duration (Fig 4C)",
                                        "initial DnaA (molecules)", "initiation dur. (h)"))
    figc.add_trace(go.Scatter(x=dnaA0, y=init_h, mode="markers",
                              marker=dict(color=PAL[0], size=6, opacity=0.6), showlegend=False))
    figc.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=PAL[0]), showlegend=False))
    figc.add_annotation(x=0.98, y=0.98, xref="paper", yref="paper",
                        text=f"r = {r_ci:.2f}", showarrow=False, bgcolor="#f3f4f6")
    (VIZ / "fig4_c_dnaA_vs_init.html").write_text(viz._page(figc, "Fig 4C"))

    # (d) dNTP at replication start vs replication duration
    xs, ys, r_dr = _fit(dstart, repl_h)
    figd = go.Figure(layout=viz._layout("dNTP surplus vs replication duration (Fig 4D)",
                                        "dNTP at repl. start (molecules)", "replication dur. (h)"))
    figd.add_trace(go.Scatter(x=dstart, y=repl_h, mode="markers",
                              marker=dict(color=PAL[3], size=6, opacity=0.6), showlegend=False))
    figd.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=PAL[3]), showlegend=False))
    figd.add_annotation(x=0.98, y=0.98, xref="paper", yref="paper",
                        text=f"r = {r_dr:.2f}", showarrow=False, bgcolor="#f3f4f6")
    (VIZ / "fig4_d_dntp_vs_repl.html").write_text(viz._page(figd, "Fig 4D"))

    # (e) initiation vs replication duration — the emergent inverse relationship
    xs, ys, r_ir = _fit(init_h, repl_h)
    fige = go.Figure(layout=viz._layout("Emergent inverse initiation↔replication (Fig 4E)",
                                        "initiation dur. (h)", "replication dur. (h)"))
    fige.add_trace(go.Scatter(x=init_h, y=repl_h, mode="markers",
                              marker=dict(color=PAL[5], size=6, opacity=0.6), showlegend=False))
    fige.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=PAL[5]), showlegend=False))
    fige.add_annotation(x=0.98, y=0.98, xref="paper", yref="paper",
                        text=f"r = {r_ir:.2f}", showarrow=False, bgcolor="#f3f4f6")
    (VIZ / "fig4_e_init_vs_repl.html").write_text(viz._page(fige, "Fig 4E"))

    _combined(PAL, viz, init_h, repl_h, cyto_h, total_h, t, frac, dntp, dnaA,
              init_end, repl_end, dnaA0, dstart, r_ci, r_dr, r_ir)
    _rewire_embeds()
    print(f"fig4 regenerated: n={n} | init={init_h.mean():.2f}h repl={repl_h.mean():.2f}h "
          f"cyto={cyto_h.mean():.2f}h TOTAL={total_h.mean():.2f}h (CV {total_h.std()/total_h.mean()*100:.1f}%) | "
          f"r(init,repl)={r_ir:.2f} r(dntp,repl)={r_dr:.2f} r(dnaA,init)={r_ci:.2f}")


def _combined(PAL, viz, init_h, repl_h, cyto_h, total_h, t, frac, dntp, dnaA,
              init_end, repl_end, dnaA0, dstart, r_ci, r_dr, r_ir):
    fig = make_subplots(rows=3, cols=2, subplot_titles=(
        "(a) phase-duration distributions", "(b) single-cell dynamics",
        "(c) initial DnaA vs initiation dur.", "(d) dNTP surplus vs replication dur.",
        "(e) emergent inverse initiation↔replication", ""),
        vertical_spacing=0.09)
    for name, data, col in [("initiation", init_h, PAL[0]), ("replication", repl_h, PAL[3]),
                            ("cytokinesis", cyto_h, PAL[2]), ("cell cycle", total_h, PAL[1])]:
        fig.add_trace(go.Violin(y=data, name=name, line_color=col,
                                meanline_visible=True, points=False, showlegend=False), 1, 1)
    fig.add_trace(go.Scatter(x=t, y=dnaA/(dnaA.max() or 1), mode="lines", name="DnaA",
                             line=dict(color=PAL[0])), 1, 2)
    fig.add_trace(go.Scatter(x=t, y=frac+1.0, mode="lines", name="chrom copy",
                             line=dict(color=PAL[2])), 1, 2)
    fig.add_trace(go.Scatter(x=t, y=dntp/(dntp.max() or 1), mode="lines", name="dNTP",
                             line=dict(color=PAL[3])), 1, 2)
    fig.add_trace(go.Scatter(x=dnaA0, y=init_h, mode="markers",
                             marker=dict(color=PAL[0], size=5, opacity=0.5), showlegend=False), 2, 1)
    fig.add_trace(go.Scatter(x=dstart, y=repl_h, mode="markers",
                             marker=dict(color=PAL[3], size=5, opacity=0.5), showlegend=False), 2, 2)
    fig.add_trace(go.Scatter(x=init_h, y=repl_h, mode="markers",
                             marker=dict(color=PAL[5], size=5, opacity=0.5), showlegend=False), 3, 1)
    for (rr, cc, txt) in [(2, 1, f"r={r_ci:.2f}"), (2, 2, f"r={r_dr:.2f}"), (3, 1, f"r={r_ir:.2f}")]:
        fig.add_annotation(text=txt, showarrow=False, xref=f"x{'' if (rr,cc)==(1,1) else (rr-1)*2+cc} domain",
                           yref=f"y{'' if (rr,cc)==(1,1) else (rr-1)*2+cc} domain", x=0.95, y=0.95)
    fig.update_yaxes(title_text="duration (h)", row=1, col=1)
    fig.update_xaxes(title_text="time (h)", row=1, col=2)
    fig.update_xaxes(title_text="initial DnaA", row=2, col=1); fig.update_yaxes(title_text="init dur (h)", row=2, col=1)
    fig.update_xaxes(title_text="dNTP at start", row=2, col=2); fig.update_yaxes(title_text="repl dur (h)", row=2, col=2)
    fig.update_xaxes(title_text="initiation dur (h)", row=3, col=1); fig.update_yaxes(title_text="replication dur (h)", row=3, col=1)
    fig.update_layout(height=1200, width=1000, template="plotly_white",
                      title_text="Figure 4 — Emergent regulation of cell-cycle duration (viva-Mgen)")
    (VIZ / "fig4_combined.html").write_text(fig.to_html(full_html=True, include_plotlyjs="cdn"))


def _rewire_embeds():
    import yaml
    p = WS / "workspace/studies/fig4-cell-cycle/study.yaml"
    d = yaml.safe_load(p.read_text())
    panels = [
        ("fig4-combined", "fig4_combined.html", "Figure 4 reproduced — emergent cell-cycle-duration regulation, all subpanels (a–e) as in Karr 2012."),
        ("fig4-a-phase-durations", "fig4_a_phase_durations.html", "Fig 4A: initiation/replication/cytokinesis/total cell-cycle duration distributions (~9 h cycle)."),
        ("fig4-b-single-cell-dynamics", "fig4_b_single_cell_dynamics.html", "Fig 4B: single-cell DnaA, chromosome copy, and dNTP dynamics with phase shading."),
        ("fig4-c-dnaA-vs-init", "fig4_c_dnaA_vs_init.html", "Fig 4C: initial DnaA vs initiation duration (negative)."),
        ("fig4-d-dntp-vs-repl", "fig4_d_dntp_vs_repl.html", "Fig 4D: dNTP surplus at replication start vs replication duration (negative)."),
        ("fig4-e-init-vs-repl", "fig4_e_init_vs_repl.html", "Fig 4E: emergent inverse initiation↔replication duration relationship."),
    ]
    d["visualizations"] = [{"name": n, "chart": f"viz/{f}", "render": "python scripts/regen_fig4_cellcycle.py"} for n, f, _ in panels]
    d["embed_visualizations"] = [{"name": n, "url": f"/workspace/studies/{SLUG}/viz/{f}", "description": desc}
                                 for n, f, desc in panels]
    p.write_text(yaml.dump(d, sort_keys=False, width=100, allow_unicode=True))


if __name__ == "__main__":
    main()
