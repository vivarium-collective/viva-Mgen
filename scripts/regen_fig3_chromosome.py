#!/usr/bin/env python3
"""Regenerate Figure 3 for the fig3-expression study from the coordinate-resolved
ChromosomeDynamicsReproductionProcess — now REAL (not "not modeled"):

  (a) chromosome protein-occupancy density around the genome  (Fig 3A)
  (b) % chromosome explored over time, by protein             (Fig 3B)
  (c) % RNA expressed over time with t50                      (Fig 3C)
  (d) DNA-pol + RNA-pol position kymograph                    (Fig 3D)
  (e) protein-protein collision-frequency matrix              (Fig 3E)
  (f) collisions vs DNA-binding density                       (Fig 3F)

Writes panel HTML + fig3_combined.html into the study viz/ dir, and rewrites the
study's visualizations/embed_visualizations to this clean set via a YAML
round-trip (no regex — see the persist_results corruption history).
"""
from __future__ import annotations

import math
import pathlib

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

WS = pathlib.Path(__file__).resolve().parents[1]
VIZ = WS / "workspace/studies/fig3-expression/viz"
SLUG = "fig3-expression"
CDN = "https://cdn.plot.ly/plotly-3.7.0.min.js"


def _run():
    from viva_mgen.core import build_core
    from viva_mgen.processes.chromosome import ChromosomeDynamicsReproductionProcess
    from viva_mgen.processes.transcription import TranscriptionReproductionProcess
    from viva_mgen.expression_defaults import DEFAULT_GENES, synthesis_rates
    core = build_core()
    ch = ChromosomeDynamicsReproductionProcess(config={"seed": 0}, core=core)
    dt, T = 30.0, 2940.0  # sample every 30 s over 49 min (paper's DNA-pol window)
    times, exp_all, exp_rna, exp_dna, dens, ncoll = [], [], [], [], [], []
    rna_snaps, dna_snaps = [], []
    out = None
    for k in range(int(T / dt)):
        out = ch.update({"rna_polymerase": 120.0, "replication_active": 1.0}, dt)
        times.append(k * dt / 60.0)  # minutes
        exp_all.append(out["fraction_explored"] * 100)
        exp_rna.append(out["percent_rnap_explored"] * 100)
        exp_dna.append(out["percent_dnap_explored"] * 100)
        dens.append(out["dna_binding_density"])
        ncoll.append(out["n_collisions"])
        rna_snaps.append(list(out["rna_pol_positions"]))
        dna_snaps.append(list(out["dna_pol_positions"]))
    occ = out["occupancy"]; nb = ch.nb
    occ_arr = np.zeros(nb)
    for i, v in occ.items():
        occ_arr[int(i)] = v
    # per-interval collisions (delta) vs density for panel F
    coll_delta = np.diff([0] + ncoll)
    # collision matrix -> pair table
    pairs = out["collisions"]
    # % RNA expressed over time (fraction of genes ever transcribed)
    txn = TranscriptionReproductionProcess(config={"seed": 3}, core=core)
    seen = set(); frac_expr = []; tmin = []
    for k in range(2400):  # 40 min at 1 s
        d = txn.update({"ntp": 1e12, "rna_pol": 100.0}, 1.0)
        for g, n in d["rna_counts"].items():
            if n > 0:
                seen.add(g)
        if k % 20 == 0:
            frac_expr.append(100 * len(seen) / len(DEFAULT_GENES)); tmin.append(k / 60.0)
    return dict(times=times, exp_all=exp_all, exp_rna=exp_rna, exp_dna=exp_dna,
                dens=dens, coll_delta=coll_delta, occ=occ_arr, nb=nb,
                rna_snaps=rna_snaps, dna_snaps=dna_snaps, pairs=pairs,
                frac_expr=frac_expr, tmin=tmin, bp_per_bin=ch.bp_per_bin)


def _page(fig, title):
    from viva_mgen.viz import _page as vpage
    return vpage(fig, title)


PAL = None


def main():
    from viva_mgen import viz
    global PAL
    PAL = viz.PALETTE
    VIZ.mkdir(parents=True, exist_ok=True)
    d = _run()
    kb = d["bp_per_bin"] / 1000.0
    genome_kb = [i * kb for i in range(d["nb"])]

    # (a) occupancy density around the genome
    (VIZ / "fig3_a_occupancy.html").write_text(viz.line_series_html(
        "Chromosome protein-occupancy density (Fig 3A)", genome_kb,
        {"bound-protein occupancy": d["occ"]}, x_title="genome position (kb)",
        y_title="cumulative bound-protein occupancy"))

    # (b) % chromosome explored over time by protein
    (VIZ / "fig3_b_explored.html").write_text(viz.line_series_html(
        "Chromosome explored over time (Fig 3B)", d["times"],
        {"all proteins": d["exp_all"], "RNA polymerase": d["exp_rna"], "DNA polymerase": d["exp_dna"]},
        x_title="time (min)", y_title="% chromosome explored"))

    # (c) % RNA expressed with t50
    t50 = None
    for t, f in zip(d["tmin"], d["frac_expr"]):
        if f >= 50 and t50 is None:
            t50 = t
    figc = go.Figure(layout=viz._layout("% RNA expressed over time (Fig 3C)", "time (min)", "% genes expressed"))
    figc.add_trace(go.Scatter(x=d["tmin"], y=d["frac_expr"], mode="lines", line=dict(color=PAL[0], width=2), name="% expressed"))
    if t50:
        figc.add_vline(x=t50, line_dash="dash", line_color=PAL[1], annotation_text=f"t50 = {t50:.1f} min")
    (VIZ / "fig3_c_rna_expressed.html").write_text(viz._page(figc, "Fig 3C"))

    # (d) polymerase position kymograph (positions vs time)
    figd = go.Figure(layout=viz._layout("Polymerase positions over time (Fig 3D)", "time (min)", "genome position (kb)"))
    xr, yr = [], []
    for t, snap in zip(d["times"], d["rna_snaps"]):
        for p in snap:
            xr.append(t); yr.append(p * kb)
    figd.add_trace(go.Scattergl(x=xr, y=yr, mode="markers", name="RNA pol",
                                marker=dict(color=PAL[0], size=3, opacity=0.35)))
    for i in range(2):
        xs = d["times"]; ys = [(snap[i] * kb if len(snap) > i else None) for snap in d["dna_snaps"]]
        figd.add_trace(go.Scatter(x=xs, y=ys, mode="markers", name=f"DNA pol {i+1}",
                                  marker=dict(color=PAL[1] if i == 0 else PAL[3], size=5)))
    (VIZ / "fig3_d_polymerase_traces.html").write_text(viz._page(figd, "Fig 3D"))

    # (e) collision-frequency matrix (mover x occupant)
    movers, occs, mat = _collision_matrix(d["pairs"])
    fige = viz.heatmap_html("Protein-protein collision frequency (Fig 3E)", mat,
                            x=occs, y=movers, x_title="displaced protein",
                            y_title="displacing protein", colorbar_title="collisions")
    (VIZ / "fig3_e_collision_matrix.html").write_text(fige)

    # (f) collisions vs DNA-binding density
    (VIZ / "fig3_f_collisions_vs_density.html").write_text(viz.scatter_html(
        "Collisions vs DNA-binding density (Fig 3F)", d["dens"], d["coll_delta"],
        x_title="DNA-binding density (fraction of genome bound)",
        y_title="collisions per interval"))

    # combined 3x2 figure
    _combined(d, genome_kb, kb, t50, movers, occs, mat)

    _rewire_embeds()
    print("fig3 chromosome panels + combined regenerated; t50 =", t50)


def _collision_matrix(pairs):
    movers, occs = set(), set()
    for key in pairs:
        m, o = key.split("||")
        movers.add(m); occs.add(o)
    movers = sorted(movers); occs = sorted(occs)
    mat = [[pairs.get(f"{m}||{o}", 0.0) for o in occs] for m in movers]
    return movers, occs, mat


def _combined(d, genome_kb, kb, t50, movers, occs, mat):
    fig = make_subplots(rows=3, cols=2, subplot_titles=(
        "(a) chromosome occupancy density", "(b) % chromosome explored",
        "(c) % RNA expressed", "(d) polymerase positions",
        "(e) collision-frequency matrix", "(f) collisions vs binding density"),
        specs=[[{"type": "xy"}, {"type": "xy"}], [{"type": "xy"}, {"type": "xy"}],
               [{"type": "heatmap"}, {"type": "xy"}]], vertical_spacing=0.09)
    fig.add_trace(go.Scatter(x=genome_kb, y=d["occ"], mode="lines", line=dict(color=PAL[0]), showlegend=False), 1, 1)
    for i, (name, ys) in enumerate([("all", d["exp_all"]), ("RNA pol", d["exp_rna"]), ("DNA pol", d["exp_dna"])]):
        fig.add_trace(go.Scatter(x=d["times"], y=ys, mode="lines", name=name, line=dict(color=PAL[i])), 1, 2)
    fig.add_trace(go.Scatter(x=d["tmin"], y=d["frac_expr"], mode="lines", line=dict(color=PAL[0]), showlegend=False), 2, 1)
    xr, yr = [], []
    for t, snap in zip(d["times"], d["rna_snaps"]):
        for p in snap:
            xr.append(t); yr.append(p * kb)
    fig.add_trace(go.Scattergl(x=xr, y=yr, mode="markers", marker=dict(color=PAL[0], size=2, opacity=0.3), showlegend=False), 2, 2)
    for i in range(2):
        fig.add_trace(go.Scatter(x=d["times"], y=[(s[i] * kb if len(s) > i else None) for s in d["dna_snaps"]],
                                 mode="markers", marker=dict(color=PAL[1], size=4), showlegend=False), 2, 2)
    fig.add_trace(go.Heatmap(z=mat, x=occs, y=movers, colorscale="Blues", showscale=False), 3, 1)
    fig.add_trace(go.Scatter(x=d["dens"], y=d["coll_delta"], mode="markers",
                             marker=dict(color=PAL[2], size=6, opacity=0.7), showlegend=False), 3, 2)
    fig.update_layout(title_text="Figure 3 — Chromosome DNA–protein interactions (viva-Mgen)",
                      paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                      font=dict(family="system-ui, sans-serif", size=11), height=1050, showlegend=True,
                      legend=dict(orientation="h", y=-0.05))
    div = "viz-fig3c"
    html = (f'<!doctype html><html><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<title>Figure 3 (viva-Mgen)</title><script src="{CDN}"></script>'
            f'<style>html,body{{margin:0;background:transparent}}#{div}{{width:100%;height:100vh;min-height:900px}}</style>'
            f'</head><body><div id="{div}"></div><script>var s={fig.to_json()};'
            f"Plotly.newPlot('{div}',s.data,s.layout,{{responsive:true,displayModeBar:false}});</script></body></html>")
    (VIZ / "fig3_combined.html").write_text(html)


def _rewire_embeds():
    import yaml
    p = WS / "workspace/studies/fig3-expression/study.yaml"
    d = yaml.safe_load(p.read_text())
    panels = [("fig3-combined", "fig3_combined.html", "Figure 3 reproduced — chromosome DNA–protein interactions, all subpanels (a–f) as in Karr 2012."),
              ("fig3-a-occupancy", "fig3_a_occupancy.html", "Fig 3A: chromosome protein-occupancy density around the genome."),
              ("fig3-b-explored", "fig3_b_explored.html", "Fig 3B: % chromosome explored over time, by protein."),
              ("fig3-c-rna-expressed", "fig3_c_rna_expressed.html", "Fig 3C: % RNA expressed over time with t50."),
              ("fig3-d-polymerase-traces", "fig3_d_polymerase_traces.html", "Fig 3D: RNA-pol + DNA-pol position kymograph."),
              ("fig3-e-collision-matrix", "fig3_e_collision_matrix.html", "Fig 3E: protein-protein collision-frequency matrix."),
              ("fig3-f-collisions-vs-density", "fig3_f_collisions_vs_density.html", "Fig 3F: collisions vs DNA-binding density.")]
    d["visualizations"] = [{"name": n, "chart": f"viz/{f}", "render": "python sims/run.py"} for n, f, _ in panels]
    d["embed_visualizations"] = [{"name": n, "url": f"/workspace/studies/{SLUG}/viz/{f}", "description": desc}
                                 for n, f, desc in panels]
    p.write_text(yaml.dump(d, sort_keys=False, width=100, allow_unicode=True))


if __name__ == "__main__":
    raise SystemExit(main())
