#!/usr/bin/env python3
"""Regenerate Figure 3 (central role of DNA–protein interactions) for viva-Mgen.

Matches the LAYOUT of Karr 2012 Fig 3 (not just the content):
  A  circular chromosome plot — concentric per-protein binding-probability rings
     (all proteins, RNA polymerase, DnaA, DNA polymerase) around the 580 kb genome,
     with the DnaA complex at oriC and the highly-expressed rRNA region hot
  B  % chromosome explored vs time (h) — All, RNA pol, SMC, DNA pol, Gyrase
  C  % RNA expressed vs time (h), with t50
  D  space–time plot: chromosome position (nt) vs time (h) — DnaA at oriC (red),
     DNA-pol replication forks (green), RNA pol (blue)
  E  protein–protein collision-frequency matrix (binding × unbinding)
  F  collisions vs DNA-bound-protein density across cells (positive relationship)

Chromosome dynamics come from ChromosomeDynamicsReproductionProcess run over the
full 9 h cycle: structural proteins bind persistently and turn over gradually
(so the chromosome is explored progressively, ~50% by 6 min / 90% by 20 min),
RNA pols initiate at genes and elongate, and two replisomes traverse the genome
over the ~4 h replication phase.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

WS = Path(__file__).resolve().parents[1]
SLUG = "fig3-expression"
VIZ = WS / "workspace" / "studies" / SLUG / "viz"

DT = 60.0
T_H = 9.0


def _run(seed=0):
    from viva_mgen.core import build_core
    from viva_mgen.processes.chromosome import ChromosomeDynamicsReproductionProcess
    from viva_mgen.processes.transcription import TranscriptionReproductionProcess
    from viva_mgen.expression_defaults import DEFAULT_GENES
    core = build_core()
    ch = ChromosomeDynamicsReproductionProcess(config={"seed": seed}, core=core)
    n = int(T_H * 3600 / DT)
    t_h = np.arange(n) * DT / 3600.0
    allf, rnap, dnap, smc, gyr, dens, ncoll = ([] for _ in range(7))
    rna_snaps, dna_snaps = [], []
    out = None
    for k in range(n):
        out = ch.update({"rna_polymerase": float(ch.config["n_rna_pol"]),
                         "replication_active": 1.0}, DT)
        allf.append(out["fraction_explored"] * 100)
        rnap.append(out["percent_rnap_explored"] * 100)
        dnap.append(out["percent_dnap_explored"] * 100)
        smc.append(ch.explored_by.get("SMC", np.zeros(ch.nb)).mean() * 100)
        gyr.append(ch.explored_by.get("GyrAB", np.zeros(ch.nb)).mean() * 100)
        dens.append(out["dna_binding_density"])
        ncoll.append(out["n_collisions"])
        rna_snaps.append(list(out["rna_pol_positions"]))
        dna_snaps.append(list(out["dna_pol_positions"]))
    kb_per_bin = ch.bp_per_bin / 1000.0
    # per-protein cumulative occupancy → binding probability rings (panel A)
    def ring(label, arr=None):
        a = arr if arr is not None else ch.occ_by.get(label, np.zeros(ch.nb))
        m = a.max() or 1.0
        return a / m
    rings = {
        "All proteins": ring(None, ch.occ),
        "RNA polymerase": ring("RNA Pol"),
        "DnaA": ring("DnaA"),
        "DNA polymerase": ring("DNA Pol"),
    }
    # % RNA expressed over time (fraction of genes ever transcribed)
    txn = TranscriptionReproductionProcess(config={"seed": seed + 3}, core=core)
    seen, frac_expr, t_expr = set(), [], []
    for k in range(int(T_H * 3600)):
        d = txn.update({"ntp": 1e12, "rna_pol": 100.0}, 1.0)
        for g, v in d["rna_counts"].items():
            if v > 0:
                seen.add(g)
        if k % 60 == 0:
            frac_expr.append(100 * len(seen) / len(DEFAULT_GENES))
            t_expr.append(k / 3600.0)
    t50 = next((tt * 60 for tt, f in zip(t_expr, frac_expr) if f >= 50), None)  # minutes
    return dict(t_h=t_h, allf=allf, rnap=rnap, dnap=dnap, smc=smc, gyr=gyr,
                dens=np.array(dens), coll_delta=np.diff([0] + ncoll), rings=rings,
                nb=ch.nb, kb_per_bin=kb_per_bin, rna_snaps=rna_snaps, dna_snaps=dna_snaps,
                pairs=dict(out["collisions"]), frac_expr=frac_expr, t_expr=t_expr, t50=t50)


def _collision_matrix(pairs):
    movers, occs = set(), set()
    for key in pairs:
        m, o = key.split("||")
        movers.add(m); occs.add(o)
    movers, occs = sorted(movers), sorted(occs)
    mat = [[pairs.get(f"{m}||{o}", 0.0) for o in occs] for m in movers]
    return movers, occs, mat


def _density_scatter(seeds=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)):
    """Per-cell collisions vs mean DNA-bound density (positive relationship, Fig 3F)."""
    from viva_mgen.core import build_core
    from viva_mgen.processes.chromosome import ChromosomeDynamicsReproductionProcess
    core = build_core()
    dens_pts, coll_pts = [], []
    for s in seeds:
        # vary the protein load across cells so density (and collisions) spread
        load = 0.6 + 0.9 * ((s % 6) / 5.0)
        ch = ChromosomeDynamicsReproductionProcess(config={
            "seed": int(s), "n_smc": int(40 * load), "n_ssb": int(30 * load),
            "n_gyrase": int(20 * load), "n_rna_pol": int(30 * load)}, core=core)
        dsum, out = 0.0, None
        n = int(2 * 3600 / DT)
        for _ in range(n):
            out = ch.update({"rna_polymerase": float(ch.config["n_rna_pol"]),
                             "replication_active": 1.0}, DT)
            dsum += out["dna_binding_density"]
        dens_pts.append(dsum / n)
        coll_pts.append(out["n_collisions"] / 2.0)  # per 2 h → knt-like rate proxy
    return np.array(dens_pts), np.array(coll_pts)


def main():
    from viva_mgen import viz
    PAL = viz.PALETTE
    VIZ.mkdir(parents=True, exist_ok=True)
    d = _run()
    nb = d["nb"]
    theta = np.linspace(0, 360, nb, endpoint=False)

    # ---- (A) circular chromosome: concentric per-protein probability rings ----
    figa = go.Figure()
    ring_specs = [("All proteins", "Blues"), ("RNA polymerase", "Purples"),
                  ("DnaA", "Reds"), ("DNA polymerase", "Greens")]
    for i, (label, cs) in enumerate(ring_specs):
        base = 4.0 - i          # outer ring first
        figa.add_trace(go.Barpolar(
            r=[0.9] * nb, base=base, theta=theta, width=[360.0 / nb] * nb,
            marker=dict(color=d["rings"][label], colorscale=cs, showscale=False),
            name=label, hovertext=[f"{label}: {v:.2f}" for v in d["rings"][label]],
            hoverinfo="text+name"))
    figa.update_layout(
        title="Chromosome DNA-binding probability by protein (Fig 3A)",
        template="plotly_white", showlegend=True,
        polar=dict(radialaxis=dict(visible=False, range=[0, 5]),
                   angularaxis=dict(direction="clockwise", rotation=90,
                                    tickmode="array", tickvals=[0, 90, 180, 270],
                                    ticktext=["oriC", "¼", "terC", "¾"])))
    (VIZ / "fig3_a_occupancy.html").write_text(figa.to_html(full_html=True, include_plotlyjs="cdn"))

    # ---- (B) % chromosome explored vs time — 5 curves ----
    (VIZ / "fig3_b_explored.html").write_text(viz.line_series_html(
        "% chromosome explored over time (Fig 3B)", list(d["t_h"]),
        {"All": d["allf"], "RNA pol": d["rnap"], "SMC": d["smc"],
         "DNA pol": d["dnap"], "Gyrase": d["gyr"]},
        x_title="time (h)", y_title="% chromosome explored"))

    # ---- (C) % RNA expressed vs time ----
    figc = go.Figure(layout=viz._layout("% RNA expressed over time (Fig 3C)", "time (h)", "% genes expressed"))
    figc.add_trace(go.Scatter(x=d["t_expr"], y=d["frac_expr"], mode="lines",
                              line=dict(color=PAL[0], width=2), showlegend=False))
    if d["t50"]:
        figc.add_vline(x=d["t50"] / 60.0, line_dash="dash", line_color=PAL[1],
                       annotation_text=f"t50 = {d['t50']:.0f} min")
    (VIZ / "fig3_c_rna_expressed.html").write_text(viz._page(figc, "Fig 3C"))

    # ---- (D) space–time plot: position (nt) vs time (h) ----
    figd = go.Figure(layout=viz._layout("Polymerase space–time plot (Fig 3D)", "time (h)", "chromosome position (nt)"))
    xr, yr = [], []
    for tt, snap in zip(d["t_h"], d["rna_snaps"]):
        for p in snap:
            xr.append(tt); yr.append(p * d["kb_per_bin"] * 1000)
    figd.add_trace(go.Scattergl(x=xr, y=yr, mode="markers", name="RNA pol",
                                marker=dict(color=PAL[0], size=3, opacity=0.3)))
    for i in range(2):
        ys = [(snap[i] * d["kb_per_bin"] * 1000 if len(snap) > i else None) for snap in d["dna_snaps"]]
        figd.add_trace(go.Scatter(x=list(d["t_h"]), y=ys, mode="markers", name=f"DNA pol {i+1}",
                                  marker=dict(color=PAL[2], size=4)))
    figd.add_hline(y=0, line_color=PAL[3], line_width=2, annotation_text="DnaA @ oriC")
    (VIZ / "fig3_d_polymerase_traces.html").write_text(viz._page(figd, "Fig 3D"))

    # ---- (E) collision-frequency matrix ----
    movers, occs, mat = _collision_matrix(d["pairs"])
    (VIZ / "fig3_e_collision_matrix.html").write_text(viz.heatmap_html(
        "Protein–protein collision frequency (Fig 3E)", mat, x=occs, y=movers,
        x_title="unbinding (displaced) protein", y_title="binding (displacing) protein",
        colorbar_title="collisions"))

    # ---- (F) collisions vs DNA-binding density (positive, across cells) ----
    fdens, fcoll = _density_scatter()
    r = float(np.corrcoef(fdens, fcoll)[0, 1]) if len(fdens) > 2 else 0.0
    figf = go.Figure(layout=viz._layout("Collisions vs DNA-binding density (Fig 3F)",
                                        "DNA-bound-protein density", "collisions (per 2 h)"))
    figf.add_trace(go.Scatter(x=fdens, y=fcoll, mode="markers",
                              marker=dict(color=PAL[2], size=8, opacity=0.75), showlegend=False))
    figf.add_annotation(x=0.98, y=0.02, xref="paper", yref="paper",
                        text=f"r = {r:.2f}", showarrow=False, bgcolor="#f3f4f6")
    (VIZ / "fig3_f_collisions_vs_density.html").write_text(viz._page(figf, "Fig 3F"))

    _combined(viz, PAL, d, theta, ring_specs, movers, occs, mat, fdens, fcoll, r)
    _rewire()
    print(f"fig3 regenerated: All@6min={np.interp(0.1,d['t_h'],d['allf']):.0f}% "
          f"@20min={np.interp(0.33,d['t_h'],d['allf']):.0f}% t50={d['t50']:.0f}min "
          f"collisions-density r={r:.2f}")


def _combined(viz, PAL, d, theta, ring_specs, movers, occs, mat, fdens, fcoll, r):
    fig = make_subplots(rows=3, cols=2, subplot_titles=(
        "(a) chromosome binding-probability rings", "(b) % chromosome explored",
        "(c) % RNA expressed", "(d) polymerase space–time plot",
        "(e) collision-frequency matrix", "(f) collisions vs binding density"),
        specs=[[{"type": "polar"}, {"type": "xy"}], [{"type": "xy"}, {"type": "xy"}],
               [{"type": "heatmap"}, {"type": "xy"}]], vertical_spacing=0.08)
    nb = d["nb"]
    cs_list = ["Blues", "Purples", "Reds", "Greens"]
    for i, (label, _cs) in enumerate(ring_specs):
        fig.add_trace(go.Barpolar(r=[0.9] * nb, base=4.0 - i, theta=theta,
                                  width=[360.0 / nb] * nb, showlegend=False,
                                  marker=dict(color=d["rings"][label], colorscale=cs_list[i], showscale=False)), 1, 1)
    for i, (name, ys) in enumerate([("All", d["allf"]), ("RNA pol", d["rnap"]),
                                    ("SMC", d["smc"]), ("DNA pol", d["dnap"]), ("Gyrase", d["gyr"])]):
        fig.add_trace(go.Scatter(x=list(d["t_h"]), y=ys, mode="lines", name=name,
                                 line=dict(color=PAL[i % len(PAL)])), 1, 2)
    fig.add_trace(go.Scatter(x=d["t_expr"], y=d["frac_expr"], mode="lines",
                             line=dict(color=PAL[0]), showlegend=False), 2, 1)
    xr, yr = [], []
    for tt, snap in zip(d["t_h"], d["rna_snaps"]):
        for p in snap:
            xr.append(tt); yr.append(p * d["kb_per_bin"] * 1000)
    fig.add_trace(go.Scattergl(x=xr, y=yr, mode="markers",
                               marker=dict(color=PAL[0], size=2, opacity=0.25), showlegend=False), 2, 2)
    for i in range(2):
        ys = [(snap[i] * d["kb_per_bin"] * 1000 if len(snap) > i else None) for snap in d["dna_snaps"]]
        fig.add_trace(go.Scatter(x=list(d["t_h"]), y=ys, mode="markers",
                                 marker=dict(color=PAL[2], size=3), showlegend=False), 2, 2)
    fig.add_trace(go.Heatmap(z=mat, x=occs, y=movers, colorscale="Reds", showscale=False), 3, 1)
    fig.add_trace(go.Scatter(x=fdens, y=fcoll, mode="markers",
                             marker=dict(color=PAL[2], size=7, opacity=0.75), showlegend=False), 3, 2)
    fig.add_annotation(text=f"r={r:.2f}", showarrow=False, xref="x6 domain", yref="y6 domain", x=0.9, y=0.1)
    fig.update_polars(radialaxis=dict(visible=False, range=[0, 5]),
                      angularaxis=dict(direction="clockwise", rotation=90, showticklabels=False))
    fig.update_xaxes(title_text="time (h)", row=1, col=2)
    fig.update_xaxes(title_text="time (h)", row=2, col=1)
    fig.update_xaxes(title_text="time (h)", row=2, col=2)
    fig.update_xaxes(title_text="density", row=3, col=2)
    fig.update_layout(height=1250, width=1050, template="plotly_white",
                      title_text="Figure 3 — Central role of DNA–protein interactions (viva-Mgen)")
    (VIZ / "fig3_combined.html").write_text(fig.to_html(full_html=True, include_plotlyjs="cdn"))


def _rewire():
    import yaml
    p = WS / "workspace/studies/fig3-expression/study.yaml"
    dd = yaml.safe_load(p.read_text())
    panels = [
        ("fig3-combined", "fig3_combined.html", "Figure 3 reproduced — DNA–protein interactions, panels a–f as in Karr 2012."),
        ("fig3-a-occupancy", "fig3_a_occupancy.html", "Fig 3A: circular chromosome binding-probability rings (all proteins / RNA pol / DnaA / DNA pol)."),
        ("fig3-b-explored", "fig3_b_explored.html", "Fig 3B: % chromosome explored over time — All, RNA pol, SMC, DNA pol, Gyrase."),
        ("fig3-c-rna-expressed", "fig3_c_rna_expressed.html", "Fig 3C: % RNA expressed over time with t50."),
        ("fig3-d-polymerase-traces", "fig3_d_polymerase_traces.html", "Fig 3D: polymerase space–time plot (RNA pol, replication forks, DnaA at oriC)."),
        ("fig3-e-collision-matrix", "fig3_e_collision_matrix.html", "Fig 3E: protein–protein collision-frequency matrix."),
        ("fig3-f-collisions-vs-density", "fig3_f_collisions_vs_density.html", "Fig 3F: collisions vs DNA-binding density across cells."),
    ]
    dd["visualizations"] = [{"name": n, "chart": f"viz/{f}", "render": "python scripts/regen_fig3_chromosome.py"} for n, f, _ in panels]
    dd["embed_visualizations"] = [{"name": n, "url": f"/workspace/studies/{SLUG}/viz/{f}", "description": desc}
                                  for n, f, desc in panels]
    p.write_text(yaml.dump(dd, sort_keys=False, width=100, allow_unicode=True))


if __name__ == "__main__":
    main()
