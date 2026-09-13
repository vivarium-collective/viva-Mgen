#!/usr/bin/env python3
"""Regenerate the ParCa (parameter-fitting) study figures for viva-Mgen.

The parameter calculator reproduces Karr 2012 FitConstants. These panels show what
it produces from the real observed knowledge-base data:
  A  fitted synthesis rate vs observed expression (per gene, log-log) — the rates
     are derived from real expression, coloured by RNA type;
  B  per-gene mRNA half-life distribution + the stable-RNA half-lives (by type);
  C  metabolic demand — the NMP composition the transcriptome implies (AT-rich);
  D  expression↔metabolism closure — the precursor supply fluxes the iPS189 network
     provides for the fitted composition.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

WS = Path(__file__).resolve().parents[1]
SLUG = "parca-parameter-fitting"
VIZ = WS / "workspace" / "studies" / SLUG / "viz"


def _data():
    from viva_mgen import parca
    from viva_mgen.kb import load_gene_expression, load_genes
    panel = parca.calculate_parameters()
    expr = load_gene_expression()
    rows = []
    for g in load_genes():
        key = (g.get("symbol") or "").strip() or g["gene_id"]
        if key not in panel:
            continue
        e = expr.get(g["gene_id"], {})
        obs = e.get("expression_37C") or e.get("expression_mean")
        rows.append((key, e.get("rna_type") or "mRNA", obs, panel[key][0], panel[key][1] / 60.0))
    return rows, parca.metabolic_demand(), parca.close_metabolic_loop()


_COLORS = {"mRNA": "#1f77b4", "tRNA": "#2ca02c", "rRNA": "#d62728", "sRNA": "#9467bd"}


def main() -> int:
    from viva_mgen import viz
    VIZ.mkdir(parents=True, exist_ok=True)
    rows, demand, loop = _data()

    # A — fitted synthesis rate vs observed expression, per RNA type (log-log)
    figA = go.Figure(layout=viz._layout(
        "Fitted synthesis rate vs observed expression (ParCa)",
        "observed expression (Weiner 2003)", "fitted synthesis rate (1/s)"))
    for rt in ("mRNA", "tRNA", "rRNA", "sRNA"):
        pts = [(o, s) for _, t, o, s, _ in rows if t == rt and o and o > 0]
        if pts:
            figA.add_trace(go.Scattergl(x=[p[0] for p in pts], y=[p[1] for p in pts],
                           mode="markers", name=rt,
                           marker=dict(color=_COLORS.get(rt, "#888"), size=5, opacity=0.6)))
    figA.update_xaxes(type="log"); figA.update_yaxes(type="log")
    (VIZ / "parca_a_synthesis_vs_expression.html").write_text(viz._page(figA, "ParCa A"))

    # B — half-life by RNA type
    by_type = {}
    for _, t, _, _, hl in rows:
        by_type.setdefault(t, []).append(hl)
    cats = [t for t in ("mRNA", "tRNA", "rRNA", "sRNA") if t in by_type]
    means = [float(np.mean(by_type[t])) for t in cats]
    figB = go.Figure(layout=viz._layout("Per-gene RNA half-life by type (ParCa)",
                                        "RNA type", "half-life (min)"))
    figB.add_trace(go.Bar(x=cats, y=means, marker_color=[_COLORS[t] for t in cats],
                          text=[f"{m:.1f}" for m in means], textposition="outside"))
    figB.update_yaxes(type="log")
    (VIZ / "parca_b_halflife_by_type.html").write_text(viz._page(figB, "ParCa B"))

    # C — metabolic NMP demand (AT-rich)
    nmp = demand.get("nmp", {})
    figC = go.Figure(layout=viz._layout("Ribonucleotide demand the expression implies (ParCa)",
                                        "ribonucleotide", "fraction of RNA demand"))
    figC.add_trace(go.Bar(x=list(nmp.keys()), y=list(nmp.values()),
                          marker_color=["#d62728", "#1f77b4", "#2ca02c", "#ff7f0e"],
                          text=[f"{v:.2f}" for v in nmp.values()], textposition="outside"))
    (VIZ / "parca_c_metabolic_demand.html").write_text(viz._page(figC, "ParCa C"))

    # D — expression↔metabolism closure supply fluxes
    figD = go.Figure(layout=viz._layout("Expression↔metabolism closure — precursor supply (ParCa)",
                                        "", "max supply flux (feasible > 0)"))
    figD.add_trace(go.Bar(x=["ribonucleotides", "amino acids"],
                          y=[loop.get("nmp_supply_flux", 0.0), loop.get("aa_supply_flux", 0.0)],
                          marker_color=["#1f77b4", "#2ca02c"],
                          text=[f"{loop.get('nmp_supply_flux', 0):.2f}",
                                f"{loop.get('aa_supply_flux', 0):.2f}"], textposition="outside"))
    (VIZ / "parca_d_closure.html").write_text(viz._page(figD, "ParCa D"))

    _combined(rows, demand, loop, cats, means, nmp)
    _rewire()
    print(f"parca figures regenerated: {len(rows)} genes, NMP A+U="
          f"{nmp.get('A', 0) + nmp.get('U', 0):.2f}, closed-loop feasible={bool(loop.get('feasible'))}")
    return 0


def _combined(rows, demand, loop, cats, means, nmp):
    fig = make_subplots(rows=2, cols=2, subplot_titles=(
        "(a) fitted synthesis vs observed expression", "(b) RNA half-life by type",
        "(c) ribonucleotide demand (AT-rich)", "(d) expression↔metabolism closure"))
    for rt in ("mRNA", "tRNA", "rRNA", "sRNA"):
        pts = [(o, s) for _, t, o, s, _ in rows if t == rt and o and o > 0]
        if pts:
            fig.add_trace(go.Scattergl(x=[p[0] for p in pts], y=[p[1] for p in pts], mode="markers",
                          name=rt, marker=dict(color=_COLORS.get(rt, "#888"), size=4, opacity=0.6)), 1, 1)
    fig.update_xaxes(type="log", row=1, col=1); fig.update_yaxes(type="log", row=1, col=1)
    fig.add_trace(go.Bar(x=cats, y=means, marker_color=[_COLORS[t] for t in cats], showlegend=False), 1, 2)
    fig.update_yaxes(type="log", row=1, col=2)
    fig.add_trace(go.Bar(x=list(nmp.keys()), y=list(nmp.values()), showlegend=False,
                         marker_color=["#d62728", "#1f77b4", "#2ca02c", "#ff7f0e"]), 2, 1)
    fig.add_trace(go.Bar(x=["ribonucleotides", "amino acids"],
                         y=[loop.get("nmp_supply_flux", 0.0), loop.get("aa_supply_flux", 0.0)],
                         showlegend=False, marker_color=["#1f77b4", "#2ca02c"]), 2, 2)
    fig.update_layout(height=820, width=1050, template="plotly_white",
                      title_text="Parameter calculator (ParCa) — Karr 2012 FitConstants reproduction (viva-Mgen)")
    (VIZ / "parca_combined.html").write_text(fig.to_html(full_html=True, include_plotlyjs="cdn"))


def _rewire():
    import yaml
    p = WS / "workspace/studies/parca-parameter-fitting/study.yaml"
    d = yaml.safe_load(p.read_text())
    panels = [
        ("parca-combined", "parca_combined.html", "ParCa summary — panels a–d: fit, half-lives, demand, closure."),
        ("parca-a-synthesis-vs-expression", "parca_a_synthesis_vs_expression.html",
         "Fitted synthesis rate vs observed expression, per RNA type (log-log)."),
        ("parca-b-halflife-by-type", "parca_b_halflife_by_type.html",
         "Per-gene RNA half-life by type — mRNA short-lived, tRNA/rRNA stable."),
        ("parca-c-metabolic-demand", "parca_c_metabolic_demand.html",
         "Ribonucleotide demand the transcriptome implies — AT-rich (A+U dominant)."),
        ("parca-d-closure", "parca_d_closure.html",
         "Expression↔metabolism closure — iPS189 precursor supply fluxes (feasible > 0)."),
    ]
    d["visualizations"] = [{"name": n, "chart": f"viz/{f}", "render": "python scripts/regen_parca.py"}
                           for n, f, _ in panels]
    d["embed_visualizations"] = [{"name": n, "url": f"/workspace/studies/{SLUG}/viz/{f}", "description": desc}
                                 for n, f, desc in panels]
    p.write_text(yaml.dump(d, sort_keys=False, width=100, allow_unicode=True))


if __name__ == "__main__":
    main()
