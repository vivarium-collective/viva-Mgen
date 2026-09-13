#!/usr/bin/env python3
"""Regenerate Figure 1 (whole-cell integration) for viva-Mgen.

Matches Karr 2012 Fig 1: the whole cell is 28+ independent submodels integrated
through shared cell variables. Panels: (a) sunburst of the submodels grouped by
process category, (b) the process ↔ cell-variable wiring matrix (integration
diagram), (c) an illustrative integrated single-cell time-course over the full
9 h cycle on one shared clock (the paper's Fig 1 is a schematic; this extra panel
shows the composed cell actually running). Hierarchy-aware: reads the new
``["cell", <compartment>, <name>]`` store paths.

Writes fig1_a..c + fig1_combined and rewires the study's visualizations via a
YAML round-trip.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

WS = Path(__file__).resolve().parents[1]
SLUG = "fig1-architecture"
VIZ = WS / "workspace" / "studies" / SLUG / "viz"

# submodel → process category (Karr Fig 1 colour groups)
CATEGORY = {
    "metabolism": "Metabolism",
    "transcription": "RNA", "rna_decay": "RNA", "transcriptional_regulation": "RNA",
    "rna_processing": "RNA", "rna_modification": "RNA", "trna_aminoacylation": "RNA",
    "translation": "Protein", "protein_decay": "Protein", "protein_processing_i": "Protein",
    "protein_translocation": "Protein", "protein_processing_ii": "Protein",
    "protein_folding": "Protein", "protein_modification": "Protein",
    "protein_activation": "Protein", "complexation": "Protein",
    "ribosome_assembly": "Protein", "terminal_organelle_assembly": "Protein",
    "replication": "DNA-chromosome", "replication_initiation": "DNA-chromosome",
    "supercoiling": "DNA-chromosome", "condensation": "DNA-chromosome",
    "segregation": "DNA-chromosome", "dna_damage": "DNA-chromosome",
    "dna_repair": "DNA-chromosome", "chromosome": "DNA-chromosome",
    "ftsz": "Cytokinesis", "cytokinesis": "Cytokinesis",
    "mass": "Other", "host_interaction": "Other",
}


def _store_keys(port_map):
    keys = set()
    for path in (port_map or {}).values():
        # new hierarchy: ["cell", <compartment>, <name>]; also accept legacy ["stores", name]
        if isinstance(path, (list, tuple)) and len(path) >= 2 and path[0] in ("cell", "stores"):
            keys.add(path[-1])
    return keys


def _wiring(doc):
    procs = [k for k, v in doc.items() if isinstance(v, dict) and v.get("_type") == "process"]
    touched = {p: _store_keys(doc[p].get("inputs")) | _store_keys(doc[p].get("outputs")) for p in procs}
    stores = []
    for p in procs:
        for s in sorted(touched[p]):
            if s not in stores:
                stores.append(s)
    mat = [[1 if s in touched[p] else 0 for s in stores] for p in procs]
    return procs, stores, mat


def _run(runtime_s=32400.0, dt=300.0):
    from viva_mgen.core import build_core
    from viva_mgen.composites import build_mgen
    from process_bigraph import Composite, gather_emitter_results
    core = build_core()
    doc = build_mgen(core, interval=dt)
    procs, stores, mat = _wiring(doc)
    sim = Composite({"state": doc}, core=core)
    sim.run(runtime_s)
    rows = [r for r in gather_emitter_results(sim)[("emitter",)] if r]
    return procs, stores, mat, rows, runtime_s, dt


def main():
    from viva_mgen import viz
    PAL = viz.PALETTE
    VIZ.mkdir(parents=True, exist_ok=True)
    procs, stores, mat, rows, runtime_s, dt = _run()

    # (a) sunburst of submodels by category
    labels, parents = ["cell"], [""]
    cats = sorted(set(CATEGORY.get(p, "Other") for p in procs))
    for c in cats:
        labels.append(c); parents.append("cell")
    for p in procs:
        labels.append(p); parents.append(CATEGORY.get(p, "Other"))
    figa = go.Figure(go.Sunburst(labels=labels, parents=parents, branchvalues="remainder"))
    figa.update_layout(title="28-submodel architecture, grouped by process category (Fig 1A)",
                       template="plotly_white", margin=dict(t=60, l=0, r=0, b=0))
    (VIZ / "fig1_a_architecture.html").write_text(figa.to_html(full_html=True, include_plotlyjs="cdn"))

    # (b) wiring matrix
    (VIZ / "fig1_b_wiring.html").write_text(viz.heatmap_html(
        "Process ↔ cell-variable wiring (Fig 1B)", z=mat, x=stores, y=procs,
        x_title="cell variable (store)", y_title="submodel (process)", colorbar_title="touches"))

    # (c) integrated dynamics over the full 9 h cycle
    n = len(rows)
    t_h = np.linspace(0.0, runtime_s / 3600.0, n) if n > 1 else np.array([0.0])
    def series(key):
        return np.array([float(r.get(key, 0.0)) for r in rows])
    def mapsum(key):
        return np.array([sum((r.get(key, {}) or {}).values()) for r in rows])
    mass = series("mass"); repl = series("replicated_fraction")
    gr = series("growth_rate")
    tot_mrna = mapsum("rna_counts"); tot_prot = mapsum("protein_counts")
    def norm(a):
        m = a.max() or 1.0
        return a / m
    (VIZ / "fig1_c_dynamics.html").write_text(viz.line_series_html(
        "Integrated cell dynamics on one shared clock (Fig 1)", t_h, {
            "mass (norm)": norm(mass),
            "growth rate (norm)": norm(gr),
            "replicated fraction": repl,
            "total mRNA (norm)": norm(tot_mrna),
            "total protein (norm)": norm(tot_prot),
        }, x_title="time (h)", y_title="normalized level"))

    # combined
    fig = make_subplots(rows=2, cols=2, specs=[[{"type": "domain"}, {"type": "heatmap"}],
                                               [{"type": "xy", "colspan": 2}, None]],
                        subplot_titles=("(a) 28-submodel architecture", "(b) process ↔ cell-variable wiring",
                                        "(c) integrated cell dynamics over the 9 h cycle"),
                        vertical_spacing=0.16, row_heights=[0.56, 0.44])
    fig.add_trace(go.Sunburst(labels=labels, parents=parents, branchvalues="remainder"), 1, 1)
    fig.add_trace(go.Heatmap(z=mat, x=stores, y=procs, colorscale="Blues", showscale=False), 1, 2)
    for i, (name, ys) in enumerate([("mass", norm(mass)), ("growth rate", norm(gr)),
                                    ("replicated fraction", repl), ("total mRNA", norm(tot_mrna)),
                                    ("total protein", norm(tot_prot))]):
        fig.add_trace(go.Scatter(x=t_h, y=ys, mode="lines", name=name,
                                 line=dict(color=PAL[i % len(PAL)])), 2, 1)
    # The 87 store labels are unreadable in this cell (and their rotated text
    # collides with panel c) — hide them here; panel b standalone keeps them.
    fig.update_xaxes(showticklabels=False, title_text="87 cell variables (see panel b)", row=1, col=2)
    fig.update_yaxes(tickfont=dict(size=7), row=1, col=2)
    fig.update_xaxes(title_text="time (h)", row=2, col=1)
    fig.update_yaxes(title_text="normalized", row=2, col=1)
    fig.update_layout(height=1150, width=1050, template="plotly_white", legend=dict(y=0.42),
                      title_text="Figure 1 — The whole-cell model integrates 28 submodels through shared cell variables (viva-Mgen)")
    (VIZ / "fig1_combined.html").write_text(fig.to_html(full_html=True, include_plotlyjs="cdn"))

    _rewire()
    print(f"fig1 regenerated: {len(procs)} processes, {len(stores)} stores, {n} time points over "
          f"{runtime_s/3600:.0f} h; replicated_fraction {repl[0]:.2f}->{repl[-1]:.2f}")


def _rewire():
    import yaml
    p = WS / "workspace/studies/fig1-architecture/study.yaml"
    d = yaml.safe_load(p.read_text())
    panels = [
        ("fig1-combined", "fig1_combined.html", "Figure 1 reproduced — whole-cell integration: 28-submodel architecture, wiring matrix, and integrated 9 h dynamics."),
        ("fig1-a-architecture", "fig1_a_architecture.html", "Fig 1A: sunburst of the 28 submodels grouped by process category."),
        ("fig1-b-wiring", "fig1_b_wiring.html", "Fig 1B: process ↔ cell-variable wiring (integration diagram)."),
        ("fig1-c-dynamics", "fig1_c_dynamics.html", "Fig 1C: integrated single-cell dynamics over the full 9 h cell cycle."),
    ]
    d["visualizations"] = [{"name": n, "chart": f"viz/{f}", "render": "python scripts/regen_fig1_architecture.py"} for n, f, _ in panels]
    d["embed_visualizations"] = [{"name": n, "url": f"/workspace/studies/{SLUG}/viz/{f}", "description": desc}
                                 for n, f, desc in panels]
    p.write_text(yaml.dump(d, sort_keys=False, width=100, allow_unicode=True))


if __name__ == "__main__":
    main()
