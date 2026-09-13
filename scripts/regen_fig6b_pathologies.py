#!/usr/bin/env python3
"""Regenerate Figure 6B (common molecular pathologies of single-gene disruptions).

Matches Karr 2012 Fig 6B: a grid whose COLUMNS are phenotype classes (WT, then a
representative disruption of each essential category) and whose ROWS are the
single-cell temporal dynamics of Growth (fg·h⁻¹), Protein (fg), RNA (fg), DNA
(fg) and Septum (nm). Each disruption's trace is overlaid on the WT reference
(grey); cells whose dynamics differ significantly from WT are shaded red — the
paper's way of showing WHICH macromolecule/process each disruption breaks.

Each column runs the whole-cell composite with the relevant submodel impaired
(reduced but genuine: knocking out a process makes its product plateau, and the
downstream cascade — e.g. no mRNA ⇒ no protein — is emergent):
  WT · Metabolic (nutrient-limited) · RNA (Δtranscription) · Protein (Δtranslation)
  · Other (Δtranslocation) · DNA (Δreplication) · Cytokinesis (Δcytokinesis)
  · Quasi-Ess (mild nutrient limitation)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

WS = Path(__file__).resolve().parents[1]
SLUG = "fig6-gene-essentiality"
VIZ = WS / "workspace" / "studies" / SLUG / "viz"

T_H = 15.0
DT = 300.0
BIRTH_MASS = 3.93
F_PROT, F_RNA, F_DNA = 0.62, 0.093, 0.169  # dry-mass fractions (fg per fg mass)

# column label (paper gene in parens) -> composite modification
COLUMNS = [
    ("WT", {}),
    ("Metabolic\n(trmK)", {"nutrient_scale": 0.08}),
    ("RNA\n(rpoE)", {"drop": "transcription"}),
    ("Protein\n(asnS)", {"drop": "translation"}),
    ("Other\n(ffh)", {"drop": "protein_translocation"}),
    ("DNA\n(dnaN)", {"drop": "replication"}),
    ("Cytokinesis\n(parC)", {"drop": "cytokinesis"}),
    ("Quasi-Ess\n(tilS)", {"nutrient_scale": 0.55}),
]
ROWS = ["Growth (fg·h⁻¹)", "Protein (fg)", "RNA (fg)", "DNA (fg)", "Septum (nm)"]


def _simulate(core, mod):
    from viva_mgen.composites import build_mgen
    from process_bigraph import Composite, gather_emitter_results
    kw = {k: v for k, v in mod.items() if k != "drop"}
    doc = build_mgen(core, interval=DT, **kw)
    if mod.get("drop") and mod["drop"] in doc:
        del doc[mod["drop"]]
    sim = Composite({"state": doc}, core=core)
    sim.run(T_H * 3600.0)
    rows = [r for r in gather_emitter_results(sim)[("emitter",)] if r]
    def ser(k):
        return np.array([float(r.get(k, 0.0)) for r in rows])
    def msum(k):
        return np.array([sum((r.get(k, {}) or {}).values()) for r in rows])
    t = np.arange(len(rows)) * DT / 3600.0
    return dict(t=t, mass=ser("mass"), prot=msum("protein_counts"),
                rna=msum("rna_counts"), repl=ser("replicated_fraction"),
                sept=ser("septum_diameter"))


def _readouts(d, ref):
    """Five paper-row readouts in representative fg / (fg·h⁻¹) / nm.

    Protein & RNA are calibrated so a WT cell roughly doubles birth content (they
    track the model's cumulative synthesis, normalised to the WT trajectory);
    DNA and Septum are direct model outputs; Growth is d(mass)/dt.
    """
    growth = np.gradient(d["mass"], d["t"])                      # fg·h⁻¹
    p_ref = ref["prot"].max() or 1.0
    r_ref = ref["rna"].max() or 1.0
    protein = BIRTH_MASS * F_PROT * (1.0 + d["prot"] / p_ref)    # ~2× birth in WT
    rna = BIRTH_MASS * F_RNA * (1.0 + d["rna"] / r_ref)
    dna = BIRTH_MASS * F_DNA * (1.0 + d["repl"])                 # doubles on replication
    return [growth, protein, rna, dna, d["sept"]]


def main():
    from viva_mgen.core import build_core
    VIZ.mkdir(parents=True, exist_ok=True)
    core = build_core()

    sims = [(_label, _simulate(core, mod)) for _label, mod in COLUMNS]
    wt = sims[0][1]
    wt_ro = _readouts(wt, wt)
    all_ro = [(_label, _readouts(d, wt)) for _label, d in sims]

    ncol, nrow = len(COLUMNS), len(ROWS)
    fig = make_subplots(rows=nrow, cols=ncol, shared_xaxes=True,
                        column_titles=[c[0] for c in COLUMNS],
                        row_titles=ROWS, horizontal_spacing=0.012, vertical_spacing=0.03)
    for ci, (label, ro) in enumerate(all_ro):
        for ri in range(nrow):
            y = ro[ri]; wy = wt_ro[ri]
            # deviation vs WT (relative RMS over the trajectory)
            denom = (np.abs(wy).max() or 1.0)
            dev = np.sqrt(np.mean((y - wy) ** 2)) / denom
            if ci > 0 and dev > 0.15:
                fig.add_vrect(x0=0, x1=T_H, fillcolor="rgba(214,39,40,0.10)",
                              line_width=0, row=ri + 1, col=ci + 1)
            if ci > 0:  # WT reference behind each disruption
                fig.add_trace(go.Scatter(x=wt["t"], y=wy, mode="lines",
                                         line=dict(color="rgba(150,150,150,0.6)", width=1),
                                         showlegend=False), ri + 1, ci + 1)
            fig.add_trace(go.Scatter(x=sims[ci][1]["t"], y=y, mode="lines",
                                     line=dict(color="#111", width=1.4),
                                     showlegend=False), ri + 1, ci + 1)
    fig.update_xaxes(showticklabels=False)
    for ci in range(ncol):
        fig.update_xaxes(showticklabels=True, title_text="time (h)", row=nrow, col=ci + 1)
    fig.update_annotations(font_size=11)
    fig.update_layout(height=900, width=1500, template="plotly_white",
                      title_text="Figure 6B — Common molecular pathologies of single-gene disruptions (viva-Mgen)")
    (VIZ / "fig6b_pathologies.html").write_text(fig.to_html(full_html=True, include_plotlyjs="cdn"))

    _rewire()
    # brief per-column summary
    for label, ro in all_ro:
        print(f"{label.replace(chr(10),' '):22s} growth_end={ro[0][-1]:.2f} prot_end={ro[1][-1]:.2f} "
              f"rna_end={ro[2][-1]:.2f} dna_end={ro[3][-1]:.2f} sept_end={ro[4][-1]:.0f}")


def _rewire():
    import yaml
    p = WS / "workspace/studies/fig6-gene-essentiality/study.yaml"
    d = yaml.safe_load(p.read_text())
    viz = d.get("embed_visualizations") or []
    # keep existing 6A panels; add/replace the 6B grid
    viz = [v for v in viz if v.get("name") != "fig6-b-pathologies"]
    viz.insert(1, {"name": "fig6-b-pathologies",
                   "url": f"/workspace/studies/{SLUG}/viz/fig6b_pathologies.html",
                   "description": "Fig 6B: molecular pathologies — Growth/Protein/RNA/DNA/Septum dynamics for WT + a representative disruption of each essential class; red = differs from WT."})
    d["embed_visualizations"] = viz
    vs = d.get("visualizations") or []
    vs = [v for v in vs if v.get("name") != "fig6-b-pathologies"]
    vs.append({"name": "fig6-b-pathologies", "chart": "viz/fig6b_pathologies.html",
               "render": "python scripts/regen_fig6b_pathologies.py"})
    d["visualizations"] = vs
    p.write_text(yaml.dump(d, sort_keys=False, width=100, allow_unicode=True))


if __name__ == "__main__":
    main()
