#!/usr/bin/env python
"""Regenerate the Visualizations for the structural-model studies (s01/s02).

Builds a cohesive, dark-themed set of interactive figures from each study's
packed-cell metadata + the shared roster, covering the molecules, their counts,
the structure-preparation process, and a full 3D interactive cell (three.js,
embedded from the workbench's parsimony viewer).

    PYTHONPATH=<worktree> python scripts/regen_structural_viz.py

Reads:  workspace/studies/<slug>/viz/3d/*.meta.json  (per-species placed counts,
        category, structure source) + viva_mgen.structural roster/counts.
Writes: workspace/studies/<slug>/viz/<name>.html     (self-contained Plotly /
        iframe HTML), registered in each study.yaml's visualizations: block.
"""
from __future__ import annotations

import csv
import glob
import json
from collections import Counter, defaultdict
from pathlib import Path

import plotly.graph_objects as go

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "viva_mgen" / "structural" / "data"

# ── theme ────────────────────────────────────────────────────────────────────
BG = "#0d1017"
PANEL = "#141925"
INK = "#e7ecf3"
DIM = "#8b98ad"
GRID = "rgba(139,152,173,0.14)"
ACCENT = "#5eb0ff"
FONT = "Inter, -apple-system, Segoe UI, Helvetica, Arial, sans-serif"

# functional-category palette (coherent with the 3D viewer's grouping)
CAT_COLORS = {
    "metabolism": "#5ed0c5", "translation": "#f6c453", "transcription": "#e88ad0",
    "RNA synthesis/maturation": "#c792ea", "protein transport/singaling": "#7ea6ff",
    "protein folding/maturation": "#8bd450", "cytokinesis/motility": "#ff9e64",
    "DNA replication/maintenance": "#ff6b8b", "MG-specific": "#b0bec5",
    "host cell interaction": "#f78c6c", "lipoprotein": "#89ddff",
    "uncharacterized": "#5b6472", "Nucleoid": "#d8c69a", "RNA": "#e0a44b",
}
DEFAULT_C = "#6b7688"


def cat_color(c):
    return CAT_COLORS.get(c, DEFAULT_C)


def layout(fig, title, subtitle="", height=460):
    fig.update_layout(
        title=dict(text=f"<b>{title}</b>" + (f"<br><span style='font-size:12px;color:{DIM}'>{subtitle}</span>" if subtitle else ""),
                   x=0.02, xanchor="left", font=dict(size=17, color=INK)),
        paper_bgcolor=BG, plot_bgcolor=BG, font=dict(family=FONT, color=INK, size=12),
        margin=dict(l=16, r=16, t=64 if subtitle else 52, b=16), height=height,
        colorway=list(CAT_COLORS.values()),
        legend=dict(bgcolor="rgba(0,0,0,0)", font=dict(color=DIM, size=11)),
    )
    fig.update_xaxes(gridcolor=GRID, zerolinecolor=GRID, linecolor=GRID, tickfont=dict(color=DIM))
    fig.update_yaxes(gridcolor=GRID, zerolinecolor=GRID, linecolor=GRID, tickfont=dict(color=DIM))
    return fig


def write(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(path), include_plotlyjs="cdn", full_html=True,
                   config={"displayModeBar": False, "responsive": True})
    # dark page background around the plot
    html = path.read_text().replace("<body>", f"<body style='margin:0;background:{BG}'>", 1)
    path.write_text(html)


# ── data ─────────────────────────────────────────────────────────────────────
def roster():
    return {r["prot_id"]: r for r in csv.DictReader(open(DATA / "mgen_proteins.csv"))}


def study_meta(slug):
    f = glob.glob(str(ROOT / "workspace/studies" / slug / "viz/3d/*.meta.json"))
    return json.load(open(f[0]))["ingredients"] if f else {}


def placed(meta):
    """(prot_id, name, category, count, structure_db) for count>0."""
    out = []
    for pid, v in meta.items():
        c = int(v.get("count", 0))
        if c > 0:
            st = (v.get("structure") or {})
            out.append((pid, v.get("display_name") or pid, v.get("category") or "uncharacterized",
                        c, st.get("db")))
    return out


# ── 3D interactive cell (three.js, embedded from the workbench viewer) ─────────
def write_3d_embed(slug, pack_rel, out: Path, label):
    pack_url = f"/workspace/studies/{slug}/viz/3d/{pack_rel}"
    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>{label} — 3D cell</title><style>
  html,body{{margin:0;height:100%;background:{BG};font-family:{FONT};color:{INK}}}
  .bar{{padding:8px 14px;font-size:13px;color:{DIM};border-bottom:1px solid {GRID}}}
  .bar b{{color:{INK}}} iframe{{border:0;width:100%;height:calc(100% - 37px);display:block}}
</style></head><body>
  <div class="bar"><b>Interactive 3D {label} cell</b> — every molecule placed at true abundance
  (drag to orbit · scroll to zoom · toggle species in the panel · <b>About</b> for details)</div>
  <iframe src="/parsimony-viewer/?file={pack_url}" allow="xr-spatial-tracking"></iframe>
</body></html>"""
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html)


# ── charts ───────────────────────────────────────────────────────────────────
def fig_inventory_sunburst(pl, label):
    by_cat = defaultdict(list)
    for pid, name, cat, c, _ in pl:
        by_cat[cat].append((name, c))
    ids, labels, parents, values, colors = ["cell"], ["Whole cell"], [""], [sum(c for *_, c, _ in pl)], [PANEL]
    for cat, items in sorted(by_cat.items(), key=lambda kv: -sum(c for _, c in kv[1])):
        cid = f"cat::{cat}"
        ids.append(cid); labels.append(cat); parents.append("cell")
        values.append(sum(c for _, c in items)); colors.append(cat_color(cat))
        for name, c in sorted(items, key=lambda x: -x[1])[:14]:
            ids.append(f"{cid}::{name}"); labels.append(name); parents.append(cid)
            values.append(c); colors.append(cat_color(cat))
    fig = go.Figure(go.Sunburst(ids=ids, labels=labels, parents=parents, values=values,
                                branchvalues="total", marker=dict(colors=colors, line=dict(color=BG, width=1)),
                                insidetextorientation="radial",
                                hovertemplate="<b>%{label}</b><br>%{value:,} molecules<extra></extra>"))
    return layout(fig, f"Molecular inventory — {label}",
                  "molecules by functional category (top species shown); click a wedge to zoom", 560)


def fig_copy_numbers(pl, label):
    top = sorted(pl, key=lambda x: -x[3])[:30][::-1]
    fig = go.Figure(go.Bar(
        x=[c for *_, c, _ in top], y=[n for _, n, *_ in top], orientation="h",
        marker=dict(color=[cat_color(cat) for _, _, cat, _, _ in top]),
        hovertemplate="<b>%{y}</b><br>%{x:,} copies<extra></extra>"))
    fig.update_xaxes(type="log", title="copies per cell (log)")
    return layout(fig, f"Copy numbers — {label}",
                  "the 30 most abundant molecular species, coloured by function", 640)


def fig_spatial_classes(pl, ros, label):
    cls = Counter()
    for pid, name, cat, c, _ in pl:
        r = ros.get(pid, {})
        if cat == "RNA":
            cls["tRNA (free RNA)"] += c
        elif r.get("dna_binding"):
            cls["DNA-bound (nucleoid)"] += c
        elif r.get("compartment") == "m":
            cls["membrane"] += c
        elif r.get("compartment") == "e":
            cls["extracellular"] += c
        else:
            cls["cytoplasm"] += c
    order = ["cytoplasm", "membrane", "extracellular", "DNA-bound (nucleoid)", "tRNA (free RNA)"]
    cols = {"cytoplasm": "#5ed0c5", "membrane": "#7ea6ff", "extracellular": "#f78c6c",
            "DNA-bound (nucleoid)": "#ff6b8b", "tRNA (free RNA)": "#e0a44b"}
    items = [(k, cls.get(k, 0)) for k in order if cls.get(k, 0)]
    fig = go.Figure(go.Bar(x=[v for _, v in items], y=[k for k, _ in items], orientation="h",
                           marker=dict(color=[cols[k] for k, _ in items]),
                           text=[f"{v:,}" for _, v in items], textposition="auto",
                           hovertemplate="<b>%{y}</b><br>%{x:,} molecules<extra></extra>"))
    return layout(fig, f"Spatial distribution — {label}",
                  "molecules by Maritan spatial class (cytoplasm · membrane · extracellular · DNA-bound · free RNA)", 400)


def fig_structure_provenance(ros, label):
    # roster-wide preparation outcome (how each species got a structure)
    prov = Counter()
    for pid, r in ros.items():
        if r.get("pdb_id"):
            prov["curated experimental (PDB/mmCIF)"] += 1
        elif r.get("uniprot"):
            prov["AlphaFold (by UniProt)"] += 1
        else:
            prov["no structure (complex/other)"] += 1
    cols = {"curated experimental (PDB/mmCIF)": "#8bd450", "AlphaFold (by UniProt)": "#7ea6ff",
            "no structure (complex/other)": "#5b6472"}
    labels = list(prov); values = [prov[k] for k in labels]
    fig = go.Figure(go.Pie(labels=labels, values=values, hole=0.58,
                           marker=dict(colors=[cols[k] for k in labels], line=dict(color=BG, width=2)),
                           textinfo="label+percent", textfont=dict(color=INK, size=11),
                           hovertemplate="<b>%{label}</b><br>%{value} species (%{percent})<extra></extra>"))
    fig.update_layout(annotations=[dict(text=f"<b>{sum(values)}</b><br>species", x=0.5, y=0.5,
                                        font=dict(color=INK, size=16), showarrow=False)])
    return layout(fig, "Structure preparation",
                  "how each of the 683 roster species was modelled (Maritan's curation workflow)", 460)


def fig_machines(pl, label):
    keys = {"RIBOSOME_70S": "70S ribosome", "RNA_POLYMERASE_HOLOENZYME": "RNA polymerase (holo)",
            "RNA_POLYMERASE": "RNA polymerase (core)", "DNA_GYRASE": "DNA gyrase",
            "GROEL_GROES": "GroEL/ES chaperonin", "DNA_POLYMERASE_HOLOENZYME": "DNA polymerase (holo)",
            "PYRUVATE_DEHYDROGENASE": "pyruvate dehydrogenase", "tRNA": "tRNA (pool)"}
    cnt = {pid: c for pid, _, _, c, _ in pl}
    items = [(lbl, cnt.get(pid, 0)) for pid, lbl in keys.items() if cnt.get(pid, 0)]
    items = sorted(items, key=lambda x: x[1])
    fig = go.Figure(go.Bar(x=[v for _, v in items], y=[k for k, _ in items], orientation="h",
                           marker=dict(color="#f6c453"), text=[f"{v:,}" for _, v in items],
                           textposition="auto", hovertemplate="<b>%{y}</b><br>%{x:,} copies<extra></extra>"))
    fig.update_xaxes(title="copies per cell")
    return layout(fig, f"Molecular machines — {label}",
                  "the large assembled complexes + the tRNA pool", 420)


def fig_pipeline(pl, ros, label):
    n_roster = len(ros)
    n_pdb = sum(1 for r in ros.values() if r.get("pdb_id"))
    n_af = sum(1 for r in ros.values() if r.get("uniprot") and not r.get("pdb_id"))
    n_none = n_roster - n_pdb - n_af
    n_placed_species = len(pl)
    n_mol = sum(c for *_, c, _ in pl)
    nodes = ["S1 roster (683 proteins)", "curated PDB", "AlphaFold", "no structure",
             "placed species", "molecules in cell"]
    ncol = ["#b0bec5", "#8bd450", "#7ea6ff", "#5b6472", "#5eb0ff", "#f6c453"]
    src = [0, 0, 0, 1, 2, 4]
    tgt = [1, 2, 3, 4, 4, 5]
    val = [n_pdb, n_af, n_none, n_pdb, n_af, n_mol]
    fig = go.Figure(go.Sankey(
        node=dict(label=nodes, color=ncol, pad=18, thickness=16, line=dict(color=BG, width=0)),
        link=dict(source=src, target=tgt, value=val,
                  color=["rgba(139,214,80,.35)", "rgba(126,166,255,.35)", "rgba(91,100,114,.3)",
                         "rgba(139,214,80,.35)", "rgba(126,166,255,.35)", "rgba(246,196,83,.35)"])))
    return layout(fig, f"Preparation pipeline — {label}",
                  f"roster → structure resolution → placement ({n_placed_species} species, {n_mol:,} molecules)", 440)


def fig_maritan_compare(pl_s01):
    cnt = {pid: c for pid, _, _, c, _ in pl_s01}
    mono = sum(c for pid, _, _, c, _ in pl_s01 if pid.endswith("_MONOMER"))
    ours = {"protein monomers": mono, "70S ribosomes": cnt.get("RIBOSOME_70S", 0),
            "RNA polymerase": cnt.get("RNA_POLYMERASE", 0) + cnt.get("RNA_POLYMERASE_HOLOENZYME", 0),
            "tRNAs": cnt.get("tRNA", 0)}
    maritan = {"protein monomers": 21377, "70S ribosomes": 69, "RNA polymerase": 77, "tRNAs": 1653}
    cats = list(ours)
    fig = go.Figure()
    fig.add_bar(name="Maritan 2022 (Table 1)", x=cats, y=[maritan[c] for c in cats], marker_color="#5b6472")
    fig.add_bar(name="viva-Mgen (s01)", x=cats, y=[ours[c] for c in cats], marker_color=ACCENT)
    fig.update_yaxes(type="log", title="copies per cell (log)")
    fig.update_layout(barmode="group")
    return layout(fig, "Reproduction vs Maritan 2022",
                  "our packed-cell counts against the published WC-MG figures (Table 1, birth cell)", 440)


def fig_baseline_vs_sim(pl_s01, pl_s02):
    def bycat(pl):
        d = Counter()
        for _, _, cat, c, _ in pl:
            if cat != "Nucleoid":
                d[cat] += c
        return d
    a, b = bycat(pl_s01), bycat(pl_s02)
    cats = sorted(set(a) | set(b), key=lambda c: -(a.get(c, 0)))
    fig = go.Figure()
    fig.add_bar(name="s01 — Maritan/WC-MG", x=cats, y=[a.get(c, 0) for c in cats], marker_color=ACCENT)
    fig.add_bar(name="s02 — viva_mgen sim", x=cats, y=[b.get(c, 0) for c in cats], marker_color="#f6c453")
    fig.update_layout(barmode="group")
    fig.update_yaxes(title="molecules")
    fig.update_xaxes(tickangle=-30)
    return layout(fig, "Published vs simulated proteome",
                  "molecules by function: the WC-MG baseline (s01) vs what viva_mgen predicts (s02)", 460)


# ── orchestration ────────────────────────────────────────────────────────────
def _pack_rel(slug):
    f = glob.glob(str(ROOT / "workspace/studies" / slug / "viz/3d/*.pack.json"))
    return Path(f[0]).name if f else None


def register(slug, entries):
    """Rewrite the study.yaml visualizations: block (name + chart path)."""
    import yaml
    sy = ROOT / "workspace/studies" / slug / "study.yaml"
    doc = yaml.safe_load(sy.read_text())
    doc["visualizations"] = [{"name": n, "chart": f"viz/{f}",
                              "render": "python scripts/regen_structural_viz.py"} for n, f in entries]
    sy.write_text(yaml.dump(doc, default_flow_style=False, sort_keys=False, allow_unicode=True, width=100))


def build_study(slug, label, ros, pl, cross):
    vdir = ROOT / "workspace/studies" / slug / "viz"
    entries = []
    pr = _pack_rel(slug)
    if pr:
        write_3d_embed(slug, pr, vdir / "3d_cell.html", label)
        entries.append(("3d-cell", "3d_cell.html"))
    charts = [
        ("inventory-sunburst", "inventory_sunburst.html", fig_inventory_sunburst(pl, label)),
        ("copy-numbers", "copy_numbers.html", fig_copy_numbers(pl, label)),
        ("spatial-classes", "spatial_classes.html", fig_spatial_classes(pl, ros, label)),
        ("molecular-machines", "molecular_machines.html", fig_machines(pl, label)),
        ("structure-prep", "structure_prep.html", fig_structure_provenance(ros, label)),
        ("preparation-pipeline", "preparation_pipeline.html", fig_pipeline(pl, ros, label)),
    ]
    for name, fn, fig in charts:
        write(fig, vdir / fn)
        entries.append((name, fn))
    for name, fn, fig in cross:
        write(fig, vdir / fn)
        entries.append((name, fn))
    register(slug, entries)
    print(f"{slug}: {len(entries)} visualizations")


def main():
    ros = roster()
    s01 = placed(study_meta("s01-maritan-baseline"))
    s02 = placed(study_meta("s02-reproduction-driven"))
    cmp_fig = ("maritan-comparison", "maritan_comparison.html", fig_maritan_compare(s01))
    dual_fig = ("baseline-vs-sim", "baseline_vs_sim.html", fig_baseline_vs_sim(s01, s02))
    build_study("s01-maritan-baseline", "Maritan baseline", ros, s01, [cmp_fig, dual_fig])
    build_study("s02-reproduction-driven", "reproduction-driven", ros, s02, [dual_fig])


if __name__ == "__main__":
    main()
