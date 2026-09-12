"""Shared interactive-visualization helpers for viva-Mgen studies.

Each helper returns a self-contained, theme-aware HTML page (pinned Plotly 3.x
CDN) that renders on a transparent surface so the host workbench card's theme
shows through; a small post-script re-themes chart chrome (font/axis/grid) on
light/dark toggle. Colors follow the dataviz skill's validated default palette.
"""

from __future__ import annotations

import json
import uuid
from typing import Sequence

import plotly.graph_objects as go

_PLOTLY_CDN = "https://cdn.plot.ly/plotly-3.7.0.min.js"

# validated categorical palette (light hexes; the post-script keeps chrome, not
# series, theme-aware — series hues read acceptably on both surfaces)
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
STATUS = {"good": "#0ca30c", "warning": "#fab219", "serious": "#ec835a", "critical": "#d03b3b"}

_POST_SCRIPT = """
(function(){
  var gd = document.getElementById('%(id)s');
  function chrome(){
    var root = document.documentElement;
    var t = root.getAttribute('data-theme');
    var dark = t === 'dark' || (t !== 'light' && window.matchMedia &&
               window.matchMedia('(prefers-color-scheme: dark)').matches);
    var font = dark ? '#ffffff' : '#0b0b0b';
    var sub  = dark ? '#c3c2b7' : '#52514e';
    var grid = dark ? '#333330' : '#e6e5e0';
    Plotly.relayout(gd, {
      'font.color': font, 'legend.font.color': sub,
      'xaxis.color': sub, 'yaxis.color': sub,
      'xaxis.gridcolor': grid, 'yaxis.gridcolor': grid,
      'xaxis.title.font.color': font, 'yaxis.title.font.color': font,
      'title.font.color': font
    });
  }
  chrome();
  if (window.matchMedia) window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', chrome);
  new MutationObserver(chrome).observe(document.documentElement, {attributes:true, attributeFilter:['data-theme']});
})();
"""


def _layout(title, x_title="", y_title=""):
    return go.Layout(
        title=dict(text=title, x=0.02, xanchor="left", font=dict(size=18)),
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="system-ui, -apple-system, Segoe UI, Roboto, sans-serif", size=13),
        xaxis=dict(title=x_title, zeroline=False, showline=True, ticks="outside"),
        yaxis=dict(title=y_title, zeroline=False, showline=True, ticks="outside"),
        margin=dict(l=64, r=28, t=54, b=52), hovermode="closest",
        legend=dict(orientation="h", y=-0.18, x=0),
    )


def _page(fig: go.Figure, title: str) -> str:
    div_id = "viz-" + uuid.uuid4().hex[:10]
    fig_json = fig.to_json()
    return f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<script src="{_PLOTLY_CDN}"></script>
<style>html,body{{margin:0;background:transparent}}#{div_id}{{width:100%;height:100vh;min-height:420px}}</style>
</head><body>
<div id="{div_id}"></div>
<script>
var spec = {fig_json};
Plotly.newPlot('{div_id}', spec.data, spec.layout, {{responsive:true, displayModeBar:false}});
{_POST_SCRIPT % {"id": div_id}}
</script></body></html>"""


def line_series_html(title, x, series: dict, x_title="", y_title="", modes=None) -> str:
    """Multi-series line chart. ``series`` maps name -> y-values (same length as x)."""
    fig = go.Figure(layout=_layout(title, x_title, y_title))
    for i, (name, y) in enumerate(series.items()):
        fig.add_trace(go.Scatter(
            x=list(x), y=list(y), mode=(modes or "lines"), name=name,
            line=dict(color=PALETTE[i % len(PALETTE)], width=2),
            marker=dict(color=PALETTE[i % len(PALETTE)], size=6),
        ))
    return _page(fig, title)


def scatter_html(title, x, y, x_title="", y_title="", color=None, text=None) -> str:
    fig = go.Figure(layout=_layout(title, x_title, y_title))
    fig.add_trace(go.Scatter(
        x=list(x), y=list(y), mode="markers", name=title,
        text=text, hovertemplate="(%{x:.3g}, %{y:.3g})<extra></extra>",
        marker=dict(color=color or PALETTE[0], size=8, opacity=0.75,
                    line=dict(width=1, color="rgba(255,255,255,0.6)")),
    ))
    return _page(fig, title)


def grouped_bar_html(title, categories: Sequence, groups: dict, x_title="", y_title="") -> str:
    """Grouped bars. ``groups`` maps series-name -> values aligned with categories."""
    fig = go.Figure(layout=_layout(title, x_title, y_title))
    for i, (name, vals) in enumerate(groups.items()):
        fig.add_trace(go.Bar(
            x=list(categories), y=list(vals), name=name,
            marker=dict(color=PALETTE[i % len(PALETTE)],
                        line=dict(width=2, color="rgba(255,255,255,0)")),
        ))
    fig.update_layout(barmode="group", bargap=0.25, bargroupgap=0.08)
    return _page(fig, title)


def donut_html(title, labels: Sequence, values: Sequence) -> str:
    fig = go.Figure(layout=_layout(title))
    fig.add_trace(go.Pie(
        labels=list(labels), values=list(values), hole=0.55, sort=False,
        marker=dict(colors=[PALETTE[i % len(PALETTE)] for i in range(len(labels))],
                    line=dict(color="rgba(255,255,255,0.9)", width=2)),
        textinfo="label+percent",
    ))
    fig.update_layout(showlegend=True)
    return _page(fig, title)


def heatmap_html(title, z, x=None, y=None, x_title="", y_title="", colorbar_title="") -> str:
    fig = go.Figure(layout=_layout(title, x_title, y_title))
    fig.add_trace(go.Heatmap(
        z=z, x=x, y=y, colorscale="Blues", reversescale=False,
        colorbar=dict(title=colorbar_title),
    ))
    return _page(fig, title)


def confusion_matrix_html(title, matrix, labels=("essential", "non-essential")) -> str:
    """2x2 confusion matrix (rows=reference, cols=model) as an annotated heatmap."""
    z = matrix
    xlab = [f"model {l}" for l in labels]
    ylab = [f"ref {l}" for l in labels]
    fig = go.Figure(layout=_layout(title, "model prediction", "reference"))
    fig.add_trace(go.Heatmap(z=z, x=xlab, y=ylab, colorscale="Blues",
                             showscale=True, colorbar=dict(title="count")))
    for i, row in enumerate(z):
        for j, v in enumerate(row):
            fig.add_annotation(x=xlab[j], y=ylab[i], text=str(int(v)),
                               showarrow=False, font=dict(size=20, color="#0b0b0b"))
    fig.update_yaxes(autorange="reversed")
    return _page(fig, title)
