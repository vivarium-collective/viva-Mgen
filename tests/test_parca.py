"""The native ParCa computes per-gene parameters from the real observed
expression (Weiner 2003, decoded from the KB), reproducing the initialization
stage of Karr 2012 FitConstants."""

import math

from viva_mgen import parca
from viva_mgen.kb import load_gene_expression, load_genes


def test_observed_expression_available():
    expr = load_gene_expression()
    assert len(expr) >= 500, "decoded observed gene expression should cover ~525 genes"
    # rRNA genes are the highest-expressed (validation of the decode + alignment)
    top = sorted(expr.items(), key=lambda kv: -(kv[1]["expression_mean"] or 0))[:3]
    assert all("rrn" in gid.lower() for gid, _ in top), [gid for gid, _ in top]


def test_parca_panel_from_real_expression():
    panel = parca.calculate_parameters()
    assert len(panel) >= 500
    # every entry is a 5-tuple of positive rates / lengths
    for key, v in panel.items():
        synth, half_s, transl, prot_half, length = v
        assert synth > 0 and half_s > 0 and transl > 0 and length > 0

    # half-lives are the real decoded per-gene KB values: rRNA/tRNA far more stable
    # than every mRNA. Classify by the REAL KB RNA type (not the name heuristic).
    expr = load_gene_expression()
    key_type = {}
    for g in load_genes():
        k = (g.get("symbol") or "").strip() or g["gene_id"]
        key_type[k] = (expr.get(g["gene_id"], {}).get("rna_type") or "").strip()
    mrna = [v[1] for k, v in panel.items() if key_type.get(k) == "mRNA"]
    rrna = [v[1] for k, v in panel.items() if key_type.get(k) == "rRNA"]
    assert rrna and max(mrna) < min(rrna), "rRNA half-life must exceed every mRNA half-life"

    # abundant EF-Tu (tuf) is synthesized well above the median (real expression signal)
    synth = sorted(v[0] for v in panel.values())
    median = synth[len(synth) // 2]
    if "tuf" in panel:
        assert panel["tuf"][0] > median


def test_decay_rates_match_halflives():
    rates = parca.mrna_decay_rates()
    panel = parca.calculate_parameters()
    for k in list(panel)[:20]:
        assert math.isclose(rates[k], math.log(2.0) / panel[k][1], rel_tol=1e-6)
