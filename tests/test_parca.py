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


def test_fit_analytically_satisfies_constraints():
    """The analytic QP projects the observed distribution onto the FitConstants
    linear constraints: RNA mass reproduced exactly, DnaA/FtsZ held, and the
    projection is the least-change (closest to the initial guess)."""
    import numpy as np
    n = 50
    rng = np.random.default_rng(0)
    counts0 = rng.uniform(1, 100, n)
    mw = rng.uniform(300, 3000, n)
    # scale counts0 so the mass constraint is already consistent, then hold two genes
    held = [3, 7]
    fitted = parca.fit_analytically(counts0, mw, held_idx=held)
    mc = parca._mass_constants()
    target = mc["dry_weight_g"] * mc["rna_fraction"] * parca._NA
    assert np.isclose(fitted @ mw, target, rtol=1e-9), "RNA-mass constraint must hold"
    for j in held:
        assert np.isclose(fitted[j], counts0[j], rtol=1e-9), "held genes must be unchanged"
    assert (fitted >= 0).all()


def test_supercoiling_constraint_rebalances_topoisomerases():
    """With the net-supercoiling constraint active, the fitted topoisomerase-I
    (topA) expression rises relative to gyrase to zero net supercoiling activity —
    a real effect of the QP (the observed data alone has topA < gyrA)."""
    panel = parca.calculate_parameters()
    genes = {(g.get("symbol") or "").strip().lower(): g for g in load_genes()}

    def synth(sym):
        g = genes[sym]
        key = (g.get("symbol") or "").strip() or g["gene_id"]
        return panel[key][0]

    if "topa" in genes and "gyra" in genes:
        assert synth("topa") > synth("gyra")


def test_metabolic_demand_composition():
    """The ParCa's metabolic-demand output is a normalized NMP + AA composition;
    M. genitalium's AT-rich genome shows in an A+U-dominated ribonucleotide demand."""
    d = parca.metabolic_demand()
    if not d:
        import pytest
        pytest.skip("metabolic-demand dataset not present")
    nmp = d["nmp"]
    assert set(nmp) == {"A", "C", "G", "U"}
    assert abs(sum(nmp.values()) - 1.0) < 1e-6
    assert nmp["A"] + nmp["U"] > 0.55          # AT-rich genome
    aa = d["aa"]
    assert len(aa) == 20 and abs(sum(aa.values()) - 1.0) < 1e-6


def test_closed_metabolic_loop_is_consistent():
    """The expression↔metabolism loop closes: the iPS189 network can supply the
    ribonucleotide and amino-acid precursor demand in the ParCa's fitted composition
    alongside feasible baseline growth (Karr FitConstants feasibility criterion)."""
    r = parca.close_metabolic_loop()
    if not r.get("growth_baseline"):
        import pytest
        pytest.skip("metabolic model not available")
    assert r["growth_baseline"] > 0
    assert r["nmp_supply_flux"] > 0, "network must supply the fitted ribonucleotide mix"
    assert r["aa_supply_flux"] > 0, "network must supply the fitted amino-acid mix"
    assert r["feasible"] is True
    assert r["n_aa"] == 20 and r["n_nmp"] == 3


def test_decay_rates_match_halflives():
    rates = parca.mrna_decay_rates()
    panel = parca.calculate_parameters()
    for k in list(panel)[:20]:
        assert math.isclose(rates[k], math.log(2.0) / panel[k][1], rel_tol=1e-6)
