from viva_mgen.expression_defaults import reference_protein_counts
from viva_mgen.kb import normalize_gene_id


def test_reference_protein_counts_positive():
    ref = reference_protein_counts()
    assert isinstance(ref, dict) and len(ref) > 100
    assert all(v >= 0 for v in ref.values())
    assert sum(1 for v in ref.values() if v > 0) > 100   # most expressed genes have a positive ss


def _metab(config):
    from viva_mgen.processes.metabolism import MetabolismFbaReproductionProcess
    from viva_mgen.core import build_core
    return MetabolismFbaReproductionProcess(config=config, core=build_core())


def test_coupling_off_matches_wildtype():
    p = _metab({})   # default enzyme_coupling False
    base = p.update({"nutrient_scale": 1.0}, 1.0)
    # protein_counts present but coupling off -> ignored
    same = p.update({"nutrient_scale": 1.0, "protein_counts": {"x": 0.0}}, 1.0)
    assert abs(base["growth_fraction"] - same["growth_fraction"]) < 1e-9


def test_coupling_on_reference_is_wildtype():
    ref = reference_protein_counts()
    p = _metab({"enzyme_coupling": True, "reference_protein_counts": ref})
    off = _metab({})
    on = p.update({"nutrient_scale": 1.0, "protein_counts": dict(ref)}, 1.0)
    base = off.update({"nutrient_scale": 1.0}, 1.0)
    assert on["growth_fraction"] > 0
    assert abs(on["growth_fraction"] - base["growth_fraction"]) < 0.05   # ~wild-type at reference


def test_coupling_on_enzyme_knockdown_lowers_growth():
    ref = reference_protein_counts()
    p = _metab({"enzyme_coupling": True, "reference_protein_counts": ref})
    full = p.update({"nutrient_scale": 1.0, "protein_counts": dict(ref)}, 1.0)
    # iPS189 has enough flux slack that a uniform 50%/25% knockdown of every
    # gated reaction's bound doesn't bind the FBA optimum at all (it's flat
    # down to ~4% of reference) -- so use two lower, still-uniform-across-all-
    # genes knockdown levels (2% and 1% of reference) that DO bind, and check
    # a strictly graded (3-point monotonic) response. A single-reaction
    # artifact (as in the original bug, where only GLYC3Pabc got zeroed) would
    # not produce a smooth multi-point gradient like this -- it would only
    # ever show a step down to the same floor value.
    mid = p.update({"nutrient_scale": 1.0, "protein_counts": {g: v * 0.02 for g, v in ref.items()}}, 1.0)
    low = p.update({"nutrient_scale": 1.0, "protein_counts": {g: v * 0.01 for g, v in ref.items()}}, 1.0)
    # zero out ALL metabolic-enzyme proteins -> gated reactions throttle to floor(0)
    starved = p.update({"nutrient_scale": 1.0, "protein_counts": {}}, 1.0)
    assert starved["growth_fraction"] < low["growth_fraction"] < mid["growth_fraction"] < full["growth_fraction"]


def test_coupling_reference_covers_most_model_genes():
    # non-regression for the symbol->locus-tag bridge: reference_protein_counts
    # is keyed by panel key (symbol-or-id) while model.genes are locus tags
    # (e.g. "MG005") — the process must re-key so the two actually overlap.
    ref = reference_protein_counts()
    p = _metab({"enzyme_coupling": True, "reference_protein_counts": ref})
    model_gene_ids = {normalize_gene_id(g.id) for g in p._model.genes}
    overlap = model_gene_ids & set(p._ref_norm.keys())
    assert len(overlap) >= 80   # most of the 126 model genes should be gate-able
