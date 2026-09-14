from viva_mgen.expression_defaults import reference_protein_counts


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
    # zero out ALL metabolic-enzyme proteins -> gated reactions throttle to floor(0)
    starved = p.update({"nutrient_scale": 1.0, "protein_counts": {}}, 1.0)
    assert starved["growth_fraction"] < full["growth_fraction"]
