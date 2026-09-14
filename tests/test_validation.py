import math
from viva_mgen.validation import reference, emergent_macro_fractions, mrna_cv, score

def test_reference_entries_valid():
    for vid in ("mass_fraction_protein", "single_cell_mrna_cv", "doubling_time_h"):
        r = reference(vid)
        assert set(r) >= {"description", "target", "band", "unit", "source", "observable"}
        assert r["band"][0] <= r["band"][1]

def test_emergent_macro_fractions_renormalizes():
    row = {"emergent_mass_fractions": {"protein": 0.6, "DNA": 0.2, "RNA": 0.1, "metabolite": 0.1}}
    f = emergent_macro_fractions(row)  # renormalize over protein+DNA+RNA (drop metabolite)
    assert abs(sum(f.values()) - 1.0) < 1e-9
    assert abs(f["protein"] - 0.6/0.9) < 1e-9 and "metabolite" not in f

def test_emergent_macro_fractions_empty():
    assert emergent_macro_fractions({"emergent_mass_fractions": {}}) == {"protein": 0.0, "DNA": 0.0, "RNA": 0.0}

def test_mrna_cv():
    assert mrna_cv([100.0, 100.0, 100.0]) == 0.0
    cv = mrna_cv([80.0, 100.0, 120.0])
    assert abs(cv - (20.0*math.sqrt(2/3)) / 100.0) < 1e-9   # population std / mean

def test_score():
    r = reference("doubling_time_h")  # target 9, band [7,11]
    assert score(9.0, r)["pass"] is True and score(9.0, r)["closeness"] == 1.0
    assert score(20.0, r)["pass"] is False and score(20.0, r)["closeness"] < 0.5

def test_emergent_macro_fractions_matches_fig2_run_wiring():
    # A tiny stand-in for the fig2-growth run.py's `rows[-1]` (the final emitter
    # row) — fig2 calls emergent_macro_fractions(rows[-1]) directly, so this
    # exercises the same code path the run.py wiring uses, without running the
    # full ~9 h simulation.
    rows = [
        {"emergent_mass_fractions": {"protein": 0.10, "DNA": 0.03, "RNA": 0.02, "metabolite": 0.05}},
        {"emergent_mass_fractions": {"protein": 0.62, "DNA": 0.169, "RNA": 0.093, "metabolite": 0.20}},
    ]
    f = emergent_macro_fractions(rows[-1])
    total = 0.62 + 0.169 + 0.093
    assert abs(f["protein"] - 0.62 / total) < 1e-9
    assert abs(f["DNA"] - 0.169 / total) < 1e-9
    assert abs(f["RNA"] - 0.093 / total) < 1e-9
    assert abs(sum(f.values()) - 1.0) < 1e-9
    # within the fig2 study.yaml bands (secondary gates)
    assert 0.55 <= f["protein"] <= 0.80
    assert 0.10 <= f["DNA"] <= 0.35
    assert 0.05 <= f["RNA"] <= 0.25
