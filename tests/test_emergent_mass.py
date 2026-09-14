import math
from viva_mgen.mass_composition import (
    rna_mass_g, protein_mass_g, dna_mass_g, metabolite_mass_g, emergent_composition,
    _N_A,
)

def test_rna_mass_one_gene():
    # one 300-nt transcript, count 1 → 300*340/N_A g
    g = rna_mass_g({"g1": 1.0}, lengths={"g1": 300.0})
    assert math.isclose(g, 300*340.0/_N_A, rel_tol=1e-9)

def test_protein_mass_one_gene():
    # count 2, 300-nt ORF → 300/3=100 aa each → 2*100*110/N_A
    g = protein_mass_g({"g1": 2.0}, lengths={"g1": 300.0})
    assert math.isclose(g, 2*100*110.0/_N_A, rel_tol=1e-9)

def test_dna_mass_matches_karr_fraction():
    # chromosome_copy 1 → ~0.6-0.7 fg (Karr DNA fraction 0.1688 * 3.93 fg ~= 0.66 fg)
    fg = dna_mass_g(1.0) * 1e15
    assert 0.55 < fg < 0.75

def test_metabolite_mass_positive():
    assert metabolite_mass_g({"atp": 1e6, "gtp": 1e6}) > 0.0

def test_emergent_composition_protein_dominant_and_sums():
    comp = emergent_composition(
        rna_counts={"g1": 100.0}, protein_counts={"g1": 50000.0},
        chromosome_copy=1.0, pools={"atp": 1e5},
        lengths={"g1": 1000.0})
    assert comp["total"] > 0
    assert abs((comp["RNA"]+comp["protein"]+comp["DNA"]+comp["metabolite"]) - comp["total"]) < 1e-30
    assert comp["protein"] == max(comp["RNA"], comp["protein"], comp["DNA"], comp["metabolite"])


def test_composite_emits_emergent_macromolecule_mass():
    from viva_mgen.composites.mgen import build_mgen
    from viva_mgen.core import build_core
    from process_bigraph import Composite
    core = build_core()
    doc = build_mgen(core=core)
    comp = Composite({"state": doc}, core=core)
    comp.run(3.0)
    phys = comp.state["cell"]["physiology"]
    assert phys["emergent_mass"] > 0.0
    fr = phys["emergent_mass_fractions"]
    assert set(fr.keys()) == {"RNA", "protein", "DNA"}
    s = sum(fr.values())
    assert abs(s - 1.0) < 1e-6 or s == 0.0
    assert "metabolite" not in fr
