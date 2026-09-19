import os
from pathlib import Path

import pytest
from viva_mgen.structural.build import mgen_ingredients, mgen_capsule, mgen_chromosome, build_mgen_pack
from viva_mgen.structural.counts import maritan_counts
from viva_mgen.structural.maritan_tables import load_proteins, load_genes


def _counts():
    return maritan_counts(load_proteins(), load_genes())


def test_ingredients_have_structure_or_sphere_and_region():
    ings = mgen_ingredients(_counts(), top_n=20)
    assert len(ings) == 20
    for ing in ings:
        assert ing.count > 0
        assert ing.region in ("interior", "surface", "fiber")
        assert (ing.structure is not None) or (ing.sphere_radius is not None)


def test_capsule_is_maritan_sphere():
    # Faithful to Maritan 2022 Table 1 (Frame 149 s): a ~144.47 nm sphere,
    # modelled as a NEAR-spherical spherocylinder (small half_len, aspect ~1)
    # whose volume conserves that sphere's volume.
    import math
    from viva_mgen.structural.build import _capsule_volume_a3, MGEN_MARITAN_SPHERE_RADIUS_A
    cap = mgen_capsule()
    assert cap.radius > 0
    aspect = (2 * cap.half_len + 2 * cap.radius) / (2 * cap.radius)
    assert 1.0 <= aspect < 1.3          # near-spherical, not elongated
    sphere_v = (4.0 / 3.0) * math.pi * MGEN_MARITAN_SPHERE_RADIUS_A ** 3
    assert _capsule_volume_a3(cap) == pytest.approx(sphere_v, rel=1e-3)


def test_chromosome_scale():
    chrom = mgen_chromosome()
    assert chrom.n_chromosomes == 1 and chrom.beads > 1000
    assert chrom.genome_csv and Path(chrom.genome_csv).exists()


def test_protein_occupancy_matches_maritan():
    # The copy-number total is calibrated so the packed protein volume occupancy
    # reproduces Maritan 2022 Table 1's reported protein volume fraction (0.144,
    # curated recipe, Frame 149 s). Guard the calibration.
    from viva_mgen.structural.build import estimate_occupancy
    occ = estimate_occupancy(maritan_counts(load_proteins(), load_genes()))
    assert occ["occupancy"] == pytest.approx(0.144, abs=0.02)


@pytest.mark.skipif(not (os.environ.get("PARSIMONY_HOME") or os.environ.get("PARSIMONY_BIN")),
                    reason="parsimony binary not configured")
def test_small_end_to_end_pack(tmp_path):
    res = build_mgen_pack(_counts(), out_dir=tmp_path, top_n=8, name="mgen_test")
    assert res["n_placed"] > 0
    assert (tmp_path).exists()
