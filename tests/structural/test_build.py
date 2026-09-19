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


def test_capsule_and_chromosome_scale():
    cap = mgen_capsule()
    assert cap.radius > 0 and cap.half_len >= cap.radius
    chrom = mgen_chromosome()
    assert chrom.n_chromosomes == 1 and chrom.beads > 1000
    assert chrom.genome_csv and Path(chrom.genome_csv).exists()


@pytest.mark.skipif(not (os.environ.get("PARSIMONY_HOME") or os.environ.get("PARSIMONY_BIN")),
                    reason="parsimony binary not configured")
def test_small_end_to_end_pack(tmp_path):
    res = build_mgen_pack(_counts(), out_dir=tmp_path, top_n=8, name="mgen_test")
    assert res["n_placed"] > 0
    assert (tmp_path).exists()
