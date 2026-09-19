"""Assemble the mgen ingredient roster + cell geometry and pack it with
:mod:`pbg_parsimony` (the reusable parsimony-engine bridge).

The organism-specific half of the structural pipeline: turn a
``{ProtID: count}`` copy-number map (Task 3, ``viva_mgen.structural.counts``)
into a list of :class:`pbg_parsimony.Ingredient`, size a spherocylinder cell
(:class:`pbg_parsimony.Capsule`) and a single circular chromosome
(:class:`pbg_parsimony.Chromosome`), and hand all three to
:func:`pbg_parsimony.build_pack`.
"""

from __future__ import annotations

import csv
import math
import tempfile
from pathlib import Path

from pbg_parsimony import Capsule, Chromosome, Ingredient, build_pack

from viva_mgen.constants import GENOME_LENGTH_BP
from viva_mgen.structural.maritan_tables import ProteinRow, load_genes, load_proteins
from viva_mgen.structural.structures import structure_ref_for, uniprot_index

DATA = Path(__file__).resolve().parent / "data"

# --- cell geometry -----------------------------------------------------------
# M. genitalium is among the smallest known free-living cells: Maritan et al.
# 2022 and the Karr et al. 2012 whole-cell model both describe a
# roughly-spherical/pleomorphic cell body ~0.3-0.4 um in diameter.
# MGEN_RADIUS_UM is the spherocylinder cap radius (~half that diameter).
MGEN_RADIUS_UM = 0.20

# MGEN_VOLUME_FL: viva_mgen.constants.CELL_INITIAL_DRY_WEIGHT_FG is 3.93 fg at
# ~30% dry-mass fraction (FRACTION_WET_WEIGHT=0.7) -> a wet weight ~13 fg;
# at a typical cytoplasmic density of ~1.1 g/mL that implies a cell volume on
# the order of a few hundredths of a fL. A bare sphere at MGEN_RADIUS_UM is
# (4/3)*pi*(0.20 um)^3 ~= 0.034 fL; MGEN_VOLUME_FL below is a round
# literature/order-of-magnitude figure (~2x that bare-sphere volume, allowing
# a modest capsule elongation) consistent with both estimates. Refine against
# Maritan et al. 2022's stated cell dimensions if/when available.
MGEN_VOLUME_FL = 0.067

# Genome contour -> bead count, following ecoli_3d's GENOME_BEADS convention
# (34,000 beads for a 4,641,652 bp E. coli genome, i.e. ~135 bp/bead).
_BP_PER_BEAD = 135.0

# --- genome-annotation CSV for pbg_parsimony -------------------------------
# pbg_parsimony's Rust engine (parsimony-core's genome.rs, ``Genome::from_csv``)
# PARSES ``Chromosome.genome_csv`` (Python only ``shutil.copy``s it through
# unread; the parsing happens downstream, in the packer binary) to seat RNAP
# at real, abundance-weighted transcription sites. It expects a specific,
# positionally-parsed layout that is NOT what mgen_genes.csv (Task 1's S2
# reader) provides:
#   line 1: ``# genome_length_bp=<N>``
#   header: ``old_locus_tag,locus_tag,start,end,strand,biotype``
#   rows:   old_locus_tag (matched against ``MG_<digits>`` tokens embedded in
#           ingredient ids, e.g. ``MG_003_MONOMER`` -> ``MG_003``), locus_tag
#           (informational; duplicated here, no second identifier available),
#           1-based start/end (bp), strand (``+``/``-``), biotype.
# This reshapes ``load_genes()`` (coord/length/direction/name) into that exact
# layout rather than passing mgen_genes.csv verbatim, which would fail to
# parse (no ``genome_length_bp`` comment, wrong column count/positions) and
# make the packer silently skip genome-driven placement.
_GENOME_CSV_HEADER = ("old_locus_tag", "locus_tag", "start", "end", "strand", "biotype")
_genome_csv_cache: str | None = None


def _pbg_genome_csv() -> str:
    """Derive pbg_parsimony's expected genome-annotation CSV from this
    workspace's own S2 gene table (``load_genes()``); write once per process
    to a cache file and reuse the path on subsequent calls."""
    global _genome_csv_cache
    if _genome_csv_cache is not None and Path(_genome_csv_cache).exists():
        return _genome_csv_cache

    genes = load_genes()
    fd, path = tempfile.mkstemp(prefix="mgen_genome_pbg_", suffix=".csv")
    with open(fd, "w", newline="") as f:
        f.write(f"# genome_length_bp={GENOME_LENGTH_BP}\n")
        writer = csv.writer(f)
        writer.writerow(_GENOME_CSV_HEADER)
        for gene in genes:
            if not gene.coord or not gene.length:
                continue
            start = int(gene.coord)
            end = start + int(gene.length) - 1
            if end <= start:
                continue
            strand = "+" if gene.direction == "Forward" else "-"
            biotype = "protein_coding" if gene.protein_monomer else gene.gtype
            writer.writerow([gene.gene_id, gene.gene_id, start, end, strand, biotype])

    _genome_csv_cache = path
    return path


def mgen_capsule() -> Capsule:
    """Spherocylinder sized to MGEN_VOLUME_FL at cap radius MGEN_RADIUS_UM."""
    return Capsule.from_volume_fl(MGEN_VOLUME_FL, radius_um=MGEN_RADIUS_UM)


def mgen_chromosome() -> Chromosome:
    """Single circular chromosome; bead count sized to GENOME_LENGTH_BP, with
    a genome-annotation CSV so the pack seats RNAP at real (abundance-
    weighted) transcription sites instead of uniformly along the fiber."""
    beads = max(1, round(GENOME_LENGTH_BP / _BP_PER_BEAD))
    return Chromosome(beads=beads, n_chromosomes=1, genome_csv=_pbg_genome_csv())


# --- molecular-weight sphere-radius proxy ------------------------------------
# For proteins with neither a curated PDB id nor a resolvable UniProt
# accession (no all-atom structure), approximate as a sphere from molecular
# weight: MW (Da) from sequence length at the standard average residue mass
# (~110 Da/aa), treated as a sphere at globular-protein density
# (~1.35 g/cm^3; e.g. Erickson 2009, Biol. Proced. Online 11:32).
#   mass (g)      = MW / Avogadro
#   volume (cm^3) = mass / density
#   radius (cm)   = (3 * volume / (4 * pi)) ** (1/3)
_AA_MW_DA = 110.0
_PROTEIN_DENSITY_G_PER_CM3 = 1.35
_AVOGADRO = 6.02214076e23
_CM_TO_ANGSTROM = 1e8


def _sphere_radius_angstrom(seq_length: int) -> float:
    mw_da = seq_length * _AA_MW_DA
    mass_g = mw_da / _AVOGADRO
    volume_cm3 = mass_g / _PROTEIN_DENSITY_G_PER_CM3
    radius_cm = (3.0 * volume_cm3 / (4.0 * math.pi)) ** (1.0 / 3.0)
    return radius_cm * _CM_TO_ANGSTROM


# Deterministic palette keyed by the S1 "function" column (used verbatim as
# Ingredient.category — see mgen_ingredients). Not exhaustive: unrecognized
# categories fall back to a neutral gray.
_CATEGORY_COLOR = {
    "DNA replication/maintenance": (0.85, 0.30, 0.30),
    "RNA synthesis/maturation": (0.30, 0.75, 0.35),
    "transcription": (0.25, 0.55, 0.85),
    "translation": (0.65, 0.35, 0.85),
    "protein folding/maturation": (0.90, 0.60, 0.20),
    "protein transport/singaling": (0.20, 0.70, 0.70),
    "metabolism": (0.55, 0.55, 0.55),
    "lipoprotein": (0.85, 0.75, 0.30),
    "cytokinesis/motility": (0.90, 0.40, 0.70),
    "host cell interaction": (0.40, 0.20, 0.70),
    "MG-specific": (0.35, 0.65, 0.45),
    "uncharacterized": (0.60, 0.60, 0.60),
}
_DEFAULT_COLOR = (0.60, 0.60, 0.60)


def _color_for(category: str) -> tuple:
    return _CATEGORY_COLOR.get(category, _DEFAULT_COLOR)


def _to_ingredient(protein: ProteinRow, count: int, uniprot_by_prot: dict) -> Ingredient | None:
    structure = structure_ref_for(protein, uniprot_by_prot)
    sphere_radius = None
    if structure is None:
        if not protein.seq_length:
            return None  # no PDB/UniProt and no MW basis -> skip (brief policy)
        sphere_radius = _sphere_radius_angstrom(protein.seq_length)
    return Ingredient(
        id=protein.prot_id,
        count=count,
        structure=structure,
        sphere_radius=sphere_radius,
        color=_color_for(protein.function),
        region="surface" if protein.compartment == "m" else "interior",
        compartment="cytoplasm",
        display_name=protein.name,
        category=protein.function,
    )


def mgen_ingredients(counts: dict, top_n: int | None = None) -> list:
    """Build the placeable-ingredient roster from a ``{ProtID: count}`` map.

    Only ``ProtID``s with ``counts[prot_id] > 0`` are considered; if ``top_n``
    is given, only the ``top_n`` highest-count of those are placed. Proteins
    with no resolvable structure AND no sequence length (no MW basis for a
    sphere proxy) are skipped.
    """
    proteins = load_proteins()
    genes = load_genes()
    uniprot_by_prot = uniprot_index(genes, proteins)  # built once, not per-protein
    by_id = {p.prot_id: p for p in proteins}

    present = [(prot_id, count) for prot_id, count in counts.items()
               if count > 0 and prot_id in by_id]
    present.sort(key=lambda kv: kv[1], reverse=True)
    if top_n is not None:
        present = present[:top_n]

    ingredients = []
    for prot_id, count in present:
        ing = _to_ingredient(by_id[prot_id], count, uniprot_by_prot)
        if ing is not None:
            ingredients.append(ing)
    return ingredients


def build_mgen_pack(counts: dict, *, out_dir, top_n: int | None = None, name: str = "mgen") -> dict:
    """Assemble ingredients + geometry and pack the cell via ``pbg_parsimony``.

    Single-membrane cell (M. genitalium has no outer membrane) -> no
    ``envelope=`` (the ``build_pack`` default, ``None``).
    """
    ingredients = mgen_ingredients(counts, top_n=top_n)
    capsule = mgen_capsule()
    chromosome = mgen_chromosome()
    return build_pack(ingredients, capsule, chromosome, out_dir=out_dir, name=name, envelope=None)
