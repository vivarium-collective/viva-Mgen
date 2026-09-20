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

from pbg_parsimony import Capsule, Chromosome, Ingredient, StructureRef, build_pack

from viva_mgen.constants import GENOME_LENGTH_BP
from viva_mgen.structural.maritan_tables import ProteinRow, load_genes, load_proteins
from viva_mgen.structural.structures import structure_ref_for, uniprot_index

DATA = Path(__file__).resolve().parent / "data"

# --- cell geometry -----------------------------------------------------------
# FAITHFUL TO MARITAN. Biological M. genitalium is flask-shaped (an ovoid body
# with a tapered terminal organelle), but Maritan et al. 2022 -- the model this
# investigation reproduces -- state explicitly: "The cell shape is approximated
# as a sphere, and the attachment organelle is omitted", modelling the first
# cell-cycle phase (one chromosome present, ~23% of the cycle). To reproduce
# Maritan faithfully we use their sphere, NOT a biological flask.
#
# Maritan Table 1, Frame 149 s (the very start of the cell cycle = a one-
# chromosome birth cell): cell radius 144.47 nm. The parsimony packer will not
# fill a zero-length cylinder (a true half_len=0 sphere packs nothing), so the
# cell is modelled as a NEAR-spherical spherocylinder with a small medial
# half-length, its cap radius solved to CONSERVE Maritan's sphere volume
# (aspect ratio ~1.1:1). This keeps the packed volume -- and therefore the
# volume-occupancy comparison below -- honest against Maritan's 0.144 protein
# fraction.
MGEN_MARITAN_SPHERE_RADIUS_A = 1444.7   # 144.47 nm, Maritan Table 1 Frame 149 s
_MGEN_HALF_LEN_A = 150.0                 # small medial length so the packer fills it

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
    """Maritan's spherical cell (Table 1, Frame 149 s: radius 144.47 nm),
    modelled as a near-spherical spherocylinder with a small medial half-length
    and cap radius solved to conserve the sphere's volume (a true half_len=0
    sphere packs nothing). Faithful to Maritan's stated sphere approximation."""
    import math
    R = MGEN_MARITAN_SPHERE_RADIUS_A
    target_v = (4.0 / 3.0) * math.pi * R ** 3
    hl = _MGEN_HALF_LEN_A
    # Solve pi*r^2*(2*hl) + (4/3)*pi*r^3 = target_v for r by bisection.
    lo, hi = 0.0, R
    for _ in range(60):
        r = 0.5 * (lo + hi)
        v = math.pi * r * r * (2.0 * hl) + (4.0 / 3.0) * math.pi * r ** 3
        if v < target_v:
            lo = r
        else:
            hi = r
    return Capsule(half_len=hl, radius=0.5 * (lo + hi))


def _capsule_volume_a3(cap: Capsule) -> float:
    """Spherocylinder volume in Angstrom^3: cylinder (length 2*half_len) + two
    hemispherical caps."""
    import math
    r = cap.radius
    return math.pi * r * r * (2.0 * cap.half_len) + (4.0 / 3.0) * math.pi * r ** 3


def _protein_volume_a3(protein) -> float:
    """Rough molecular volume (Angstrom^3) from residue count via the standard
    protein specific volume: MW ~= 110 Da/residue, V ~= 1.21 A^3/Da
    (Harpaz/Gerstein). Returns 0 when the residue count is unknown (complexes
    without a seq_length; their subunits are counted separately)."""
    aa = getattr(protein, "seq_length", None)
    if not aa:
        return 0.0
    return float(aa) * 110.0 * 1.21


def estimate_occupancy(counts, top_n=None) -> dict:
    """Estimate macromolecular volume occupancy (crowding) of the packed cell:
    sum of placed monomers' molecular volumes divided by the capsule volume.

    Monomers carry a residue count (seq_length) → a volume; complexes are
    covered through their subunit monomers (counted at their own abundance), so
    to avoid double counting only monomer ingredients contribute here. Returns
    ``{occupancy, macromol_volume_a3, cell_volume_a3, monomers_counted}`` — a
    labeled estimate, not a mesh-exact figure. Typical bacterial cytoplasmic
    crowding is a volume fraction of ~0.2-0.4."""
    proteins = load_proteins()
    by_id = {p.prot_id: p for p in proteins}
    cell_v = _capsule_volume_a3(mgen_capsule())
    macromol_v = 0.0
    counted = 0
    for prot_id, n in counts.items():
        p = by_id.get(prot_id)
        if p is None or getattr(p, "kind", "") != "monomer":
            continue
        v = _protein_volume_a3(p)
        if v > 0 and n > 0:
            macromol_v += v * float(n)
            counted += 1
    return {
        "occupancy": (macromol_v / cell_v) if cell_v else 0.0,
        "macromol_volume_a3": macromol_v,
        "cell_volume_a3": cell_v,
        "monomers_counted": counted,
    }


def mgen_chromosome() -> Chromosome:
    """Single circular supercoiled chromosome, RENDERED as a visible nucleoid.

    Bead count is sized to GENOME_LENGTH_BP; ``segment`` gives each bead a real
    dsDNA mesh (RCSB 1BNA, B-DNA) so the fiber actually draws (without it the
    chromosome is defined but invisible), tinted tan and coiled by ``supercoil``.
    ``genome_csv`` seats RNAP at real (abundance-weighted) transcription sites
    instead of uniformly along the fiber. Bead geometry (spacing 135 A ~= 40 bp,
    radius 12 A) follows ecoli_3d's convention."""
    beads = max(1, round(GENOME_LENGTH_BP / _BP_PER_BEAD))
    return Chromosome(
        beads=beads, spacing=135.0, bead_radius=12.0, n_chromosomes=1,
        genome_csv=_pbg_genome_csv(),
        segment=StructureRef("pdb", "1BNA"),
        supercoil={"radius": 90.0, "pitch": 130.0, "domains": 200},
        color=(0.85, 0.75, 0.45),
    )


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
    region, principal_vector = _placement(protein)
    return Ingredient(
        id=protein.prot_id,
        count=count,
        structure=structure,
        sphere_radius=sphere_radius,
        color=_color_for(protein.function),
        region=region,
        compartment="cytoplasm",
        display_name=protein.name,
        category=protein.function,
        principal_vector=principal_vector,
    )


def _placement(protein: ProteinRow) -> tuple[str, tuple | None]:
    """Positional class + orientation for a protein, following Maritan's spatial
    model (their Table 1: cytoplasm / DNA-bound / membrane / extracellular):

    * DNA-binding proteins (S1 ``dna_binding`` = dsDNA/ssDNA) -> ``region="fiber"``
      so the packer seats them ON the supercoiled nucleoid at binding sites,
      not free in the bulk (Maritan's "NAP" class).
    * Membrane (compartment ``m``) -> ``region="surface"`` with a
      ``principal_vector`` so the packer orients them to the membrane normal.
    * Extracellular (compartment ``e``) -> ``region="surface"`` (unoriented): the
      single-membrane capsule has no true exterior, so surface-seat them (mostly
      surface adhesins) rather than mis-placing them in the cytoplasm.
    * Everything else (cytoplasm ``c``) -> ``region="interior"``.
    """
    if protein.dna_binding:
        return "fiber", None
    if protein.compartment == "m":
        return "surface", (0.0, 0.0, 1.0)   # orient TM axis to the membrane normal
    if protein.compartment == "e":
        return "surface", None
    return "interior", None


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
