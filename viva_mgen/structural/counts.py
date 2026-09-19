"""Copy-number providers for structural-model ingredients.

Two providers produce ``{ProtID: int}`` copy-number maps (S1 ``ProtID`` keys,
same key space as :mod:`viva_mgen.structural.maritan_tables`) for Task 4 to
attach to each structural ingredient:

- ``maritan_counts`` — BASELINE (WC-MG-flavoured) copy numbers.
  Step-0 discovery confirmed ``datasets/karr_reference_values.json`` holds
  only aggregate MASS FRACTIONS (protein/DNA/RNA share of dry mass) — no
  per-gene protein copy numbers live there. The baseline instead derives a
  per-monomer RELATIVE abundance from ``datasets/karr_gene_expression.csv``
  (``expression_mean``; where that is absent/zero, falls back to
  ``synthesis_rate * half_life_min``), then rescales those relative
  abundances so they sum to ``TOTAL_PROTEIN_MOLECULES`` (see basis note
  below) and rounds to int. This is a DOCUMENTED PROXY for WC-MG protein
  copy numbers, not a published per-protein count — see the fidelity
  caveat below.

- ``mgen_sim_counts`` — VARIANT copy numbers from this repo's own
  simulation: a completed ``viva_mgen`` run's ``protein_counts`` emitter
  store (per-gene ``map[float]`` deltas, see
  ``viva_mgen/processes/translation.py``), accumulated to final per-gene
  totals.

Both map gene-level counts onto S1 monomer ``ProtID``s via the Task-1
``GeneRow.protein_monomer`` join, then lift monomer counts to COMPLEX
counts (including complexes composed of other complexes, e.g. the DNA
polymerase holoenzyme) with ``complex_count``, which parses the S1
"Complex Biosynthesis" string (terms like ``(2.0)MG_003_MONOMER'`` joined
by ``+``) and returns the limiting-subunit stoichiometric count.

TOTAL_PROTEIN_MOLECULES basis: CALIBRATED to Maritan et al. 2022's reported
volume occupancy. Maritan Table 1 (Frame 149 s, the one-chromosome birth cell)
gives a protein volume fraction of 0.144 inside a 144.47 nm sphere. Setting the
cell-wide protein total to 39,000 makes this proxy's packed protein occupancy
(sum of monomer molecular volumes / cell volume; see
viva_mgen.structural.build.estimate_occupancy) reproduce that 0.144 figure. As
a sanity check this is the right order of magnitude for a minimal-genome cell
(~1/30-1/50 the volume of E. coli, whose ~2-4x10^6 protein molecules/cell,
BioNumbers BNID 100088, scale down to tens of thousands; cf. Maritan's ~21,000
monomers across compartments at Frame 149 s). It remains a documented,
occupancy-calibrated estimate, not a published per-protein copy-number total.

FIDELITY CAVEAT: ``maritan_counts`` is an EXPRESSION-DERIVED PROXY for WC-MG
protein copy numbers (relative abundances built from RNA-level expression
data, scaled to an assumed cell-wide protein total) — it is NOT the
published per-protein copy numbers from Karr et al. 2012's whole-cell model
output (no such per-gene table ships with this repo's datasets). Treat
absolute values as order-of-magnitude estimates; the relative ranking
across monomers is the more defensible signal.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

from viva_mgen.structural.maritan_tables import GeneRow, ProteinRow

_DATASETS_DIR = Path(__file__).resolve().parents[2] / "datasets"
_KARR_EXPRESSION_CSV = _DATASETS_DIR / "karr_gene_expression.csv"

# See "TOTAL_PROTEIN_MOLECULES basis" in the module docstring.
# 39,000 calibrates the packed protein occupancy to Maritan's reported 0.144.
TOTAL_PROTEIN_MOLECULES = 39_000

# S1 "Complex Biosynthesis" terms look like "(2.0)MG_003_MONOMER'" (or, in the
# real supplement export, "(2.0)MG_003_MONOMER'+ '(2.0)MG_004_MONOMER'" — a
# stray leading/trailing "'" around each term after splitting on "+").
_BIOSYNTHESIS_TERM_RE = re.compile(r"\(([\d.]+)\)\s*'?([A-Za-z0-9_]+)'?")


def _to_float(raw) -> float:
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return 0.0
    return value if math.isfinite(value) else 0.0


def _gene_to_monomer(genes: list[GeneRow]) -> dict[str, str]:
    """``gene_id`` (``"MG_001"``) -> monomer ``ProtID`` (``"MG_001_MONOMER"``)."""
    return {g.gene_id: g.protein_monomer for g in genes if g.protein_monomer}


def _load_expression_relative() -> dict[str, float]:
    """``gene_id`` -> relative expression proxy from the karr expression table.

    Uses ``expression_mean``; where that is absent/zero/NaN, falls back to
    ``synthesis_rate * half_life_min``. Genes where both are unusable map to 0.
    """
    out: dict[str, float] = {}
    with _KARR_EXPRESSION_CSV.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            value = _to_float(row.get("expression_mean"))
            if not value:
                synth = _to_float(row.get("synthesis_rate"))
                half_life = _to_float(row.get("half_life_min"))
                value = synth * half_life if synth and half_life else 0.0
            out[row["gene_id"]] = max(value, 0.0)
    return out


def complex_count(biosynthesis: str, monomer_counts: dict[str, int]) -> int:
    """Limiting-subunit copy count for a complex from its S1 "Complex
    Biosynthesis" string, e.g. ``"(2.0)MG_003_MONOMER'+(2.0)MG_004_MONOMER'"``.

    Returns ``min(monomer_counts[subunit] // stoich)`` over every subunit
    term; ``0`` if the string has no parseable terms, any stoichiometry is
    non-positive, or any subunit is missing (or has count 0) in
    ``monomer_counts``.
    """
    terms = _BIOSYNTHESIS_TERM_RE.findall(biosynthesis or "")
    if not terms:
        return 0
    limiting = None
    for stoich_str, subunit in terms:
        stoich = float(stoich_str)
        if stoich <= 0:
            return 0
        available = int(monomer_counts.get(subunit, 0) // stoich)
        if limiting is None or available < limiting:
            limiting = available
    return max(limiting, 0)


def _resolve_complexes(proteins: list[ProteinRow], counts: dict[str, int]) -> None:
    """Lift ``counts`` (monomer counts already populated) to cover every
    complex in ``proteins``, in place. Handles complexes whose subunits are
    themselves complexes (e.g. the DNA polymerase holoenzyme) by resolving
    in dependency order; a subunit id that never appears in ``proteins`` at
    all counts as missing (0) rather than blocking resolution forever.
    """
    complexes = [p for p in proteins if p.kind != "monomer"]
    all_prot_ids = {p.prot_id for p in proteins}
    pending = {p.prot_id: p.biosynthesis for p in complexes}
    changed = True
    while pending and changed:
        changed = False
        for prot_id, biosynthesis in list(pending.items()):
            subunits = [s for _, s in _BIOSYNTHESIS_TERM_RE.findall(biosynthesis or "")]
            if all(s in counts or s not in all_prot_ids for s in subunits):
                counts[prot_id] = complex_count(biosynthesis, counts)
                del pending[prot_id]
                changed = True
    for prot_id in pending:  # unresolved (cyclic) references: no defensible count
        counts[prot_id] = 0


def maritan_counts(proteins: list[ProteinRow], genes: list[GeneRow]) -> dict[str, int]:
    """Baseline (WC-MG-flavoured) copy numbers for every ``ProtID`` in
    ``proteins`` — monomers from expression-derived relative abundances
    (see module docstring), complexes lifted via ``complex_count``.
    """
    gene_to_monomer = _gene_to_monomer(genes)
    expression = _load_expression_relative()

    monomer_relative: dict[str, float] = {}
    for gene_id, monomer_id in gene_to_monomer.items():
        rel = expression.get(gene_id, 0.0)
        monomer_relative[monomer_id] = monomer_relative.get(monomer_id, 0.0) + rel

    total_relative = sum(monomer_relative.values())
    scale = TOTAL_PROTEIN_MOLECULES / total_relative if total_relative > 0 else 0.0

    counts: dict[str, int] = {}
    for protein in proteins:
        if protein.kind == "monomer":
            rel = monomer_relative.get(protein.prot_id, 0.0)
            counts[protein.prot_id] = int(round(rel * scale))

    _resolve_complexes(proteins, counts)
    return counts


def _load_protein_counts_from_run_dir(run_dir) -> dict[str, float]:
    """Best-effort accumulation of a completed run's ``protein_counts``
    emitter store from its on-disk (Parquet) output under ``run_dir``.

    Not exercised by unit tests: no live ``viva_mgen`` run fixture exists
    in this repo yet, and the brief explicitly scopes this provider's
    tests to the pre-loaded-mapping path. Sums the ``protein_counts.<gene>``
    (flattened) or ``protein_counts`` (nested dict) column, whichever is
    present, across every row of every ``*.parquet`` file found under
    ``run_dir``. Raises ``FileNotFoundError`` if no parquet files are found
    or no ``protein_counts`` data is present, with a message pointing
    callers at the pre-loaded-mapping path instead.
    """
    import pandas as pd

    run_dir = Path(run_dir)
    parquet_files = sorted(run_dir.rglob("*.parquet"))
    if not parquet_files:
        raise FileNotFoundError(
            f"no *.parquet files under {run_dir}; pass an already-loaded "
            "{gene_id: count} mapping instead"
        )

    totals: dict[str, float] = {}
    found_protein_counts = False
    for path in parquet_files:
        frame = pd.read_parquet(path)
        nested_col = "protein_counts" if "protein_counts" in frame.columns else None
        flattened_cols = [c for c in frame.columns if c.startswith("protein_counts.")]
        if nested_col:
            found_protein_counts = True
            for row in frame[nested_col]:
                if not row:
                    continue
                for gene_id, delta in dict(row).items():
                    totals[gene_id] = totals.get(gene_id, 0.0) + _to_float(delta)
        elif flattened_cols:
            found_protein_counts = True
            for col in flattened_cols:
                gene_id = col[len("protein_counts."):]
                totals[gene_id] = totals.get(gene_id, 0.0) + frame[col].apply(_to_float).sum()

    if not found_protein_counts:
        raise FileNotFoundError(
            f"no protein_counts column found under {run_dir}; pass an "
            "already-loaded {gene_id: count} mapping instead"
        )
    return totals


def mgen_sim_counts(
    run_dir_or_state, proteins: list[ProteinRow], genes: list[GeneRow]
) -> dict[str, int]:
    """Variant (this repo's own simulation) copy numbers.

    ``run_dir_or_state`` accepts either:

    - a ``dict`` of already-accumulated final per-gene protein counts
      (``{gene_id: count}``, e.g. summed ``protein_counts`` emitter deltas
      from a completed run) — the tested path; or
    - a path (``str``/``Path``) to a completed run's directory, best-effort
      loaded via ``_load_protein_counts_from_run_dir`` (untested — see its
      docstring).

    Final per-gene totals are mapped gene -> monomer ``ProtID`` via the
    Task-1 join, rounded to int, then lifted to complexes with
    ``complex_count``.
    """
    if isinstance(run_dir_or_state, dict):
        gene_totals = run_dir_or_state
    else:
        gene_totals = _load_protein_counts_from_run_dir(run_dir_or_state)

    gene_to_monomer = _gene_to_monomer(genes)
    counts: dict[str, int] = {}
    for gene_id, total in gene_totals.items():
        monomer_id = gene_to_monomer.get(gene_id)
        if not monomer_id:
            continue
        counts[monomer_id] = counts.get(monomer_id, 0) + int(round(_to_float(total)))

    _resolve_complexes(proteins, counts)
    return counts
