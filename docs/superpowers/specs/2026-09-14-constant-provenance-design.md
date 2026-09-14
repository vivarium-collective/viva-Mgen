# Constant provenance + supplement sourcing for viva-Mgen (closing gap #5)

Date: 2026-09-14
Status: design, pre-implementation

## Problem

Some reduced-process kinetic constants are order-of-magnitude physiological values,
not fitted Karr values, because the KB `parameters.json` carries **empty** param
dicts for 11 processes (DNADamage, HostInteraction, MacromolecularComplexation,
ProteinActivation, ProteinFolding, ProteinModification, RNAModification,
RibosomeAssembly, TerminalOrganelleAssembly, TranscriptionalRegulation,
tRNAAminoacylation). Gap #5 asks: source real values from the Karr supplement /
literature where they exist; where they genuinely don't, document each as
irreducible. Today the fidelity status is scattered across per-process docstrings;
there is no single, auditable record of which constants are real vs. estimated.

## Key ruling

Deliver **constant provenance** as the honest closure: a single, auditable record
classifying every kinetic constant in the reduced processes by source tier, plus a
sourcing pass against the accessible Karr supplement
(`workspace/references/papers/Karr2012_supplementary_information.pdf`) that either
wires a real value or records why the constant is irreducible from available
sources. This makes the reproduction's constant-fidelity fully transparent — which
IS the gap-#5 deliverable ("a sourced constants table + wiring, or a documented
'irreducible' note per constant"). New literature values that require external
access we cannot fabricate; those constants are marked `irreducible` with the
specific value/table that would be needed.

Source tiers:
- `real_kb` — value comes from `datasets/karr_parameters.json` (fitted Karr constant).
- `real_supplement` — value sourced from the Karr 2012 supplement (cite section/table).
- `order_of_magnitude` — a physiological estimate; the KB/supplement lacks it.
- `irreducible` — genuinely not available from the KB or the accessible supplement; note what external source would supply it.

## What v1 builds

### 1. Provenance audit → `datasets/constant_provenance.json` + `viva_mgen/provenance.py`

`viva_mgen/provenance.py` with `audit_constants() -> dict`: import every process
class (via `all_process_classes()` incl. the allocator), read each `config_schema`
entry whose default is a finite number (a kinetic constant), and record
`{process: {constant: {value, source_tier, kb_match (bool), note}}}`. `kb_match`
is computed by comparing the config default to any numeric value in that process's
`karr_parameters.json` dict (a value present there ⇒ `real_kb`). Non-`real_kb`
numeric constants default to `order_of_magnitude` and are refined by the curated
overrides in step 2.

`scripts/extract_constant_provenance.py` writes `datasets/constant_provenance.json`
= `{"source": ..., "n_constants": N, "by_tier": {...counts...}, "constants": <audit>}`.

### 2. Curated overrides + supplement sourcing (`viva_mgen/provenance.py` `_CURATED`)

A curated dict in `provenance.py` overrides the auto-classification for constants
where the truth is known from this project's earlier faithfulness work + a
supplement read:
- Constants confirmed sourced from the supplement get `real_supplement` + a citation note.
- Constants confirmed unavailable get `irreducible` + a note naming what's missing
  (e.g. "ProteinFolding per-protein chaperone rate matrix — Karr KB dict empty;
  supplement gives no per-protein folding rate; would need the source Minton/
  chaperone-kinetics dataset").
- The rest stay `order_of_magnitude` with a one-line justification.
The supplement PDF is read (Task 2) to fill/confirm these; any clearly-stated
value is ALSO wired into the corresponding process `config_schema` default (with a
comment citing the supplement) and its tier set to `real_supplement`.

### 3. `docs/CONSTANT_PROVENANCE.md`

A human-readable summary generated-or-written from the audit: a table per source
tier (counts + the constants in each), so a reader sees at a glance how much of the
kinetic parameterization is real-KB vs supplement vs estimated vs irreducible. This
is the transparency artifact.

### 4. Resources + fidelity note

Add `datasets/constant_provenance.json` to the Resources tab. Note in each
`order_of_magnitude`/`irreducible` process's fidelity line that its provenance is
recorded in the constant-provenance audit (avoid restating per-constant).

## Tests (`tests/test_provenance.py`)

- `audit_constants()` covers every process (keys == all_process_classes names that have numeric config constants); each entry has value/source_tier/kb_match/note.
- `real_kb` classification is correct for a known case: e.g. DNASupercoiling `gyrase_rate=1.2` matches `karr_parameters.json` ⇒ `real_kb`; ProteinFolding `spontaneous_rate` has no KB entry ⇒ not `real_kb`.
- Every constant has a source_tier in the allowed set; no `real_kb` constant is also curated `irreducible` (consistency).
- The generated `datasets/constant_provenance.json` exists, is valid JSON, and its `by_tier` counts sum to `n_constants`.

## Out of scope (v1)

- Fabricating literature kinetic values we cannot access (marked `irreducible` instead).
- Re-fitting constants; changing simulation behavior (except wiring a genuinely-sourced supplement value, which improves fidelity and is noted).

## Fidelity impact

Every kinetic constant in the reduced model gains an explicit, auditable provenance
(`real_kb` / `real_supplement` / `order_of_magnitude` / `irreducible`), and any
value the accessible Karr supplement provides for a currently-estimated constant is
sourced and wired. This closes gap #5 in the honest sense: the reproduction now
states exactly which constants are real and which are estimates, with the
irreducible ones documented and cited. `FIDELITY_GAPS.md` gap #5 marked done.
