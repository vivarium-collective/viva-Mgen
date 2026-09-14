# Per-protein maturation routing for viva-Mgen (gap #4, phase 1: sub-item 4d + ProcessingI/Translocation)

Date: 2026-09-14
Status: design (autonomous ruling), pre-implementation

## Problem

Gap #4 is the cluster of within-submodel de-reductions. The protein-maturation
pipeline currently applies each step to ALL proteins in its input map, because the
per-protein classification wasn't carried:
- **ProcessingII** applies only the signal peptidase; the lipoprotein-specific
  diacylglyceryl-transferase (Lgt) step is skipped because no per-protein
  lipoprotein flag was available (this is FIDELITY_GAPS item 4d).
- **ProcessingI** applies N-terminal Met cleavage indiscriminately, though only a
  subset of proteins are actually Met-cleaved.
- **Translocation** would translocate any protein, though only secretory/lipoprotein
  (and membrane) proteins are translocated.

But the Karr KB **does** carry the classification per ProteinMonomer:
`signalSequenceType` ∈ {None (448), secretory (20), lipoprotein (14)} and
`nTerminalMethionineCleavage` ∈ {0 (447), 1 (35)}. So this is decodable, not
irreducible.

## Key ruling (autonomous)

Decode the real per-protein maturation classification from the KB and route each
maturation step to the correct protein subset:
- ProcessingI Met-cleavage → the 35 proteins flagged `nTerminalMethionineCleavage=1`.
- ProcessingII Lgt transferase → the 14 `lipoprotein` proteins (item 4d).
- Translocation → the `secretory` + `lipoprotein` proteins (the 34 with a signal sequence).

This replaces "applies to all / representative" with genuine KB-grounded routing.
Low figure risk: the maturation-form stores (nascent/process_i_done/translocated/
processed_ii) are a parallel representation and do not feed `protein_counts`
(which translation produces); restricting them changes the maturation
representation's fidelity, not the proteome counts the figures read. The other
gap-#4 sub-items (4a complexation network solve, 4b per-RNase RNA processing,
4c per-AA tRNA charging) remain staged in FIDELITY_GAPS.md.

## What v1 builds

### 1. Decode classification → dataset (`kb_decode.py` + `datasets/karr_protein_maturation.json`)

Add `decode_protein_maturation(mat) -> {gene_id: {signal_type, met_cleavage, signal_length}}`
to `kb_decode.py`: for each ProteinMonomer, read `signalSequenceType`
(str), `nTerminalMethionineCleavage` (0/1 → bool), `signalSequenceLength` (float),
and map its `MG_###_MONOMER` id to the gene id `MG_###`. A new
`scripts/extract_kb_maturation.py` writes `datasets/karr_protein_maturation.json`:
`{"maturation": {gene_id: {signal_type, met_cleavage, signal_length}}}` (only the
non-default entries need storing, but storing all decoded is fine).

Add `kb.load_protein_maturation() -> dict` (gene_id → classification), and a
helper `kb.maturation_by_panel_key() -> dict` keying it by the expression-panel
key (symbol-or-id, via `load_genes()`), so the protein processes — which operate on
panel-keyed maps — can look up each protein's class.

### 2. Route the maturation steps (`viva_mgen/processes/protein.py`)

- **ProteinProcessingI**: gains config `met_cleavage_genes` (default = the set of
  panel keys with `met_cleavage=True`, from the loader). The deformylase step still
  applies to all nascent proteins (all formyl-Met initiated); the **Met-cleavage**
  (aminopeptidase) throughput now applies only to flagged proteins — i.e. a protein
  NOT in the set is still deformylated and passed through, but is not counted
  against the aminopeptidase limit. Concretely: partition `nascent` into
  cleaved-subset and not; the aminopeptidase enzyme limit gates only the subset;
  both subsets move nascent→process_i_done (deformylation applies to all), so the
  overall pipeline still advances every protein. Keep the real KB deformylase 38 /
  aminopeptidase 6 rates.
- **ProteinProcessingII**: gains config `lipoprotein_genes` (default = the 14
  lipoprotein panel keys). Now applies the Lgt diacylglyceryl-transferase step to
  the lipoprotein subset (its slow KB rate 0.0165 s⁻¹ gating only those), while the
  signal peptidase (11 s⁻¹) applies to all translocated proteins as before. So
  lipoproteins are limited by min(peptidase, transferase-for-lipoproteins);
  non-lipoproteins only by the peptidase — the correct per-class behavior the
  earlier contract said it couldn't do.
- **ProteinTranslocation**: gains config `translocated_genes` (default = secretory +
  lipoprotein panel keys). Only these are translocation substrates; a
  non-secretory/non-lipoprotein protein in `process_i_done` is passed through
  untranslocated (moved to `translocated` store unchanged OR left — pick: pass
  through so the pipeline doesn't stall, but do not consume GTP for it). GTP is
  consumed only for the real translocation substrates.

Each config defaults to the decoded set (loaded lazily), so the routing is real by
default but overridable/testable.

### 3. Fidelity notes

Update the ProcessingI / ProcessingII / Translocation `description` fidelity lines:
Met-cleavage / lipoprotein-Lgt / translocation now act on the REAL per-protein
subsets decoded from the KB (`signalSequenceType`, `nTerminalMethionineCleavage`),
not on all proteins.

### 4. Resources

Add `datasets/karr_protein_maturation.json` to the dashboard Resources tab
(`data_sources.py`, knowledge-base category).

## Tests

- Unit: `decode_protein_maturation` / loaders — 14 lipoproteins, 20 secretory, 35
  met-cleavage (assert counts); MG_067 is a lipoprotein; keys map to panel keys.
- Process: ProcessingII with a lipoprotein vs non-lipoprotein input — the
  lipoprotein's throughput is limited by the slow transferase, the non-lipoprotein's
  is not (higher). ProcessingI: a met-cleavage protein vs not. Translocation: a
  secretory protein is translocated (consumes GTP), a cytoplasmic one is not.
- Non-regression: full suite passes; composite builds + runs.

## Out of scope (v1)

- 4a complexation network steady-state solve, 4b per-RNase RNA processing, 4c
  per-AA tRNA charging — staged.
- Membrane-protein topology / SRP-vs-Sec pathway distinction beyond the
  secretory/lipoprotein/none split.

## Fidelity impact

The protein-maturation pipeline now routes each step to the REAL protein subset
from the KB (14 lipoproteins get the Lgt transferase, 35 proteins get N-terminal
Met cleavage, 34 signal-sequence proteins get translocated), closing FIDELITY_GAPS
item 4d and refining ProcessingI/Translocation — replacing the "applies to all"
reduction with genuine per-protein classification.
