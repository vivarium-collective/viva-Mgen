# Maritan et al. 2022 — Building Structural Models of a Whole Mycoplasma Cell

- **BibTeX key:** `maritan2022structural`
- **PDF:** `papers/Maritan2022_structural_models_whole_mycoplasma.pdf`
- **Supplementary tables:** `papers/Maritan2022_supplementary/`
  - `S1_MG_proteins_ingredients.xlsx` — all protein ingredients (monomers +
    complexes): type, MW, monomer length, DNA-interaction info, complex
    biosynthesis, functional annotation, compartment, chosen structural
    representation (curated vs automated).
  - `S2_MG_genes.xlsx` — every MG gene: genomic coordinates, length, direction,
    type (+ codon/amino acid for tRNAs), transcription unit, **essentiality**,
    functional annotation, and (for protein-coding genes) UniProt ref + monomer
    name. Directly comparable to `datasets/genes.csv` and the `fig6` essentiality
    reference set.
  - `S3_membrane_protein_assignment.xlsx` — protein-monomer → compartment
    assignments in the WC-MG and cytoplasmic models, plus localization/signal-
    peptide predictions (BUSCA, SignalP-5.0, SOSUI, Phobius, PSORTb).
- **Cite:** Maritan, Autin, Karr, Covert, Olson & Goodsell. *J. Mol. Biol.* **434** (2022) 167351.
- **DOI:** [10.1016/j.jmb.2021.167351](https://doi.org/10.1016/j.jmb.2021.167351)

## Why it's here

The **structural companion** to the Karr 2012 whole-cell model this workspace
reproduces. Karr and Covert are co-authors; the 3D model is built directly on
the WC-MG model's molecular content, abundances, and localization (pulled from
WholeCellKB / WholeCellSimDB / WholeCellViz). It shows what the same data the
`viva_mgen` processes consume looks like as an integrated 3D cell.

## What it does

Presents the **first 3D structural model of an entire *M. genitalium* cell**
(3D-WC-MG), built with the **CellPACK** suite (Mesoscope → CellPACKgpu),
including all MG proteins, DNA, RNA, and the cell membrane at defined
life-cycle time points.

- Genome: circular ~580,076 bp encoding **525 genes** → **482 mRNA, 3 rRNA,
  4 sRNA, 36 tRNA** — same gene inventory the reproduction's expression panel
  draws on (cf. `datasets/genes.csv`, 525 genes).
- Nucleoid built with **LatticeNucleoid** (10 bp/bead DNA, nascent RNA, free
  RNA, RNA polymerase, ribosomes, nucleoid-associated proteins).
- Ingredient structures from a mix of **experimental + homology models**
  (AlphaFold2 noted as promising but, at time of writing, effective mainly for
  well-folded monomers — <half the ingredient list; 216 monomers, 116
  homooligomers, 150 heteromeric assemblies).
- Two data-gathering workflows compared: a **manually-curated** vs an
  **automated homology-driven** pipeline; both adequate for mesoscale
  properties (crowding, volume occupancy).
- Model-quality analysis estimates the **regularization** needed before these
  serve as starting points for atomic molecular-dynamics simulations.

## Relevance to this workspace

- Independent, orthogonal validation that the WC-MG molecular inventory (genes,
  transcript/protein counts, complex composition) is self-consistent enough to
  pack a whole cell — supports the fidelity claims of `fig1-architecture`,
  `fig3-expression`.
- Potential source for **crowding / volume-occupancy** acceptance bands if the
  reproduction ever adds a spatial or macromolecular-density observable.
