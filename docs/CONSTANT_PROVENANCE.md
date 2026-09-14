# Constant provenance (gap #5)

Every numeric kinetic constant declared in a process `config_schema` — 83
constants across the reduced-model submodels — is classified into one of four
source tiers by `viva_mgen.provenance.audit_constants()` and written to
`datasets/constant_provenance.json` by `scripts/extract_constant_provenance.py`.
This page is the transparency artifact: it shows how much of the
parameterization is real Karr 2012 knowledge-base (KB) data, how much is
explicitly stated in the Karr 2012 supplemental information (SI), how much is
a defensible order-of-magnitude physiological/engineering estimate, and how
much is genuinely irrecoverable from the sources available to this
reproduction.

## The four tiers

- **`real_kb`** — the constant's default is tied, by an explicit code-comment
  citation, to a specific value in the real Karr 2012 knowledge base
  (`datasets/karr_parameters.json`, consumed via `viva_mgen.kb`). Verified by
  `_REAL_KB` in `viva_mgen/provenance.py`: the (class, constant) pair names
  the exact KB key and expected value, and the constant is only tagged
  `real_kb` if the config default still matches it.
- **`real_supplement`** — the constant's value is explicitly stated in the
  Karr 2012 SI (the paper text, figures, or figure legends) and cited by
  section/figure. Curated in `_CURATED`.
- **`order_of_magnitude`** — no KB or SI value exists; the default is a
  defensible physiological estimate, stoichiometric fact (e.g. N ATP per
  reaction event), or engineering choice (e.g. a discretization resolution or
  a Hill-law softening of an originally Boolean rule).
- **`irreducible`** — the Karr KB parameter dict for that process is
  confirmed empty (`viva_mgen.kb.karr_process_params` returns `{}`), the SI is
  silent, and the constant stands in for a genuine per-species / per-enzyme /
  per-reaction rate **matrix** that a single reduced-model scalar structurally
  cannot recover — not merely an unmeasured number, but one for which no
  single "real" value exists to source without the full per-reaction Karr
  Data S1 kinetic-parameter tables (see "What's actually in the supplement
  file" below).

Only non-`real_kb` constants are curated by hand in `_CURATED`
(`viva_mgen/provenance.py`); `real_kb` is derived automatically from the
`_REAL_KB` map plus a numeric-match check against the current config default.

## What's actually in the supplement file

The only Karr 2012 supplement file present in this workspace is
`workspace/references/papers/Karr2012_supplementary_information.pdf` (16
pages). It was read in full for this audit. It is the *Cell* main-text
article plus two SI figures (S1, S2) and one page of "Extended Experimental
Procedures" — **not** the ~1,900-parameter appendix. That page states
explicitly:

> "Chapters 2 and 3 of Data S1 discuss the mathematics of each cellular
> variable and process submodel... Table S3 details the *M. genitalium*
> reconstruction..."

i.e. the actual per-constant kinetic-parameter tables live in a separate
document (**Data S1**) that is not included in this repository. Fetching it
was out of scope for this audit (see the Task 2 brief). Consequently, an
exhaustive read of the available 16 pages found exactly **one** explicit,
unambiguous numeric match to an existing constant:
`MassGrowthReproductionProcess.cell_cycle_length_s` — Figure 2A states "Mean
τ = 9.0 h" (32,400 s), which matches the existing default exactly.

Everything else in the visible text (kcat ranges for Tdk/Pdp in Figure 7,
ATP/GTP synthesis rates in Figure 5, collision frequencies in Figure 3, gene
essentiality counts in Figure 6, etc.) is either a *model-predicted result*
being validated against experiment rather than an input parameter, or too
ambiguous (e.g. a figure axis limit) to cite as an unambiguous source value
per the "do not invent values" rule — so none of it was used to reclassify
any other constant.

## By tier

### real_kb (11)

Tied to an explicit Karr KB value via `_REAL_KB` in `viva_mgen/provenance.py`.

| Process.constant | Value | KB source |
| --- | --- | --- |
| DNASupercoilingReproductionProcess.gyrase_rate | 1.2 | `gyraseActivityRate` |
| DNASupercoilingReproductionProcess.atp_per_act | 2.0 | `gyraseATPCost` |
| ProteinProcessingIReproductionProcess.deformylase_specific_rate | 38.0 | `deformylaseSpecificRate` |
| ProteinProcessingIReproductionProcess.aminopeptidase_specific_rate | 6.0 | `methionineAminoPeptidaseSpecificRate` |
| ProteinProcessingIIReproductionProcess.peptidase_specific_rate | 11.0 | `lipoproteinSignalPeptidaseSpecificRate` |
| ProteinProcessingIIReproductionProcess.transferase_specific_rate | 0.0165 | `lipoproteinDiacylglycerylTransferaseSpecificRate` |
| ProteinTranslocationReproductionProcess.translocase_specific_rate | 2.71e12 | `translocaseSpecificRate` |
| ProteinTranslocationReproductionProcess.gtp_per_monomer | 2.0 | `SRP_GTPUsedPerMonomer` |
| FtsZPolymerizationReproductionProcess.activation_rate | 1.1 | `activationFwd` |
| FtsZPolymerizationReproductionProcess.dissociation_rate | 0.01 | `activationRev` |
| ReplicationReproductionProcess.dna_pol_rate | 200.0 | 2x `dnaPolymeraseElongationRate` (100) |

### real_supplement (1)

Explicitly stated in the Karr 2012 SI main text/figures and cited by
section/figure.

| Process.constant | Value | SI citation |
| --- | --- | --- |
| MassGrowthReproductionProcess.cell_cycle_length_s | 32400.0 s (9.0 h) | Figure 2A: "Mean τ = 9.0 h" |

### irreducible (10)

Karr KB parameter dict confirmed empty for that process, SI silent, and the
constant stands in for a per-species/per-enzyme rate matrix the reduced model
structurally cannot recover from a single scalar. Each note names the
external source (the missing Data S1 chapter or an equivalent primary-data
table) that would be needed to source a real value.

| Process.constant | Value | What would be needed |
| --- | --- | --- |
| ProteinFoldingReproductionProcess.spontaneous_rate | 0.05 | Data S1 per-protein folding-kinetics table (the notFolding-set Gillespie rates) |
| ProteinFoldingReproductionProcess.chaperone_rate | 0.002 | Data S1 per-protein chaperone-substrate kinetics table |
| ProteinModificationReproductionProcess.modification_specific_rate | 6.0 | Data S1 per-enzyme modification-rate table |
| RNAModificationReproductionProcess.modification_rate | 0.4 | Data S1 RNA-modification kinetics table (13 enzymes x 86 modifications) |
| RNAModificationReproductionProcess.enzyme_kcat | 2.0 | Data S1 per-enzyme modification kcat table |
| TRNAAminoacylationReproductionProcess.synthetase_kcat | 20.0 | Data S1 synthetase-kinetics table (20 AA x 37 tRNA/tmRNA reactions) |
| DNADamageReproductionProcess.base_rate | 0.001 | Data S1 site-specific spontaneous lesion-rate table (vulnerable-motif sites) |
| DNADamageReproductionProcess.agent_rate | 0.01 | Data S1 per-agent dose-response damage-rate table |
| ProteinActivationReproductionProcess.default_k | 1.0 | Data S1 per-protein activation-rule (Boolean/metabolite-threshold) table |
| TerminalOrganelleAssemblyReproductionProcess.threshold | 1.0 | Data S1 per-reaction localization stoichiometry matrix |

### order_of_magnitude (61)

No KB or SI value exists. Each is either a defensible physiological/
stoichiometric estimate or an engineering choice (discretization resolution,
softened-Boolean modeling parameter). A subset carries a `_CURATED`
justification note (see `viva_mgen/provenance.py`) without changing tier;
the rest are plain physiological estimates.

| Process.constant | Value |
| --- | --- |
| AllocatorProcess.pool_cap | 1e9 |
| ChromosomeCondensationReproductionProcess.bind_rate | 0.001 |
| ChromosomeCondensationReproductionProcess.initial_fraction | 0.0 |
| ChromosomeDynamicsReproductionProcess.dna_pol_rate_bp | 20.0 |
| ChromosomeDynamicsReproductionProcess.genome_length_bp | 580070.0 |
| ChromosomeDynamicsReproductionProcess.n_bins | 580.0 |
| ChromosomeDynamicsReproductionProcess.n_dnab | 4.0 |
| ChromosomeDynamicsReproductionProcess.n_dnan | 4.0 |
| ChromosomeDynamicsReproductionProcess.n_gyrase | 20.0 |
| ChromosomeDynamicsReproductionProcess.n_rna_pol | 50.0 |
| ChromosomeDynamicsReproductionProcess.n_smc | 280.0 |
| ChromosomeDynamicsReproductionProcess.n_ssb | 30.0 |
| ChromosomeDynamicsReproductionProcess.n_tf | 12.0 |
| ChromosomeDynamicsReproductionProcess.n_topo | 10.0 |
| ChromosomeDynamicsReproductionProcess.rna_gene_span_bins | 20.0 |
| ChromosomeDynamicsReproductionProcess.rna_pol_rate_bp | 50.0 |
| ChromosomeDynamicsReproductionProcess.struct_bind_tau_s | 1200.0 |
| ChromosomeDynamicsReproductionProcess.struct_collision_frac | 0.12 |
| ChromosomeDynamicsReproductionProcess.struct_turnover_per_s | 0.008 |
| ChromosomeSegregationReproductionProcess.complete_threshold | 0.999 |
| ChromosomeSegregationReproductionProcess.initial_fraction | 0.0 |
| ChromosomeSegregationReproductionProcess.segregation_rate | 0.000258... |
| CytokinesisReproductionProcess.contraction_rate_nm_per_s | 0.0517... |
| CytokinesisReproductionProcess.initial_cell_width_nm | 200.0 |
| CytokinesisReproductionProcess.replication_complete_threshold | 0.999 |
| CytokinesisReproductionProcess.ring_subunits_for_full_rate | 320.0 |
| DNADamageReproductionProcess.n_bins | 580.0 |
| DNARepairReproductionProcess.atp_per_repair | 4.0 |
| DNARepairReproductionProcess.repair_rate | 0.01 |
| DNASupercoilingReproductionProcess.genome_length_bp | 580070.0 |
| DNASupercoilingReproductionProcess.initial_sigma | 0.0 |
| DNASupercoilingReproductionProcess.relaxed_bp_per_turn | 10.5 |
| DNASupercoilingReproductionProcess.setpoint | -0.06 |
| FtsZPolymerizationReproductionProcess.initial_ftsz_gtp | 0.0 |
| FtsZPolymerizationReproductionProcess.initial_ring | 0.0 |
| HostInteractionReproductionProcess.cooperativity | 4.0 |
| HostInteractionReproductionProcess.half_max_fraction | 0.7 |
| MassGrowthReproductionProcess.density_g_per_ml | 1.1 |
| MassGrowthReproductionProcess.initial_mass_fg | 3.93 |
| MetabolismFbaReproductionProcess.amino_acid_base_supply | 1e6 |
| MetabolismFbaReproductionProcess.enzyme_coupling_cap | 1.0 |
| MetabolismFbaReproductionProcess.enzyme_coupling_floor | 0.0 |
| MetabolismFbaReproductionProcess.ntp_base_supply | 1e6 |
| ProteinActivationReproductionProcess.hill_n | 2.0 |
| ProteinFoldingReproductionProcess.atp_per_fold | 7.0 |
| ProteinModificationReproductionProcess.atp_per_modification | 1.0 |
| RNAProcessingReproductionProcess.processing_rate | 3.168... |
| ReplicationInitiationReproductionProcess.complex_threshold | 30.0 |
| ReplicationInitiationReproductionProcess.initial_complex | 0.0 |
| ReplicationInitiationReproductionProcess.recruit_rate | 0.00231... |
| ReplicationReproductionProcess.dnaA_synthesis_rate | 0.00208... |
| ReplicationReproductionProcess.dnaA_threshold | 30.0 |
| ReplicationReproductionProcess.dntp_synthesis_rate | 37.25... |
| ReplicationReproductionProcess.genome_length_bp | 580070.0 |
| ReplicationReproductionProcess.init_synth_fraction | 0.08 |
| ReplicationReproductionProcess.initial_dnaA | 0.0 |
| ReplicationReproductionProcess.initial_dntp | 0.0 |
| RibosomeAssemblyReproductionProcess.gtp_per_complex | 2.0 |
| TRNAAminoacylationReproductionProcess.atp_per_charge | 1.0 |
| TranscriptionReproductionProcess.rna_pol_reference | 100.0 |
| TranslationReproductionProcess.gtp_per_protein | 600.0 |

## Reading this table

- Two of the brief's named empty-KB-dict processes —
  `TranscriptionalRegulationReproductionProcess` and
  `MacromolecularComplexationReproductionProcess` — have **no** numeric
  config constants at all (their config schemas hold only maps and a `seed`),
  so there's nothing for this audit to classify; their fidelity notes already
  state that their real KB content is the TF regulatory network / complex
  stoichiometry (both wired from the real KB elsewhere, see
  `datasets/karr_tf_regulation.json` and `datasets/karr_complexes.json`).
- The `irreducible` list is short (10 constants) relative to the 72
  originally order-of-magnitude constants: most order-of-magnitude constants
  in this reproduction are stoichiometric facts (N ATP/GTP per reaction
  event), initial conditions, or engineering choices (bin counts, pool caps)
  rather than missing rate measurements — those correctly stay
  `order_of_magnitude`, not `irreducible`.
- Regenerate this data (not this prose) with
  `PYTHONPATH=<worktree> <venv>/bin/python scripts/extract_constant_provenance.py`.
