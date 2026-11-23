# Study protocol — v1.8 (15 Oct 2026)

**Working title:** Evidence-Graded Molecular Benchmarks for Colorectal Cancer-on-a-Chip Models Derived from Single-Cell and Spatial Transcriptomics

> **Amendment to v1.7.** Adds a **fifth compartment (lymphoid / CD8 T cells)** as an additive extension. All sections of v1.7 and earlier remain in force unchanged. **Nothing in the existing four-compartment analysis is re-run, re-tuned or re-frozen.** `frozen_panel_v2.rds` (SHA-256 `a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb`) is not modified; the lymphoid panel is frozen to a **separate** file with its own hash, so every existing number remains byte-identical and reproducible.
>
> **§§14.1–14.8 were written and committed before any CD8 target was computed and before any CD8 analysis touched GSE132465.** §14.9 is the execution record, added afterwards.

## §14 (new) — Lymphoid compartment

### 14.1 Rationale

The four-compartment design of §2 was fixed by correspondence to the chip's four cell sources. It omits the lymphoid compartment, which is a major constituent of the colorectal tumour microenvironment and is of particular relevance to mismatch-repair-deficient tumours (§2 as amended in v1.7), where cytotoxic T-cell infiltration is a dominant feature. The omission is a scope limitation of the chip, not a judgement that the biology is unimportant. This amendment tests whether the framework extends to that compartment.

### 14.2 Cell definition

**Discovery (GSE178341).** CD8 T cells are cells from **tumour specimens** (`SPECIMEN_TYPE == "T"`) with `clMidwayPr == "TCD8"`. γδ T cells (`Tgd`), NK cells, ILC, `TZBTB16` and CD4 T cells are **excluded**, and are retained only as context populations for the spatial reference, where they already appear within the existing T/NK/ILC reference type.

Pre-analysis counts (label-level only, no expression read): 20,860 CD8 cells in tumour specimens; **62 donors with ≥1 CD8 cell, 58 donors with ≥20**, median 199 cells per donor. Donor power therefore exceeds that of the existing fibroblast and endothelial compartments.

**Validation (GSE132465), fixed before the validation model was run.** `Cell_type == "T cells" & Cell_subtype == "CD8+ T cells"`, within tumour-class, QC-passing cells. The cohort's T-cell annotation separates CD4, Treg, Th17, Tfh, γδ, NK and Unknown as distinct subtypes, all of which are excluded, so the validation compartment matches the discovery definition. Label-level counts before modelling: 23 donors with ≥1 CD8 cell, **20 with ≥20**, median 113.

### 14.3 Programs (pre-specified, fixed before any target was computed)

| Program | Compartment | Gene set | Set size | Genes blocked by existing evaluation pools | Free |
|---|---|---|---|---|---|
| `LYM_CYTOTOX` | Lymphoid | `GOBP_T_CELL_MEDIATED_CYTOTOXICITY` | 60 | 4 (7%) | 56 |
| `LYM_IL2_STAT5` | Lymphoid | `HALLMARK_IL2_STAT5_SIGNALING` | 199 | 50 (25%) | 149 |

Source: MSigDB via `msigdbr`, same version as recorded in `02_metadata/genesets/geneset_version.csv`. Two programs, consistent with the §9 limit of at most two per compartment.

**Conventional marker baseline for `LYM_IL2_STAT5`, pre-specified before modelling:** IL2RA, IL2RB, CISH, SOCS1, SOCS2, PIM1, MYC, BCL2L1. Recorded in `04_processed_data/programs/conventional_markers_lymphoid.csv`.

**Candidates considered and rejected, recorded for transparency:** `HALLMARK_ALLOGRAFT_REJECTION` (200 genes, 24% blocked) was rejected not for overlap with existing evaluation pools — an overlap check showed it is no worse than `IL2_STAT5` — but because its content is itself largely cytotoxic-T-cell biology, making it near-redundant with `LYM_CYTOTOX`. `REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL` (191 genes, 6% blocked) was rejected as too broad a match for "cytotoxic effector", being largely KIR/HLA/CD interaction molecules. No exhaustion program is defined, as no adequately covered MSigDB set exists for it.

### 14.4 Rules carried over unchanged from v1.0–v1.7

All thresholds and criteria are inherited, not re-set:

- Minimum **20 cells** per donor × specimen × compartment pseudobulk profile.
- Program **retained only if ≥0.80** of its genes are expressed in the lymphoid pseudobulk (`COVER_MIN`), with ≥10 evaluation and ≥5 predictor genes. **A program failing this gate is dropped, as `STR_ECM` was.**
- Evaluation / predictor gene partition: non-overlapping, ~50/50, with genes that are evaluation genes of any existing program excluded from the lymphoid predictor pool.
- Discovery: relaxed elastic net (α 0.5), 1–4 genes per program, OLS refit, donor-grouped 5-fold × 10-repeat cross-validation balanced by 10x chemistry, all filtering and scaling inside training folds.
- Baselines: discovery-cohort mean, equal-weight marker score, conventional markers, larger model.
- Primary metric: MAE ratio against the discovery-mean baseline, with correlation and calibration as supporting measures.
- Evidence categories: the v1.3 rule, unaltered. **Spatial evidence cannot promote a program that failed independent single-cell validation.**

### 14.5 Freezing and independence

The lymphoid panel is written to `04_processed_data/frozen_panel/frozen_panel_LYM_v1.rds` with its own recorded SHA-256. Scripts verifying `frozen_panel_v2.rds` continue to do so unchanged.

**Disclosure of partial blinding.** The lymphoid panel was selected entirely in the discovery cohort (GSE178341) and was never tuned, filtered or re-selected using GSE132465, so its external validation is formally independent. However, GSE132465 had already been used to validate seven programs under the four-compartment analysis, and the investigators were therefore not naive to that cohort's general behaviour. This constitutes a loss of *soft* blinding and **must be stated explicitly in the Methods**. The lymphoid result is to be reported as independent validation with this qualification attached, and not as a pre-registered first look.

### 14.6 Spatial assessment

No new deconvolution is required: T/NK/ILC is already one of the nine cell types in the RCTD reference, and criterion A already adjusts for its proportion. Criteria A and B are applied to the lymphoid programs exactly as specified in v1.5 §10d, including the ≥2/3-of-patients rule and the gene-coverage evaluability rule.

### 14.7 Chip correspondence — decision deferred, must be fixed before Table 2

The chip configuration of §2 contains **no lymphoid compartment**; THP-1 is monocytic, not lymphoid. Two framings are permissible and exactly one must be chosen before Table 2 is finalised:

- **(a) Chip gains a lymphoid compartment** (primary PBMC-derived CD8 T cells; Jurkat is not an acceptable substitute). Table 2 gains a fifth block with the usual interpretation limits.
- **(b) Chip does not gain one.** The lymphoid result is reported as a **documented gap in the chip**: a reproducible patient program that the current four-cell configuration cannot represent at all. Table 2 records the compartment with chip correspondence "none in the configuration considered".

This choice affects presentation only. It has no effect on any analysis, and the analysis proceeds before it is made.

### 14.8 Anticipated outcome, stated in advance

The lymphoid programs may fail validation, as four of the seven existing programs did. **A negative or Exploratory result is to be reported in full and is not grounds for revising the gene sets, the cell definition or the thresholds.** Any such revision would require a further untouched cohort before the word "validation" could be used again.

---

## §14.9 Execution record (added 15 Oct 2026, after running §§14.2–14.5)

### Pseudobulk

Built by `05_analysis_scripts/29_lymphoid_pseudobulk.R`, mirroring `08_GSE178341_pseudobulk.R` exactly (same 3-MAD-per-`batchID` QC on log₁₀ counts, log₁₀ genes and mitochondrial fraction; same `primary_fraction` restriction to unsorted/LiveMACS channels; same aggregation). The script re-derived QC independently and **verified that per-channel removal counts reproduce script 08 exactly across all 128 channels**, which establishes that the lymphoid compartment was built under identical QC rather than merely similar QC.

Eligible donors after QC:

| Minimum cells | `Lymphoid_CD8` (primary fraction) | `Lymphoid_CD8_allfrac` (sensitivity) |
|---|---|---|
| 10 | 56 | 60 |
| **20 (primary)** | **53** | 57 |
| 30 | 51 | 57 |
| 50 | 43 | 53 |

Lymphoid pseudobulk hash: `8d1f2c92aa37a5d1aac8637e937b1f09`.

### Coverage gate — outcome

Applied by `05_analysis_scripts/30_lymphoid_coverage.R` using the script 09 rule verbatim (CPM ≥ 1 in ≥25% of eligible donors; 53 donors).

| Program | Gene set | Set size | Expressed | Coverage | Retained |
|---|---|---|---|---|---|
| `LYM_CYTOTOX` | `GOBP_T_CELL_MEDIATED_CYTOTOXICITY` | 59 | 46 | **0.78** | **FALSE** |
| `LYM_IL2_STAT5` | `HALLMARK_IL2_STAT5_SIGNALING` | 200 | 168 | 0.84 | TRUE |

**`LYM_CYTOTOX` is dropped** under the §14.4 coverage rule, two genes short of the threshold (48 of 59 required). This is the same disposition applied to `STR_ECM` (coverage 0.774) in the original four-compartment analysis.

**No substitute gene set was sought.** `REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL` (191 genes) would plausibly have cleared the threshold, and was explicitly rejected at pre-specification. Adopting it after observing that the pre-specified set failed would constitute selecting a gene set on the basis of whether it passes, and is forbidden by §14.8. This decision is recorded here precisely because the alternative was available and tempting.

**Interpretation of the failure, for the Discussion.** Cytotoxic effector transcripts are sparsely detected in donor-level CD8 pseudobulk under the standing expression criterion. The lymphoid compartment is therefore represented in this study by an *activation and survival signalling* program, not by a cytotoxic effector program, and no claim about cytotoxicity is supported by these data.

### Gene partition (script 31)

Coverage 0.84; **84 evaluation genes, 65 predictor genes**; 19 candidate predictors removed because they are evaluation genes of an existing program. Lymphoid partition hash `a3235d265b2eaad9335322541c41dc34`. The original `gene_partition.tsv` was verified unchanged by recomputing its recorded hash. **None of the 26 frozen four-compartment panel genes appear in the lymphoid evaluation pool.**

### Discovery cross-validation (script 32)

53 donors, 5-fold × 10 repeats. Final panel: **CTSZ, IKZF4, P2RX4, TNFRSF4**. Lymphoid frozen panel SHA-256 `7dc3256d97477a796fcacd1df71021740888f7978a01ea36caf70787da7f513f`; `frozen_panel_v2.rds` verified unchanged at the end of the run.

| Method | MAE | MAE vs mean | Spearman | Calibration slope |
|---|---|---|---|---|
| `panel_refit` | 0.220 | **0.826** | 0.63 | 0.70 |
| `panel_shrunk` | 0.249 | 0.938 | 0.57 | 3.82 |
| `larger_model` | 0.155 | 0.582 | 0.84 | 1.05 |
| `equal_weight` | 0.212 | 0.799 | 0.63 | 0.75 |
| `conventional` | 0.274 | 1.032 | 0.08 | 0.26 |
| `mean` | 0.266 | 1.000 | — | — |

**Two points to carry into the Results, both unfavourable and both to be reported.** (i) The equal-weight baseline (0.799) slightly outperforms the fitted panel (0.826). The same ordering already occurs in three of the seven existing programs (apoptosis, hypoxia, TGF-β), though with smaller margins (~0.007 against 0.027 here); the general implication is that the panels' value lies in *which* genes are selected rather than in the fitted weights. (ii) The uncapped larger model (0.582) is substantially better than the 4-gene panel, so the cost of compression is larger in this compartment than elsewhere.

The conventional-marker failure (1.032, ρ 0.08) is consistent with, not exceptional to, this compartment: conventional markers also perform at or near the baseline for angiogenesis (1.018), apoptosis (0.966), TGF-β (0.965) and hypoxia (0.896), with IFN-γ (0.538) the only clear exception.

### Independent validation in GSE132465 (script 33)

4,410 CD8 cells; **20 donors** with ≥20 cells; all 4 panel genes present; evaluation-gene coverage 0.988.

| Scaling | n | MAE ratio | 95% CI | Spearman | Calibration slope | Status |
|---|---|---|---|---|---|---|
| **Frozen (primary)** | 20 | **0.646** | 0.47–0.868 | 0.579 | 0.751 | **Supported** |
| Within-cohort (secondary) | 20 | 1.03 | 0.50–2.15 | 0.43 | 0.559 | Not supported |

`LYM_IL2_STAT5` therefore meets the criterion for independent single-cell support. Its MAE ratio ranks second of the eight programs, behind `END_TNFA_NFKB` (0.44) and ahead of `EPI_HYPOXIA` (0.81). The within-cohort scaling result is Not supported, which is the same profile as both other independently supported programs (`END_TNFA_NFKB` and `EPI_HYPOXIA` are also FALSE on within-scaling), so it does not alter the category.

**Caveats for the write-up:** 20 donors is the smallest validation set in the study, and the upper confidence bound (0.868) lies closer to 1 than that of `END_TNFA_NFKB`.

### Status

`LYM_IL2_STAT5` is **independently supported in single-cell data**, pending spatial assessment under §14.6, which determines whether it is categorised with or without spatial evidence.

## Change log (continued from v1.7)

| Date | Change | Reason |
|---|---|---|
| 2026-10-15 | v1.8: lymphoid (CD8) compartment added as an additive fifth compartment, with two pre-specified programs, a separate frozen panel, inherited thresholds, and a disclosure of partial blinding with respect to GSE132465 | The four-compartment design omitted a major TME compartment, an omission made more visible by the v1.7 MSI/MSS framing. Feasibility confirmed before pre-specification: 58 donors with ≥20 CD8 cells, candidate gene sets largely unblocked by existing evaluation pools. Existing analysis untouched |
| 2026-10-15 | §14.9 execution record: pseudobulk built with QC verified identical to script 08 (128 channels), 53 eligible donors; `LYM_CYTOTOX` dropped at coverage 0.78; `LYM_IL2_STAT5` retained at 0.84; no substitute gene set sought | Records the coverage outcome and the decision not to substitute, under §14.8 |
| 2026-10-15 | §14.9 extended: partition (84 evaluation / 65 predictor genes), discovery CV (panel CTSZ, IKZF4, P2RX4, TNFRSF4; MAE ratio 0.826), and independent validation in GSE132465 (MAE ratio 0.646, 0.47–0.868, Supported, 20 donors) | Records the decisive pre-specified results, including the two unfavourable baseline comparisons |
