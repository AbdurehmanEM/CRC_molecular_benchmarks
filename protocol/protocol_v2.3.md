# Protocol amendment v2.3 — external comparator analysis and grading-rule sensitivity

**Date fixed:** 23 November 2026, after 73a (discovery only) and before 73b (external scoring).
**Comparator file:** `04_processed_data/frozen_panel/frozen_comparators_v2.3.rds`
**SHA-256:** `f9e6bb28187dfb275652f6c3e78cfb7441ee9a75592dca3b08a2e9d765b23b7f`
**Approved by:** Abdurehman Eshete Mohammed, Soni Shaikh, Yuksel Cetin (25 November 2026)

## Rationale
Bioinformatics Advances asks that new methods be compared with existing approaches on real data. The frozen panels were compared externally only with the discovery-mean baseline. Script 61 (alternative panels) tested how stable the panel genes are, not how the panel compares with a standard marker-selection method.

## Scripts
`05_analysis_scripts/revision/73a_freeze_comparators.R` (discovery only) and `73b_score_comparators.R` (external scoring). Both use `50_revision_config.R` and `50_revision_helpers.R`, and both require the script-52 reproduction gate to have passed (PASS, 11 Nov 2026). They also check the hashes of the frozen panels v2 (a14ccb45…) and LYM_v1 (7dc3256d…). 73a refuses to run if the comparator file already exists. 73b refuses to run unless the comparator hash matches.

## Pre-specified comparators (GSE178341 only)
- **Pool (both comparators).** The pool is the program's predictor genes from the gene partition that are present, by symbol, in every main validation cohort for that compartment. Only gene presence is used, never validation expression values.
- **C1_DE (standard single-cell marker route).** C1 uses donor-level pseudobulk edgeR quasi-likelihood on raw discovery counts. Each compartment is compared with all other profiled compartments (epithelial, stromal, myeloid, endothelial, CD8), with donor as a blocking factor, using eligible donors only (≥ 20 cells). C1 takes the 4 pool genes with the highest positive logFC. The number of those 4 that also meet FDR < 0.05 and logFC > 1 is reported. The strict filter is not used for selection because the pools are small (8 genes for END_ANGIO, 20 for STR_TGFB).
- **C2_CORR.** C2 takes the 4 pool genes with the highest Spearman correlation with the discovery target, using all discovery donors.
- **Score.** Each comparator's score is an equal-weight mean z-score using discovery means and SDs, calibrated by OLS on discovery donors. It is stored in the frozen-panel format and evaluated with the same code as the panels.

## 73a execution record (23 Nov 2026)
- **Discovery donors with ≥ 20 cells:** epithelial 60, stromal 27, myeloid 60, endothelial 40, CD8 53. The C1 matrix had 43,113 genes × 240 profiles from 60 donors.
- **Genes passing FDR < 0.05 and logFC > 1 (genome-wide):** epithelial 1,474, stromal 1,051, myeloid 1,343, endothelial 1,347, CD8 1,042.

| Program | Pool | C1_DE genes | C1 meeting FDR<0.05 & logFC>1 | C2_CORR genes |
|---|---|---|---|---|
| EPI_HYPOXIA | 70 | HOXB9, SLC2A1, CASP6, NEDD4L | 4 | GBE1, HDLBP, FOSL2, RBPJ |
| EPI_APOPTOSIS | 52 | ERBB3, CASP6, CCND1, EBP | 4 | CASP7, BCL10, ANXA1, CDK2 |
| STR_TGFB | 20 | TGIF1, TGFBR1, NCOR2, SMURF1 | 0 | SMURF1, SPTBN1, BMPR2, SMURF2 |
| MYE_IFNG | 74 | FGL2, HLA-DQA1, HLA-DRB1, HLA-DMA | 4 | STAT1, IFIT3, HLA-B, XAF1 |
| MYE_INFLAM | 60 | CLEC5A, MARCO, C5AR1, TLR2 | 4 | NAMPT, TLR2, LYN, C5AR1 |
| END_ANGIO | 8 | SLCO2A1, POSTN, MSX1 (only 3 with logFC > 0) | 3 | MSX1, TNFRSF21, SLCO2A1, LPL |
| END_TNFA_NFKB | 58 | PLPP3, FJX1, ACKR3, PLK2 | 4 | KDM6B, ETS2, TNFAIP3, B4GALT5 |
| LYM_IL2_STAT5 | 63 | ICOS, CST7, HOPX, BATF | 4 | CTSZ, IKZF4, PRKCH, AHR |

**Notes fixed before external scoring:**
1. C2 overlaps heavily with the frozen panels. It is identical in gene set for MYE_IFNG and EPI_APOPTOSIS, and shares 3 of 4 genes for MYE_INFLAM, STR_TGFB, END_ANGIO and END_TNFA_NFKB. For those programs, C2 against the panel mainly tests equal versus fitted weights, not gene choice. C1_DE is the primary contrast for gene choice.
2. For STR_TGFB, none of the C1 genes is strongly fibroblast-specific (0 of 4 meet the strict filter), because TGF-β pathway genes are broadly expressed.
3. END_ANGIO C1 has 3 genes, because only 3 of the 8 pool genes are enriched in endothelial cells. It is used as is and not topped up.

*Correction (25 Nov 2026): note 1 misstated the overlap. C2 is identical in gene set to the frozen panel for MYE_IFNG and MYE_INFLAM; it shares 3 of 4 genes for EPI_APOPTOSIS, STR_TGFB, END_ANGIO and END_TNFA_NFKB, 2 of 4 for LYM_IL2_STAT5 and none for EPI_HYPOXIA (Supplementary Table S8g, part C). The interpretation is unchanged.*

## Evaluation (73b)
- Cohorts, targets, frozen scaling, donor eligibility (≥ 10 donors) and the missing-gene rule are identical to the published analysis. 73b checks that the panel rows reproduce the script-52 ratios exactly. GSE144735 is descriptive only.
- For each program and cohort, 73b reports the MAE ratio against the discovery mean, a 95% donor bootstrap CI (2000 replicates, SEED 20261111), Spearman correlation and status, using the same rule as script 52.
- **Paired comparison.** 73b computes log MAE(panel) − log MAE(comparator) on the same bootstrap resamples, so values below 0 favour the panel.
- **Pooling.** DerSimonian–Laird pooling is applied to the main cohorts. Programs evaluable in only one cohort (k = 1) are reported with that cohort's bootstrap values.
- **Interpretation is fixed in advance.** The comparison is descriptive. It does not change any panel, evidence category or headline result. All rows are reported, including those where a comparator has the lower error.

## Grading-rule sensitivity (addresses the post hoc v1.3 rule; uses existing results only)
The final categories are re-derived under three alternative rules:
1. Single-cell support requires that the pooled upper 95% CI be < 1 **and** that at least 2 cohorts be individually Supported.
2. The within-cohort scaling (secondary) status is used in place of frozen scaling.
3. Spatial support requires criterion A alone, and separately criteria A and B with ≥ 3/4 of patients.

A table reports the category of each program under each rule. Any program whose category changes is stated in the Results.

## 73b execution record (23 Nov 2026)
The panel rows reproduced the script-52 values (check built into 73b). Pooled MAE ratios are against the discovery mean. "Panel/C1" and "Panel/C2" are MAE(panel)/MAE(comparator) on paired bootstrap resamples, where values below 1 favour the panel.

| Program | k (donors) | Panel | C1_DE | C2_CORR | Panel/C1 (95% CI) | Panel/C2 (95% CI) |
|---|---|---|---|---|---|---|
| END_TNFA_NFKB | 1 (11) | 0.444 | 0.994 | 0.557 | **0.447 (0.383–0.501)** | **0.798 (0.755–0.848)** |
| EPI_HYPOXIA | 3 (45) | 0.807 | 0.825 | 0.658 | 1.000 (0.856–1.168) | 1.185 (0.664–2.116) |
| EPI_APOPTOSIS | 3 (45) | 0.670 | 1.264 | 0.731 | 0.552 (0.255–1.198) | 0.910 (0.678–1.221) |
| STR_TGFB | 2 (23) | 0.749 | 0.952 | 0.916 | 0.805 (0.515–1.259) | 0.799 (0.620–1.028) |
| MYE_IFNG | 3 (44) | 1.621 | 0.985 | 1.567 | 1.629 (0.424–6.262) | 1.029 (0.946–1.120) |
| MYE_INFLAM | 3 (44) | 0.771 | 1.095 | 0.788 | 0.726 (0.515–1.023) | 0.988 (0.808–1.209) |
| END_ANGIO | 1 (11) | 2.113 | 1.171 | 1.547 | **1.804 (1.182–3.466)**, C1 lower | 1.366 (0.974–1.864) |
| LYM_IL2_STAT5 | 2 (32) | 0.377 | 0.940 | 1.269 | 0.401 (0.099–1.617) | 0.299 (0.031–2.863) |

LYM_IL2_STAT5: GSE166555 was added after the first 73b run, once `50c_save_lymphoid_validation.R` had saved its CD8 log-CPM (fixed 23 Nov: `on.exit()` reset the working directory when the script was sourced; the original is kept as `_before_wdfix`). The comparators were not changed. The panel's pooled ratio (0.377) reproduces the published 0.376 (0.128–1.111). The first run (GSE132465 only: panel 0.646, C1 0.800, panel/C1 0.808) is superseded.

Summary: the panel had the lower point error than C1 in 5 of 8 programs, equal error for hypoxia, and higher error for IFN-γ and angiogenesis. The CI excluded 1 only for END_TNFA_NFKB (favouring the panel over both comparators) and END_ANGIO (favouring C1, although C1 was itself above 1 at 1.171).

## Grading-rule sensitivity — execution record (23 Nov 2026, script 74)
The published rule reproduced Supplementary Table S11 for all 8 programs (gate built into script 74). The inputs were workbook V14 sheets S8, S8b, S10c and S11, and no new data were analysed.

| Program | Published | R1 (≥2 cohorts) | R2 (within-cohort) | R3a (A only) | R3b (≥3/4 patients) |
|---|---|---|---|---|---|
| EPI_HYPOXIA | Single-cell | Single-cell | **Exploratory** | **Top** | Single-cell |
| EPI_APOPTOSIS | Exploratory | Exploratory | **Top** | Exploratory | Exploratory |
| STR_TGFB | Exploratory | Exploratory | Exploratory | Exploratory | Exploratory |
| MYE_IFNG | Exploratory | Exploratory | **Top** | Exploratory | Exploratory |
| MYE_INFLAM | Exploratory | Exploratory | **Top** | Exploratory | Exploratory |
| END_ANGIO | Unsupported | Unsupported | Unsupported | Unsupported | Unsupported |
| END_TNFA_NFKB | Top | **Exploratory** | **Exploratory** | Top | Top |
| LYM_IL2_STAT5 | Exploratory | Exploratory | Exploratory | Exploratory | Exploratory |

Five of 8 programs change under at least one rule. The spatial variants change only hypoxia. Single-cohort support (R1) and the choice of scaling (R2) decide the TNF/NF-κB grade. TGF-β, CD8 IL2/STAT5 and angiogenesis are unchanged under every rule. Reported in Supplementary Table S11b.
