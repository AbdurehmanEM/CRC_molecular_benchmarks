# Final code review — CRC Biomarker Project (19 Oct 2026)

**Scope.** All 49 scripts in `05_analysis_scripts/` (00–41, including 12/12b, 18b, 28 and the `_before_*` backups, which were not reviewed in detail because `00_run_all.R` ignores them). Every file parses without error under R 4.3. The revised `00_run_all.R` was run on a stub project, and the script-17 cell-selection logic was tested by simulation. The real data were not available, so numerical outputs were not re-run.

**Verdict.** The core statistical design is implemented correctly: raw-count pseudobulk, leakage control, frozen-panel validation, meta-analysis and spatial block bootstrap. No error was found that changes a reported single-cell validation number.

The problems found are of three kinds:
- **Pipeline breakers.** Five defects stop a clean re-run, and in one case destroy a script. A sixth item is a missing manual extraction step.
- **A file no script writes.** The evidence table is built from a file that no script generates.
- **One spatial-reference sampling bias.** It needs a diagnostic run before you know whether it matters.

Corrected scripts are listed at the end.

---

## A. Pipeline breakers (fixed)

| # | Script | Problem | Fix |
|---|---|---|---|
| A1 | `04_download_GSE178341_GSE144735.R` | **The file is not a download script.** It contains `setwd("F:/…")`, `file.create("05_analysis_scripts/06_audit_GSE178341.R")` and `file.edit(…)`. `file.create()` **empties an existing file**, so on a clean machine stage 4 deletes the contents of script 06 and downloads nothing. The original download code is lost, as happened with script 05 earlier. | Rebuilt from the file names used by scripts 06/07/08/12/17/29. If `02_metadata/*_checksums.csv` from the original download exist, it checks the files against them instead of overwriting. |
| A2 | 20, 21, 40 ← 23 | **No script writes `spatial_support_combined/spatial_support_summary_main.csv`**, but script 20 (evidence categories), 21 and 40 (figures) read it. The evidence table therefore depends on a file made by hand or by a deleted script, and a clean re-run fails at stage 24. | Script 23 now writes it from the `All spatial cohorts (main)` rows. Script 25 now compares the existing file against those rows. |
| A3 | 17 (+ setup) | Script 17 installs spacexr from **GitHub (dmcable/spacexr)**, which provides `create.RCTD()`/`run.RCTD()`. Scripts 18 and 22 call the **Bioconductor** API `createRctd()`/`runRctd()`. On a clean machine stage 18 fails. | GitHub install removed. Setup installs Bioconductor spacexr and checks that `createRctd` exists. |
| A4 | 25 | Hard-coded `setwd("F:/CRC_Biomarker_Project 2/…")`. On any other machine stage 26 errors and, by default, **stops the pipeline before all lymphoid stages (27–39)**. | Removed. Also: calibration-slope column name corrected (`cal_slope`, previously skipped silently); stability file made explicit (see B4). |
| A5 | 00_run_all → 12b | Stage 12 runs `12b_KUL_from_saved.R`, which reads `GSE144735_counts_sparse.rds`. Only script 12 creates that file, and script 12 is not in the stage table (it also ends right after loading). A clean re-run fails at stage 12. | 12b now builds the count matrix itself if it is missing (loader copied verbatim from script 12, so results are unchanged). Script 12 can be removed from the archive. |
| A6 | 17, 22 | The Valdeolivas archives must be extracted to `04_processed_data/Valdeolivas_Visium/unzipped/SN*/`, and the GSE280318 position files to `04_processed_data/GSE280318/`. **No script does either.** Script 17 refers to "the unzip loop", which is not in the project. | **Not fixed** (archive layout unknown). Add the loop as `05b_unzip_spatial.R` or write the exact commands in the README. |
| A7 | 00_run_all | (a) The log was **overwritten** at every stage (`writeLines` to a file name), so `run_all_log.tsv` held only the last stage. (b) Six sentinels were directories created at the *start* of a script, so a crashed stage counted as done. Stage 3's sentinel was never written (it reran every time), and stage 7's sentinel was a folder made by script 08. (c) Stage 25's sentinel `07_figures` already exists after scripts 03 and 11, so **21_figures never ran**. When forced, it overwrote the 8-program figures under the same `*_updated` names. (d) Stages ran with `source()` in one shared session, which is how the `Matrix` defect in script 21 stayed hidden. | Log appends. Sentinels are end-of-script files. Stage 25 removed (40 supersedes it, and 41 sources 40). Stages run in a fresh `Rscript` process by default (`ISOLATE <- TRUE`), with per-stage `.out`/`.err` logs. Tested: isolation, failure reporting, skipping and backup exclusion all work. |

## B. Issues that can affect results or interpretation (action needed)

**B1. The RCTD reference may favour donors stored early in the h5 file (script 17). Run `17b_check_reference_donor_balance.R`.**
Script 17 does the following:
1. Samples 1,200 cells per type at random.
2. **Re-sorts them by h5 position.**
3. Keeps the first 800 that pass QC.

If cells are ordered by patient or channel in the h5 file, donors near the end are under-sampled. In a simulation with 20 donors stored in contiguous blocks, donors 17–20 received **0 cells** and donor 16 received 2; the fixed version gave 39–41 cells to every donor. This affects every reference type with more than 800 eligible cells, and through it all RCTD weights, criterion A adjustment and criterion B quartiles.
- If 17b shows a clear negative correlation between a donor's selection rate and its h5 position, run `FIX_SELECTION_ORDER <- TRUE` in script 17. This writes a separate `_randomsel.rds`. Then repeat 18 → 18b → 19 → 22 → 23 → 36–38 with that reference **as a sensitivity analysis**, and report it.
- Do not quietly replace the primary reference after the results are known. Add a protocol amendment.

**B2. Criterion B is close to circular wherever it uses composition quartiles.**
The expected-high and comparison spots are defined by the RCTD proportion of the same cell type whose genes make up the target. This covers:
- endothelial programs in every cohort;
- the myeloid fallback rule;
- every compartment in GSE280318 and GSE226997 (script 22);
- the lymphoid T/NK rule.

Spots with more endothelial cells show higher endothelial-gene expression whatever the program activity. RCTD also uses many of the same genes to estimate the proportions. B is therefore mostly a check of composition and anatomy, not of program activity.

In addition, script 23 pools two different kinds of B for epithelial and stromal programs: pathologist labels in Valdeolivas and composition quartiles in the extra cohorts. **Action:** state in Methods and Discussion which rule each patient's B used, and describe composition-quartile B as "anatomical plausibility". Consider reporting B separately for each rule.

**B3. Compartment definitions differ between discovery and validation.**

| Compartment | Discovery | Validation cohorts |
|---|---|---|
| Epithelial | stem/TA-like cE01–cE03 only | CMS-labelled cells (SMC/KUL); TC1–4 (Uhlitz); all epithelial clusters (GSE200997) |
| Myeloid | Mono + Macro | Lee myeloids minus cDC/Unknown; **Uhlitz "Monocyte 1/2" only**; marker clusters (GSE200997) |

This shift in the cell population can explain some of the epithelial offsets.

**Action:** check whether Uhlitz has macrophage labels that were left out:
`fread("03_original_data/GSE166555/GSE166555_meta_data.tsv.gz")[grepl("mac|mono|myel|DC", minor_cell_type, ignore.case=TRUE), .N, by=minor_cell_type]`
Then disclose the definitions in a supplementary table.

**B4. Supplementary Table S7b in script 27 came from the wrong file.**
`list.files("06_results/panel", pattern = "stab")[1]` picks `selection_stability.csv` (v1) or `selection_stability_lymphoid.csv`, never `_v2`. Script 41, which supersedes 27, reads the files explicitly and is correct. **Action:** make sure the submitted `Supplementary_Tables.xlsx` was produced by 41. Mark 27 as superseded or remove it. Script 25 had the same bug and is fixed.

**B5. RCTD labels in the extra cohorts are assumed, not verified (script 22).**
If RCTD returns weights without dimnames, script 22 assigns cell-type names in alphabetical order and spot names in QC order. The marker check that follows is written to a file but does not stop the run, unlike the `labels_verified` guard used in 18b and 19. **Action:** confirm `own_is_max` is TRUE for every row of `06_results/spatial_support_extra/RCTD_marker_check_extra.csv`.

**B6. GSE166555 barcodes could be duplicated (scripts 13, 34).**
Per-sample matrices are column-bound and then matched to metadata with `match()`. If a barcode occurs in two samples, the second sample's cells get the first sample's metadata. A `stopifnot` guard has been added. **Action:** re-run 13 and 34. Results should be identical. If the guard stops the run, the current GSE166555 results are wrong.

## C. Items to disclose in Methods / Limitations

1. The expression filter (CPM ≥ 1 in ≥ 25% of donors) and the evaluation/predictor split were made on **all discovery donors, outside the CV folds**. Both are unsupervised and do not use the target, so any optimism is mild, but the protocol says all filtering happens inside training data.
2. **No doublet filtering** was applied in any cohort (scDblFinder was installed but never called). Either remove the doublet claim or apply it.
3. **QC differs by data type:** fixed thresholds for the RCTD reference (≥500 UMI, ≥200 genes, <25% mt) and for spots (≥500 UMI, ≥250 genes, mt 3-MAD), and 3-MAD per channel or sample for all single-cell pseudobulk.
4. **Gene identifiers:** the frozen object says "match by Ensembl ID where available", but every validation and spatial cohort is matched **by symbol**. Duplicated symbols are summed in GSE132465 but kept as `make.unique` copies in GSE200997 and the reference.
5. **Meta-analysis Spearman SE (scripts 16, 35):** `1/sqrt(n−3)` is the Pearson value; for Spearman, about `sqrt(1.06/(n−3))` (Fieller) is more accurate, so pooled-ρ CIs are about 3% too narrow. The primary MAE-ratio pooling is unaffected.
6. **TMM sensitivity scenario (script 24):** normalisation factors use all donors (outside folds). This applies to a sensitivity analysis only.
7. The **RCTD "Myeloid" type** includes DCs and granulocytes, and malignant epithelium covers cE01–03 only. Both are known items.

## D. Minor and figure items

- **Figure 2 (forest) warning, now explained.** `annotate("rect", xmin = 0, …)` on a log-10 axis gives log10(0) = −Inf, which produces the warning. It is harmless: the shading still reaches the left edge, and no interval is lost.
- **Figure S2/3A error bars** are the per-repeat min/max of the panel MAE divided by the *mean* baseline MAE, not the range of per-repeat ratios. Reword the caption, or save per-repeat ratios in script 10.
- **Hard-coded numbers to check against the final results:**
  - Figure 1: Valdeolivas row "n = 6" and caption "6 patients, 12 sections", while the spatial analysis has 13 patients from 3 cohorts.
  - Figure S3A title: "12 sections". Script 18b averages over *every* section in `rctd/`, including the excluded patient A798015.
  - Figure 6 caption: "endothelial: single eligible cohort GSE132465".
  - `COH` labels in figure 5: n = 6/4/3.
- **`28_CD8_feasibility.R` crashes:** `readRDS()` is called on `gene_partition.sha`. The script is exploratory and not in the pipeline. Remove it or fix it, and check that the protocol doesn't quote its "blocked genes" counts.
- **Latent edge cases (no effect on current outputs):**
  - Script 11: within-cohort `scale()` has no zero-variance guard; 13 and 15 have one.
  - Script 03: `mad_flag()` has no zero-MAD guard; 08 has one.
  - Scripts 13 and 15 label a non-finite ratio differently ("Not evaluable" vs "Not supported").
  - Script 39 has no fallback when the pooled ratio is NA; 20 has one.
  - Script 10: the λ fallback can return more than 4 genes or NA. Confirm every frozen panel has 1–4 genes.
- **Scripts that rewrite other scripts:** 27 and 41 edit 21 and 40 in place (axis step). This is now a no-op; delete step 1.
- **Version pinning:** script 05 downloads Zenodo record 7739700 as `versions/latest`. Record the version ID that was actually used.
- **Archive cleanup:** exclude `_before_*` backups and script 12, or keep them in an `archive/` subfolder.

## E. Verified correct

- **Raw counts.** Integer counts are asserted at load for every cohort. Pseudobulk sums raw UMIs, and log-CPM is computed only after aggregation. No normalised matrix is treated as counts.
- **h5 reading and pseudobulk.** h5 block reading is correct (0-based indices, `indptr` windows), as is indicator-matrix aggregation. Donor-level pooling of repeat specimens is correct (C130/C171 TA+TB; Uhlitz t1/t2).
- **QC.** 3-MAD QC is applied per channel or sample. Script 29 reproduces script 08's per-channel removals exactly.
- **Gene partition.** Evaluation and predictor genes are disjoint, enforced by `stopifnot` in both partitions. The lymphoid partition respects the existing evaluation genes.
- **Panel development.** Target z-parameters, predictor scaling, glmnet tuning, gene selection, the equal-weight top-4 baseline and conventional-marker calibration all happen inside training folds. Folds are stratified by chemistry. The discovery-mean baseline is exactly 0 on the frozen scale, so `mean(|y|)` is the right MAE denominator.
- **Frozen validation.**
  - The SHA-256 is checked before every downstream use.
  - No refitting or reselection.
  - A missing panel gene makes the program not evaluable.
  - Evaluation-gene coverage ≥ 80% and ≥ 10 donors are required.
  - Bootstrap CIs are percentile; the pooled SMC + KUL analysis is stratified by cohort.
- **Meta-analysis.** DerSimonian–Laird τ² and I² are correct, the log-ratio SE is derived from the bootstrap CI, and cohorts with <10 donors are excluded from the main analysis.
- **Spatial.**
  - Composition adjustment drops Mast to avoid exact collinearity of sum-to-one weights.
  - An 8×8 block bootstrap handles spatial dependence.
  - Patient-level SE averages bootstrap replicates across sections.
  - Pooling is random-effects with a ≥ 2/3-of-patients rule.
  - LOPO now uses the 13-patient file.
  - The excluded non-neoplastic patient (A798015) is applied consistently.
- **Package loading.** No other library-attachment dependencies were found (the script 21 `Matrix` defect is already fixed).

## Corrected scripts

| File | Change |
|---|---|
| `00_run_all.R` | sentinels, appended log, isolated stages, stage 25 removed, manual-step notes |
| `00_setup_and_versions.R` | complete package list, Bioconductor spacexr check, package source recorded |
| `04_download_GSE178341_GSE144735.R` | rebuilt; verifies the originals are unchanged |
| `12b_KUL_from_saved.R` | builds KUL counts if missing; duplicate-barcode guard |
| `13_validation_GSE166555.R`, `34_lymphoid_extra_cohorts.R` | duplicate-barcode guard; extract to `04_processed_data` |
| `17_spatial_reference_RCTD.R` | Bioconductor spacexr required; opt-in fixed selection to a separate file; sessionInfo |
| `17b_check_reference_donor_balance.R` | new diagnostic for B1 |
| `23_spatial_combined.R` | writes `spatial_support_summary_main.csv` |
| `25_audit_manuscript_numbers.R` | no setwd; correct column and stability file; `_main.csv` consistency check |

## Run order on your machine

1. Copy the corrected scripts into `05_analysis_scripts/`. Keep the old ones in `archive/`.
2. Run `23_spatial_combined.R`, then `25_audit_manuscript_numbers.R`, and read section 10. If it says **DIFFERENT**, re-run 20 → 39 → 41 and re-check every evidence category and spatial number in the manuscript.
3. Run `17b_check_reference_donor_balance.R` (B1).
4. Re-run 13 and 34 (B6). The outputs should be identical (check with `26_compare_outputs.R`).
5. Check `RCTD_marker_check_extra.csv` (B5) and the Uhlitz myeloid labels (B3).
6. Add the spatial extraction step (A6). Then a full clean-room run with `RUN_FROM <- 11` (keeping the supplied frozen panel) is meaningful.

---

## Follow-up (19 Oct 2026, new PC, D: drive)

- The new PC now has a complete package library (25 packages). Its versions are nearly identical to the original runs (Seurat 5.6.0 vs 5.5.1, edgeR 4.10.6 vs 4.10.5). The original record is kept as `06_results/software_versions_original_runs.csv`.
- **A2 resolved:** the original `spatial_support_summary_main.csv` matches the main rows written by the corrected script 23 exactly, for all 7 programs (n, effect, CI, A/B counts, label). The evidence table was built from correct values, so 20 → 39 → 41 do not need re-running.
- Audit (corrected script 25) confirms: frozen panel SHA-256, discovery donors 60/27/60/40, 15,904 Valdeolivas spots, 26-gene chip table, evidence categories (END_TNFA_NFKB highest; EPI_HYPOXIA single-cell only; END_ANGIO unsupported; 4 exploratory).
- **New check for the manuscript:** the sensitivity analysis has **78 fits, 75 below 1**, with 3 exceptions: v2-only angiogenesis 1.115, v2-only TGF-β 1.204, partition seed 3 TGF-β 1.078. The audit header expected "82/85". Any manuscript text saying 82/85 must be corrected to 75/78.
- Still open: 17b donor-balance diagnostic (B1); re-run 13 and 34 with the barcode guard (B6); `own_is_max` check in `RCTD_marker_check_extra.csv` (B5; minimum own-marker ρ is −0.08 in GSM7089856_P2 and 0.06 in GSM8594568_P2CRC); Uhlitz macrophage labels (B3); spatial unzip step (A6).

### Follow-up 2 (20 Oct 2026)
- **B6 closed:** script 34 ran with the duplicate-barcode guard and reproduced protocol v1.9 exactly (GSE166555 12 donors, 0.214 [0.122-0.298] Supported; KUL 2 eligible donors, not evaluable). GSE166555 has no duplicated barcodes, so script 13 is unaffected.
- **B1 confirmed:** 17b shows h5-order bias for every reference type (selection rate vs h5 position rho -0.47 to -0.73; TVD 0.21-0.36). Donors represented vs eligible: B_Plasma 47/62, Endothelial 47/60, Fibroblast 44/58, Malignant_epithelial 55/62, Mast 35/58, Myeloid 47/62, Normal_epithelial 34/36, Pericyte_SMC 47/59, T_NK 43/62. Sensitivity analysis script `42_spatial_sensitivity_random_reference.R` added: it rebuilds the reference with random selection and re-runs 18, 18b, 19, 22, 23, 36-38 into `*_randomsel` and `06_results/sens_randomref/` (redirection tested on the real script texts; no published path is written). The published reference remains primary unless a protocol amendment says otherwise.
- **B5:** 7 flags (`own_is_max == FALSE`): Normal_epithelial in 4 tumour sections (proportion 0.5-10%), Pericyte_SMC in 2 Visium HD sections (own rho 0.75-0.83, likely beaten by the Fibroblast column through ACTA2), Myeloid in GSE226997 P4 (rho 0.28). No evidence of permuted labels. Disclose in the Figure S3 legend, especially the weak myeloid deconvolution in GSE226997 P4.
- **B3:** GSE166555 has only "Monocyte 1/2" (2,075 + 1,216 cells) and no macrophage label, so its myeloid compartment is monocytes only. Disclose.

### Follow-up 3 (20 Oct 2026): spatial sensitivity analysis with the random-selection reference (script 42)
- Corrected reference donors per type: B_Plasma 52, Endothelial 58, Fibroblast 53, Malignant_epithelial 62, Mast 55, Myeloid 60, Normal_epithelial 34, Pericyte_SMC 56, T_NK 61 (published: 35-55). Script 18b verified the labels (each marker set matches its own column: TRUE).
- All spatial cohorts (main), published vs corrected reference: END_ANGIO 0.015 vs 0.016; END_TNFA_NFKB 0.203 vs 0.203; EPI_APOPTOSIS 0.156 vs 0.157; EPI_HYPOXIA 0.055 vs 0.055; LYM_IL2_STAT5 0.146 vs 0.148; MYE_IFNG 0.241 vs 0.243; MYE_INFLAM 0.147 vs 0.141; STR_TGFB 0.110 vs 0.111. CIs nearly identical.
- **One label change:** EPI_APOPTOSIS moved from Spatial support to Partial. Criterion B fell from 9/13 to 8/13 (threshold 9); the flip is in an extra-cohort patient scored by the composition-quartile rule. Criterion A unchanged (13/13). Valdeolivas pathologist-label B unchanged (4/6).
- **No evidence category changes** (EPI_APOPTOSIS is Exploratory on its single-cell result). Recommendation: keep the published reference as primary (pre-specified); report this as a sensitivity analysis, disclose the apoptosis label change, and record it in protocol v2.1.
- Implementation note: in the single-session runner, SingleCellExperiment (attached by script 18) masked edgeR::cpm; the runner now qualifies bare cpm() calls. Not an issue in 00_run_all.R, where each stage runs in its own process.

### Follow-up 4 (20 Oct 2026): publication figures (script 43)
- `43_publication_figures.R` redraws all figures in one visual system: fixed compartment colours (passes adjacent-pair colour-vision checks, always with direct labels), a sequential green ramp for evidence, 7 pt text at 178 mm, and every count read from the result files. Outputs go to `07_figures/publication/*_aem.{pdf,tiff,png}`.
- Final set and proposed numbering: **Fig 1** study design + cohort coverage (spatial rows now 6/3/4 = 13 patients); **Fig 2** validation forest (frozen and within-cohort, numeric columns, pooled diamonds); **Fig 3** spatial support (example maps with Spearman rho; criteria A and B per patient); **Fig 4** evidence scorecard; **Fig S1** discovery donor variation; **Fig S2** panel baselines + selection stability; **Fig S5** spatial reference sensitivity (new). S3/S4 still come from script 41 (RCTD marker check; spatial by cohort group).
- Figures checked against real data on the author's PC: S1, S2, 2 and 3. Fixes applied: Fig 2 column titles; S2 dashed line behind bars and fewer ticks. Open choice: the Fig 3A example shows END_TNFA_NFKB with rho 0.19 in section SN048_A416371_Rep1 (IFN-g rho 0.57); consider another program or section for illustration.
- Manuscript call-outs and legends must be updated to this numbering (legends for Fig 1 [now includes the workflow panel], Fig 4 [scorecard] and new Fig S5).
