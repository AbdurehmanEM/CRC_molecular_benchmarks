# Revision analyses - R code, step by step

Code for the 11 analyses marked **[PENDING point n]** in `Biomarker_Manuscript_V2_tracked.docx`.
Written 11 November 2026. Tested end to end on a simulated copy of the project with the real
`Supplementary_Tables.xlsx` (R 4.3.3; glmnet 4.1-8, metafor 4.4-0, edgeR 4.0.16, data.table 1.14.10).
Script 60 (spacexr) and the clean-room part of 63 could not be tested here (see "Untested parts").

## 0. Setup (once)

1. **Extract** `F:\CRC_Biomarker_Project 2.rar`. R cannot read files inside a `.rar`.
   The config assumes `D:/Manuscripts/CRC_Biomarker_Project 2/CRC_Biomarker_Project`; change `PROJECT_ROOT` if you extract elsewhere.
2. Copy this folder to `05_analysis_scripts/revision/`.
3. Install missing packages:
   ```r
   install.packages(c("data.table", "Matrix", "readxl", "openxlsx", "glmnet", "metafor", "digest"))
   BiocManager::install(c("edgeR", "spacexr", "SpatialExperiment"))   # spacexr only for script 60
   ```
4. Open `50_revision_config.R`. **It is the only file you edit.** Set:
   - `SUPP_XLSX`: your final supplementary workbook;
   - `PB_FILES`: the pseudobulk **raw count** files your pipeline wrote, one per cohort x compartment
     (genes x donors, all genes; discovery rownames = matrix row IDs as in S5/S6, validation = symbols).
     A file can be a matrix or a `list(counts = , n_cells = )`;
   - `DISCOVERY_FEATURES`: a CSV mapping every discovery row ID to its symbol (from the GSE178341 h5 features).

## Run order

| Step | Script | Reviewer point | Needs | Runtime |
|---|---|---|---|---|
| 1 | `51_inspect_inputs.R` | all | config | seconds |
| 1b | `50c_save_lymphoid_validation.R` (optional) | saves GSE166555 CD8 log-CPM so the gate can check it | script 34 | as script 34 |
| 2 | `52_reproduce_published.R` | **gate** | pseudobulk | ~1 min |
| 3 | `55_meta_analysis_sensitivity.R` | 6 (REML/Hartung-Knapp, leave-one-cohort-out) | workbook only | seconds |
| 4 | `57_multiplicity_fdr.R` | 15 (BH-FDR, Holm) | workbook (+ step 2 for bootstrap p) | seconds |
| 5 | `58_spatial_table_sensitivity.R` | 11, 14 (thresholds, without GSE226997, criterion-B rule) | workbook only | seconds |
| 6 | `54_leave_one_donor_out.R` | 4 | step 2 | seconds |
| 7 | `56_calibration_transport.R` | 8 | step 2 | seconds |
| 8a | `50b_make_feature_map.R` | discovery row ID -> symbol map | GSE178341 h5 + rhdf5/hdf5r | seconds |
| 8 | `53_null_random_panels.R` | 7 (most important) | pseudobulk + feature map | ~0.5-2 h |
| 9 | `61_alternative_panels.R` | 13 | pseudobulk | minutes |
| 10 | `62_mmr_site_strata.R` | 26 | pseudobulk | minutes |
| 11 | `60_independent_rctd_reference.R` | 10 | spacexr, GSE132465 raw files, exported sections | hours |
| 12 | `59_spatial_rerun.R` | 9, 10 | exported sections (section 0 of the script) | ~1 h |
| 13 | `63_reproducibility_clean_room.R` | 27 | digest; clean room optional (~12 h) | minutes / hours |
| 14 | `64_collect_revision_tables.R` | all | whatever has run | seconds |

Run everything with `source("00_run_revision.R")` from the `revision` folder (each script runs in a fresh R
process; logs go to `logs/`). Or run the scripts one at a time in the order above.

**The gate (step 2) matters.** It recomputes every published validation result (S8 and S8b) from your
pseudobulk files with the same helper functions the new analyses use. If it fails, scripts 53, 54, 56, 61
and 62 refuse to run, because their numbers would not be comparable to the published ones. Read
`52_gate_comparison.csv`: `donors_match = FALSE` means the wrong file or donors not filtered to >=20 cells;
large ratio differences usually mean discovery rownames are symbols instead of row IDs. If ratios match
and only the CIs differ (different bootstrap seed), widen `GATE_TOL_CI` and say so in the Methods.
Spatial script 59 has its own gate against S10a (patient coefficients within 0.01).

## Outputs

- `06_results/revision_2026-11/*.csv`: one or more per script, prefixed with the script number.
- `07_figures/revision_2026-11/56_calibration_plots.pdf`.
- `08_tables/revision_2026-11/Supplementary_Tables_revision.xlsx`: new sheets S1d, S7c, S7d, S8c-S8f,
  S10e-S10g, S14 (amendment log) and S16 (random-panel nulls).
- `08_tables/revision_2026-11/numbers_for_manuscript.txt`: the values for each [PENDING] placeholder.

Nothing published is overwritten: all outputs go to new `revision_2026-11` folders, the frozen
panels are only read (their hashes are checked in step 1), and the seed for new analyses is 20261111.

## What each analysis does (for the Methods)

- **53 (point 7).** Three nulls on the frozen discovery scale: (A) fixed random panels of the same size,
  each gene matched to a frozen gene on discovery mean log-CPM decile x SD tercile, from non-program
  genes, fitted by OLS on discovery donors; (B) the full relaxed elastic net run on a random matched pool
  the size of the program's predictor pool, the fair null for the selection procedure; (C) negative-control
  targets: Spearman correlation of the frozen prediction with random gene-set scores matched to the
  evaluation genes. Empirical one-sided p = (1 + #null as good or better) / (1 + B).
- **54 (point 4).** Leave-one-donor-out MAE ratios, most influential donor, jackknife CI.
- **55 (point 6).** DL (reproduces S8b), REML, REML + Hartung-Knapp, Paule-Mandel + HK, HK ad hoc;
  95% prediction intervals; leave-one-cohort-out.
- **56 (point 8).** Offset, calibration intercept/slope, and MAE ratios after leave-one-out offset and
  linear recalibration against a leave-one-out cohort-mean baseline.
- **57 (point 15).** One-sided p for MAE ratio < 1 per cohort and pooled; BH and Holm adjustment.
- **58 (points 11, 14; code review B2).** Reproduces S10c from S10a/S10b, then varies the patient
  threshold (1/2, 2/3, 3/4, all), uses a stricter criterion A (patient lower 95% limit > 0), removes
  GSE226997 or each spatial cohort, and splits criterion B by region rule.
- **59 (points 9, 10).** Spot-level criterion A re-implemented: block sizes 4/8/12/16, negative-control
  targets, a residualised model (+ log UMI + matched negative-control score), and alternative weights.
- **60 (point 10).** RCTD reference from GSE132465 (nine types, <=800 cells per type, round-robin over
  donors), RCTD re-run per section, weights added for script 59.
- **61 (point 13).** Next-best and single-swap panels from S7b selection frequencies (post hoc, exploratory).
- **62 (point 26).** Published CV repeated with selection inside folds; out-of-fold MAE ratios by MMR
  status and tumour site with bootstrap CIs and ratio-of-ratios.
- **63 (point 27).** SHA-256 input manifest, package versions, optional clean-room re-run from stage 11
  with CSV-by-CSV comparison.

## Untested parts - check these first

- **Script 60**: the `createRctd()` / `runRctd()` arguments and the `assay(out, "weights")` extraction
  follow the Bioconductor spacexr API but could not be run here. Copy the exact calls from your
  script 18 (which work) if they differ. Check the printed cell-type mapping table before using the weights.
- **Script 59, section 0**: you must export one `.rds` per section from your script 19/22 objects
  (template in the script). Then check the spatial gate; if it fails, adjust `SPATIAL_LOGNORM` (and,
  if needed, the within-section scoring in `section_scores()`) to match script 19.
- **Script 63, clean room**: assumes `00_run_all.R` contains a line `RUN_FROM <- ...` and that your
  pipeline scripts use paths relative to the project root.
