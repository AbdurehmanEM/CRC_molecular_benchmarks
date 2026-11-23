###############################################################################
# 50_revision_config.R  -  THE ONLY FILE YOU NEED TO EDIT
#
# Revision analyses requested in the pre-submission review (11 Nov 2026).
# Every revision script (51-64) sources this file first.
#
# Before running anything:
#   1. Extract "CRC_Biomarker_Project 2.rar" (R cannot read inside a .rar).
#   2. Set PROJECT_ROOT below to the extracted folder.
#   3. Fill in PB_FILES (pseudobulk count files written by your pipeline).
#   4. Run 51_inspect_inputs.R and read its report.
#   5. Run 52_reproduce_published.R. It is a GATE: if your published MAE
#      ratios are not reproduced, stop and fix the loader or options here
#      before trusting any new result.
###############################################################################

## ---- 1. Project location ----------------------------------------------------
# Windows path with forward slashes. The archive must be extracted first.
PROJECT_ROOT <- Sys.getenv("CRC_PROJECT_ROOT", unset = "D:/Manuscripts/CRC_Biomarker_Project 2/CRC_Biomarker_Project")

## ---- 2. Output folders (created if missing; nothing published is overwritten)
OUT_DIR  <- file.path(PROJECT_ROOT, "06_results", "revision_2026-11")
TAB_DIR  <- file.path(PROJECT_ROOT, "08_tables",  "revision_2026-11")
FIG_DIR  <- file.path(PROJECT_ROOT, "07_figures", "revision_2026-11")

## ---- 3. Supplementary workbook (source of the frozen panels, partition and
##         published numbers used by the reproduction gate)
SUPP_XLSX <- file.path(PROJECT_ROOT, "10_submission", "Supplementary_Tables_aem.xlsx")
# If your final workbook lives elsewhere (e.g. "Supplementary_Tables.xlsx"), change it here.

## ---- 4. Frozen panel files (only used to confirm the hash has not changed) --
FROZEN_FILES <- c(
  four_compartment = file.path(PROJECT_ROOT, "04_processed_data", "frozen_panel", "frozen_panel_v2.rds"),
  lymphoid         = file.path(PROJECT_ROOT, "04_processed_data", "frozen_panel", "frozen_panel_LYM_v1.rds"))
FROZEN_SHA256 <- c(
  four_compartment = "a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb",
  lymphoid         = "7dc3256d97477a796fcacd1df71021740888f7978a01ea36caf70787da7f513f")

## ---- 5. Pseudobulk inputs (set for YOUR files, 11 Nov 2026) -----------------
# Discovery: raw summed UMI counts, list(counts = dgCMatrix genes x donors, n_cells = ...);
#            log-CPM is computed here and donors with < MIN_CELLS cells are dropped.
# Validation: log-CPM matrices already written by your validation scripts
#            (11/12b/13/15/30), one per compartment, symbols as rownames.
PB_DIR <- file.path(PROJECT_ROOT, "04_processed_data")
PB_FILES <- list(
  GSE178341 = list(type = "counts",
                   Epithelial  = "GSE178341/pseudobulk_Epithelial.rds",
                   Stromal     = "GSE178341/pseudobulk_Stromal.rds",
                   Myeloid     = "GSE178341/pseudobulk_Myeloid.rds",
                   Endothelial = "GSE178341/pseudobulk_Endothelial.rds",
                   Lymphoid    = "GSE178341/pseudobulk_Lymphoid_CD8.rds"),
  GSE132465 = list(type = "logcpm_nested",          # list(logcpm = list(<cmp> = matrix), ncells = list(...))
                   Epithelial  = "GSE132465/validation_pseudobulk_logcpm.rds",
                   Stromal     = "GSE132465/validation_pseudobulk_logcpm.rds",
                   Myeloid     = "GSE132465/validation_pseudobulk_logcpm.rds",
                   Endothelial = "GSE132465/validation_pseudobulk_logcpm.rds",
                   Lymphoid    = "GSE132465/validation_pseudobulk_logcpm_lymphoid.rds"),
  GSE166555 = list(type = "logcpm_list",            # list(<cmp> = matrix)
                   Epithelial  = "GSE166555_validation_pseudobulk_logcpm.rds",
                   Stromal     = "GSE166555_validation_pseudobulk_logcpm.rds",
                   Myeloid     = "GSE166555_validation_pseudobulk_logcpm.rds",
                   Endothelial = "GSE166555_validation_pseudobulk_logcpm.rds",
                   Lymphoid    = NA),   # not saved by script 34 (see 50c_save_lymphoid_validation.R)
  GSE200997 = list(type = "logcpm_list",
                   Epithelial  = "GSE200997/GSE200997_validation_pseudobulk_logcpm.rds",
                   Stromal     = "GSE200997/GSE200997_validation_pseudobulk_logcpm.rds",
                   Myeloid     = "GSE200997/GSE200997_validation_pseudobulk_logcpm.rds",
                   Endothelial = "GSE200997/GSE200997_validation_pseudobulk_logcpm.rds",
                   Lymphoid    = NA),   # T and NK not separable (excluded by design)
  GSE144735 = list(type = "logcpm_list",
                   Epithelial  = "GSE144735/validation_pseudobulk_logcpm.rds",
                   Stromal     = "GSE144735/validation_pseudobulk_logcpm.rds",
                   Myeloid     = "GSE144735/validation_pseudobulk_logcpm.rds",
                   Endothelial = "GSE144735/validation_pseudobulk_logcpm.rds",
                   Lymphoid    = NA)    # 2 donors only (not evaluable)
)
# Element name used for the lymphoid compartment inside the nested GSE132465 file:
LYMPHOID_ELEMENT <- c("Lymphoid", "CD8", "Lymphoid_CD8")
# If 50c_save_lymphoid_validation.R has written extra lymphoid files, they are used automatically:
EXTRA_LYMPHOID_DIR <- file.path(PB_DIR, "revision_lymphoid_validation")

# Discovery log-CPM: computed on eligible donors only (FALSE, protocol default) or on all
# donors before filtering (TRUE). Library-size prior scaling makes these differ slightly;
# if the gate fails with small ratio differences, try TRUE.
DISCOVERY_LOGCPM_ALL_DONORS <- FALSE

DISCOVERY   <- "GSE178341"
VALIDATION  <- c("GSE132465", "GSE166555", "GSE200997")   # main (>=10 donors rule applies)
DESCRIPTIVE <- c("GSE144735")                            # shown, never pooled

# Discovery feature table (row_id, symbol) for every discovery gene; written by
# 50b_make_feature_map.R from the GSE178341 h5 file. Needed by scripts 53 and 61.
DISCOVERY_FEATURES <- file.path(PROJECT_ROOT, "04_processed_data", "GSE178341", "GSE178341_feature_map.csv")
GSE178341_H5_DIR   <- file.path(PROJECT_ROOT, "03_original_data", "GSE178341")

## ---- 6. Fixed analysis rules (protocol v0.6-v0.9; do NOT change) -----------
MIN_CELLS      <- 20      # cells per donor x compartment profile
MIN_DONORS     <- 10      # donors per cohort for evaluation / pooling
MIN_EVAL_COVER <- 0.80    # evaluation-gene coverage
PRIOR_COUNT    <- 2       # log2-CPM prior count (library size only, no TMM)
N_BOOT         <- 2000    # donor bootstrap resamples
SEED           <- 20261111  # seed for all NEW revision analyses

## ---- 7. Reproduction-gate tolerances ----------------------------------------
GATE_TOL_RATIO <- 0.005   # |recomputed - published| MAE ratio, primary (frozen) analysis
GATE_TOL_RATIO_SECONDARY <- 0.02   # same, within-cohort (secondary) analysis
GATE_TOL_CI_REL <- 0.10   # CI limits: tolerance = 10% of published CI width + 0.01 (bootstrap seeds differ)
GATE_STRICT    <- as.logical(Sys.getenv("CRC_GATE_STRICT", unset = "TRUE"))  # TRUE: scripts 53+ refuse to run if the gate failed

## ---- 8. Revision-analysis settings --------------------------------------------
N_NULL_FIXED   <- 1000    # script 53, random fixed panels
N_NULL_SELECT  <- 200     # script 53, selection from random matched pools
N_NULL_TARGET  <- 1000    # script 53, negative-control targets
CV_FOLDS       <- 5
CV_REPEATS     <- 10

## ---- 9. Spatial inputs (scripts 59 and 60 only) -----------------------------
# A per-section .rds file with list(counts = genes x spots (symbols),
#   coords = data.frame(spot, array_row, array_col), patient, cohort,
#   weights = spots x cell-type RCTD proportions (published reference)).
# Write these once from your script 19/22 objects (see the template in
# 59_spatial_rerun.R, section 0). One file per section, in this folder:
SPATIAL_SECTION_DIR <- file.path(PROJECT_ROOT, "04_processed_data", "revision_spatial_sections")
SPATIAL_EXCLUDE_PATIENTS <- c("A798015")      # non-neoplastic (protocol v1.4)
SPATIAL_BLOCK_SIZES <- c(4, 8, 12, 16)        # 8 = published
SPATIAL_N_BOOT      <- 1000
SPATIAL_N_NEGCTRL   <- 100
# Options that must match script 19 (the gate in 59 checks them against S10a):
SPATIAL_LOGNORM     <- "log2cpm"   # "log2cpm" (prior 2) or "log1p_cp10k"
SPATIAL_DROP_TYPE   <- "Mast"      # omitted proportion (collinearity)

# Independent RCTD reference (script 60): GSE132465 single cells
GSE132465_COUNTS <- file.path(PROJECT_ROOT, "03_original_data", "GSE132465",
                              "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz")
GSE132465_ANNOT  <- file.path(PROJECT_ROOT, "03_original_data", "GSE132465",
                              "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz")
REF_CELLS_PER_TYPE <- 800

## ---- 10. Optional validation-cohort MMR status (script 62) ------------------
# CSV with columns cohort, donor, mmr ("MMRd"/"MMRp"); NA to skip validation strata.
VALIDATION_MMR_FILE <- NA

## ---- 11. Clean-room re-run (script 63) ---------------------------------------
CLEAN_ROOM_DIR <- "D:/CRC_clean_room"
CLEAN_ROOM_RUN_FROM <- 11
