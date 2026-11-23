###############################################################################
# 64_collect_revision_tables.R - STEP 14: assemble outputs for the manuscript
#
# Writes 08_tables/revision_2026-11/Supplementary_Tables_revision.xlsx with one
# sheet per new analysis (only for analyses that have been run), plus S1d and
# S14, and a plain-text file listing the numbers needed for each [PENDING]
# placeholder in Biomarker_Manuscript_V2_tracked.docx.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision()
suppressPackageStartupMessages(library(openxlsx))

rd <- function(n) { f <- file.path(OUT_DIR, paste0(n, ".csv")); if (file.exists(f)) fread(f) else NULL }
sheets <- list(
  S1d_patients_by_program = {
    s1c <- read_supp("S1c_validation_patients")
    setnames(s1c, c("Program code", "Cohort", "Patients (n)"), c("program", "cohort", "n"))
    w <- dcast(s1c, program ~ cohort, value.var = "n"); w[, total := rowSums(.SD, na.rm = TRUE), .SDcols = -1]; w },
  S7c_alternative_panels    = rd("61_alternative_panels"),
  S7d_mmr_site_strata       = rd("62_mmr_site_strata"),
  S8c_calibration_transport = rd("56_calibration_transport"),
  S8d_leave_one_donor_out   = rd("54_leave_one_donor_out"),
  S8e_meta_sensitivity      = rd("55_meta_analysis_sensitivity"),
  S8f_multiplicity          = rd("57_multiplicity_fdr"),
  S10e_spatial_table_sens   = rd("58_spatial_table_sensitivity"),
  S10f_spatial_rerun        = rd("59_spatial_rerun_pooled"),
  S10g_spatial_negcontrol   = rd("59_spatial_negative_control"),
  S14_amendment_log = data.table(
    version = c("v1.0", "v1.3", "v1.4", "v1.5", "v1.8", "v1.9", "v2.1"),
    date = c("2026-10-09", "2026-10-13", "2026-10-13", "2026-10-13", "2026-10-15", "2026-10-15", "2026-10-20"),
    change = c("DerSimonian-Laird pooling of cohorts with >=10 donors", "Evidence-category assignment rule; chip assay rule",
               "A798015 excluded as non-neoplastic (existing rule)", "GSE280318 and GSE226997 added; composition-quartile criterion B",
               "CD8 T-cell compartment added (separate frozen panel)", "CD8 validation extended to all cohorts; CD8 spatial parameters",
               "Post hoc spatial sensitivity analysis with random-selection RCTD reference"),
    data_already_seen = c("GSE132465 single-cohort results", "All single-cell and Valdeolivas spatial results", "Valdeolivas spatial run",
                          "Valdeolivas spatial results", "Four-compartment validation results", "GSE132465 CD8 result", "All results"),
    disposition = c("Post-result; reported with timing", "Post-result; grading descriptive", "Eligibility rule unchanged",
                    "Exploratory; Valdeolivas-only reported", "Soft blinding lost; disclosed", "All cohorts reported; pooled is headline",
                    "Sensitivity only; no category change"))
)
f53 <- rd("53_null_random_panels_summary"); if (!is.null(f53)) sheets$S16_null_random_panels <- f53
sheets <- sheets[!vapply(sheets, is.null, logical(1))]

wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#D9E2F3", border = "bottom")
for (nm in names(sheets)) {
  addWorksheet(wb, substr(nm, 1, 31))
  writeData(wb, substr(nm, 1, 31), sheets[[nm]], headerStyle = hs)
  setColWidths(wb, substr(nm, 1, 31), cols = seq_len(ncol(sheets[[nm]])), widths = "auto")
  freezePane(wb, substr(nm, 1, 31), firstRow = TRUE)
}
xf <- file.path(TAB_DIR, "Supplementary_Tables_revision.xlsx")
saveWorkbook(wb, xf, overwrite = TRUE); msg("wrote", xf)

## ---- numbers for the [PENDING] placeholders ------------------------------------
L <- c("NUMBERS FOR [PENDING] PLACEHOLDERS (generated", format(Sys.time()), ")", "")
x <- sheets$S8d_leave_one_donor_out
if (!is.null(x)) { r <- x[program == "END_TNFA_NFKB" & cohort == "GSE132465"]
  if (nrow(r)) L <- c(L, sprintf("Point 4: TNF/NF-kB LODO ratios %.2f-%.2f; jackknife 95%% CI %.2f-%.2f; most influential %s (ratio without: %.2f); all <1: %s",
                                 r$loo_min, r$loo_max, r$jk_lo, r$jk_hi, r$most_influential, r$ratio_without_it, r$all_loo_below_1)) }
x <- sheets$S8e_meta_sensitivity
if (!is.null(x)) { r <- x[scaling == "frozen" & cohorts == "main" & analysis == "all cohorts" & method == "REML + HK"]
  for (i in seq_len(nrow(r))) L <- c(L, sprintf("Point 6: %s REML+HK pooled %.2f (%.2f-%.2f), k=%d, supported=%s",
                                                r$program[i], r$est[i], r$lo[i], r$hi[i], r$k[i], r$supported[i])) }
x <- sheets$S8c_calibration_transport
if (!is.null(x)) { r <- x[n >= MIN_DONORS]
  L <- c(L, sprintf("Point 8: calibration slopes %.2f to %.2f; offsets %.2f to %.2f; median frozen ratio %.2f vs LOO-linear-recalibrated %.2f",
                    min(r$cal_slope), max(r$cal_slope), min(r$offset), max(r$offset), median(r$ratio_frozen), median(r$ratio_linear_loo))) }
x <- sheets$S8f_multiplicity
if (!is.null(x)) { r <- x[family == "per cohort"]
  L <- c(L, sprintf("Point 15: %d of %d per-cohort tests supported unadjusted; %d at BH-FDR 5%%",
                    sum(r$supported_unadjusted), nrow(r), sum(r$supported_FDR05))) }
x <- sheets$S16_null_random_panels
if (!is.null(x)) for (i in seq_len(nrow(x))) if (grepl("validation|mean log", x$setting[i]))
  L <- c(L, sprintf("Point 7: %s | %s | %s: observed %.3f, null median %.3f (95%% %.3f-%.3f), p=%.3f",
                    x$program[i], x$null[i], x$setting[i], x$observed[i], x$null_median[i], x$null_q025[i], x$null_q975[i], x$p_empirical[i]))
x <- sheets$S10e_spatial_table_sens
if (!is.null(x)) { r <- x[(label_changed | category_changed)]
  L <- c(L, "Points 11/14: scenarios changing a spatial label or evidence category:",
         if (nrow(r)) sprintf("  %s | %s: %s; %s", r$scenario, r$program, r$spatial_label, r$evidence_category) else "  none") }
x <- sheets$S10f_spatial_rerun
if (!is.null(x)) L <- c(L, "Points 9/10: pooled criterion A by block size / model / weights in S10f (compare block 8 with published S10c).")
x <- sheets$S7c_alternative_panels
if (!is.null(x)) L <- c(L, "Point 13: alternative panels in S7c (compare next-best and swap rows with 'frozen').")
x <- sheets$S7d_mmr_site_strata
if (!is.null(x)) L <- c(L, "Point 26: subgroup CV MAE ratios and ratio-of-ratios in S7d.")
writeLines(L, file.path(TAB_DIR, "numbers_for_manuscript.txt"))
cat(L, sep = "\n")
