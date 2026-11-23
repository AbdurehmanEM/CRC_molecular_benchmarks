# =============================================================================
# 00_run_all.R  -  master script for the CRC biomarker prioritization pipeline
#
# Usage (from the project root):
#   source("05_analysis_scripts/00_run_all.R")          # run everything not yet done
#   RUN_FROM <- 9; source(".../00_run_all.R")           # resume from stage 9
#   RUN_ONLY <- c(19, 22, 23); source(".../00_run_all.R")
#   SKIP_DOWNLOADS <- TRUE; source(".../00_run_all.R")  # data already present
#
# Every stage is skipped automatically if its sentinel output already exists,
# so the script is safe to re-source after an interruption.
# =============================================================================

PROJECT_ROOT <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
setwd(PROJECT_ROOT)
stopifnot(dir.exists("05_analysis_scripts"))
SCRIPTS <- "05_analysis_scripts"
LOGDIR  <- file.path("06_results", "run_logs"); dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)
if (!exists("SKIP_DOWNLOADS")) SKIP_DOWNLOADS <- FALSE
if (!exists("RUN_FROM"))       RUN_FROM <- 0
if (!exists("RUN_ONLY"))       RUN_ONLY <- NULL
if (!exists("FORCE"))          FORCE <- FALSE
# ISOLATE = TRUE (default) runs every stage in a fresh `Rscript` process, so a script
# that forgets library() fails here instead of silently borrowing packages or objects
# from an earlier stage (the 21_figures/Matrix defect). Per-stage output goes to
# 06_results/run_logs/stage_<id>.out / .err. Set ISOLATE <- FALSE to source in-session.
if (!exists("ISOLATE"))        ISOLATE <- TRUE

# ---- stage table ------------------------------------------------------------
# id     : stage number. NOTE: ids match script prefixes only up to 18; from id 19
#          onward they are offset (id 19 = script 18b, id 20 = script 19, ...).
#          RUN_FROM / RUN_ONLY use the id column, not the script number.
# Sentinels revised 19 Oct 2026: each is now a file written at the END of its
# script, so a stage that crashes part-way is not marked done. Previous sentinels
# for stages 1, 2, 4, 6, 7, 18 were directories created at the start of the script;
# stage 3 pointed at a path the script never writes (it reran every time); stage 7
# pointed at a folder created by script 08.
# script : filename prefix; resolved by pattern so minor renames still work
# needs  : "net" = downloads from the internet; "ram" = needs >= 16 GB RAM
# mins   : rough runtime on the development machine
# out    : sentinel file; if it exists the stage is skipped unless FORCE = TRUE
S <- function(id, script, out, mins, needs = "") data.frame(id, script, out, mins, needs,
                                                            stringsAsFactors = FALSE)
STAGES <- rbind(
 S( 0, "00_setup",                "06_results/software_versions.csv",                                   2, ""),
 S( 1, "01_download_GSE132465",   "02_metadata/GSE132465_checksums.csv",                                20, "net"),
 S( 2, "02_pilot_load",           "04_processed_data/GSE132465/GSE132465_full_counts_sparse.rds",       10, ""),
 S( 3, "03_pilot_qc",             "06_results/GSE132465/cell_meta_with_qc_and_compartment.csv.gz",       5, ""),
 S( 4, "04_download_GSE178341",   "02_metadata/GSE178341_checksums.csv",                                45, "net"),
 S( 5, "05_download_spatial",     "02_metadata/GSE166555_checksums.csv",                                90, "net"),
 S( 6, "06_audit_GSE178341",      "06_results/GSE178341_audit/h5_checks.tsv",                            5, ""),
 S( 7, "07_GSE178341_compartment","06_results/GSE178341_audit/donor_compartment_counts_ruleB_unsorted_only.csv", 5, ""),
 S( 8, "08_GSE178341_pseudobulk", "04_processed_data/GSE178341/donor_metadata_tumour.csv",               40, "ram"),
 S( 9, "09_programs_and_gene",    "04_processed_data/programs/gene_partition.tsv",                       10, ""),
 S(10, "10_panel_development",    "04_processed_data/frozen_panel/frozen_panel_v2.rds",                  30, ""),
 S(11, "11_validation_GSE132465", "06_results/validation_GSE132465/validation_results.csv",              20, ""),
 S(12, "12b_KUL_from_saved",      "06_results/validation_GSE144735/validation_results_KUL_and_pooled.csv",20, ""),
 S(13, "13_validation_GSE166555", "06_results/validation_GSE166555/validation_results.csv",              20, ""),
 S(14, "14_GSE200997_own_annot",  "04_processed_data/GSE200997/GSE200997_cell_assignment.rds",           45, "ram"),
 S(15, "15_validation_GSE200997", "06_results/validation_GSE200997/validation_results.csv",              15, ""),
 S(16, "16_cross_cohort_meta",    "06_results/cross_cohort_meta/meta_analysis_results.csv",               5, ""),
 S(17, "17_spatial_reference",    "04_processed_data/spatial_reference/GSE178341_reference_counts.rds",  30, "ram"),
 S(18, "18_spatial_RCTD",         "06_results/spatial_RCTD_section_summary.csv",                        120, ""),
 S(19, "18b_RCTD_label_check",    "06_results/spatial_RCTD_marker_check.csv",                            15, ""),
 S(20, "19_spatial_support",      "06_results/spatial_support/spatial_support_summary.csv",              25, ""),
 S(21, "22_spatial_extra",        "06_results/spatial_support_extra/criterionA_by_section_extra.csv",    90, ""),
 S(22, "23_spatial_combined",     "06_results/spatial_support_combined/spatial_support_summary_combined.csv", 2, ""),
 S(23, "24_discovery_sensitivity","06_results/sensitivity_discovery/discovery_sensitivity_results.csv",  60, ""),
 S(24, "20_evidence_and_chip",    "08_tables/Table_chip_assay_proposal.csv",                              2, ""),
 # Stage 25 (21_figures) removed 19 Oct 2026: its sentinel "07_figures" already
 # exists after scripts 03/11, so it never ran in a clean pipeline; and when forced
 # it overwrote the 8-program figures of stage 38 with 7-program versions under the
 # same *_updated file names. Stage 38 (40_figures_with_lymph) supersedes it.
 S(26, "25_audit_manuscript",     "__always_run__",                               3, ""),
 S(27, "29_lymphoid_pseudobulk", "04_processed_data/GSE178341/pseudobulk_Lymphoid_CD8.rds",            60, "ram"),
 S(28, "30_lymphoid_coverage",   "06_results/programs/lymphoid_coverage.csv",                            2, ""),
 S(29, "31_lymphoid_partition",  "04_processed_data/programs/gene_partition_lymphoid.tsv",               2, ""),
 S(30, "32_lymphoid_panel",      "04_processed_data/frozen_panel/frozen_panel_LYM_v1.rds",              20, ""),
 S(31, "33_lymphoid_validation", "06_results/validation_GSE132465/validation_results_lymphoid.csv",     15, ""),
 S(32, "34_lymphoid_extra",      "06_results/validation_lymphoid/validation_results_extra_cohorts.csv", 25, "ram"),
 S(33, "35_lymphoid_meta",       "06_results/validation_lymphoid/meta_analysis_results_lymphoid.csv",    2, ""),
 S(34, "36_lymphoid_spatial_val","06_results/spatial_support_lymphoid/spatial_support_summary.csv",     25, ""),
 S(35, "37_lymphoid_spatial_ext","06_results/spatial_support_lymphoid/criterionA_by_section_extra.csv", 20, ""),
 S(36, "38_lymphoid_spatial_com","06_results/spatial_support_lymphoid_combined/spatial_support_summary_combined.csv", 2, ""),
 S(37, "39_evidence_and_chip",   "08_tables/Table_chip_assay_proposal_8.csv",                            2, ""),
 S(38, "40_figures_with_lymph",  "07_figures/Figure6_evidence_summary_updated.pdf",                     10, ""),
 S(39, "41_build_submission",    "10_submission/figures/Figure1.pdf",                                    15, ""),
 S(40, "44_build_supplementary", "10_submission/Supplementary_Tables_aem.xlsx",                          2, "")
)

# ---- helpers ----------------------------------------------------------------
resolve <- function(prefix) {
  f <- list.files(SCRIPTS, pattern = paste0("^", prefix, ".*\\.R$"), full.names = TRUE)
  f <- f[!grepl("_before_|_backup|~$", f)]          # ignore patch backups
  if (length(f) == 0) return(NA_character_)
  f[which.min(nchar(f))]
}
done <- function(out) out != "__always_run__" && (file.exists(out) || dir.exists(out))

cat("\n================ CRC biomarker pipeline ================\n")
cat("Project root :", PROJECT_ROOT, "\n")
cat("R version    :", R.version.string, "\n")
cat("Started      :", format(Sys.time()), "\n")
cat("Skip downloads:", SKIP_DOWNLOADS, "| Force rerun:", FORCE, "\n\n")

plan <- STAGES
if (!is.null(RUN_ONLY)) plan <- plan[plan$id %in% RUN_ONLY, ]
plan <- plan[plan$id >= RUN_FROM, ]
if (SKIP_DOWNLOADS)     plan <- plan[plan$needs != "net", ]

cat(sprintf("%-4s %-34s %-9s %s\n", "id", "script", "approx", "status"))
for (i in seq_len(nrow(plan))) {
  st <- plan[i, ]; f <- resolve(st$script)
  status <- if (is.na(f)) "SCRIPT NOT FOUND" else if (!FORCE && done(st$out)) "skip (output exists)" else "to run"
  cat(sprintf("%-4s %-34s %-9s %s\n", st$id, st$script, paste0("~", st$mins, "m"), status))
}
cat("\nTotal estimated runtime if all run:", sum(plan$mins), "minutes\n\n")

if (!exists("DRY_RUN")) DRY_RUN <- FALSE
if (DRY_RUN) { cat("DRY_RUN = TRUE - nothing executed.\n"); } else {

t_all <- Sys.time(); failed <- character(0)
for (i in seq_len(nrow(plan))) {
  st <- plan[i, ]; f <- resolve(st$script)
  if (is.na(f)) { warning("Script not found for stage ", st$id, " (", st$script, ")"); failed <- c(failed, st$script); next }
  if (!FORCE && done(st$out)) next
  cat("\n---- stage", st$id, ":", basename(f), "----\n"); flush.console()
  t0 <- Sys.time()
  if (ISOLATE) {
    rs  <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    lo  <- file.path(LOGDIR, sprintf("stage_%02d.out", st$id)); le <- sub("\\.out$", ".err", lo)
    rc  <- system2(rs, c("--no-save", "--no-restore", shQuote(f)), stdout = lo, stderr = le)
    ok  <- identical(as.integer(rc), 0L)
    if (!ok) { message("STAGE ", st$id, " FAILED (exit ", rc, "). Last lines of ", le, ":")
               if (file.exists(le)) message(paste(tail(readLines(le, warn = FALSE), 15), collapse = "\n")) }
  } else {
    ok <- tryCatch({ source(f, echo = FALSE, local = new.env(parent = globalenv())); TRUE },
                   error = function(e) { message("STAGE ", st$id, " FAILED: ", conditionMessage(e)); FALSE })
  }
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  cat("---- stage", st$id, if (ok) "OK" else "FAILED", "in", mins, "min ----\n")
  # append (writeLines() to a file name overwrote the log at every stage)
  cat(sprintf("%s\tstage %s\t%s\t%s\t%s min\n", format(Sys.time()), st$id, basename(f),
              if (ok) "OK" else "FAILED", mins),
      file = file.path(LOGDIR, "run_all_log.tsv"), append = TRUE)
  if (!ok) { failed <- c(failed, basename(f))
    if (!exists("CONTINUE_ON_ERROR") || !CONTINUE_ON_ERROR) { message("Stopping. Set CONTINUE_ON_ERROR <- TRUE to carry on."); break } }
}
cat("\n================ finished in",
    round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), "minutes ================\n")
if (length(failed)) cat("FAILED STAGES:", paste(failed, collapse = ", "), "\n") else cat("All requested stages completed.\n")
writeLines(capture.output(sessionInfo()), file.path(LOGDIR, paste0("sessionInfo_run_all_", format(Sys.Date()), ".txt")))
}

# -----------------------------------------------------------------------------
# MANUAL STEPS THAT ARE NOT AUTOMATED (see README.md)
#   * Valdeolivas Visium archives (stage 5) must be extracted to
#     04_processed_data/Valdeolivas_Visium/unzipped/SN*/ before stage 17; script 17
#     stops otherwise. No pipeline script performs this extraction - add it
#     (e.g. as 05b_unzip_spatial.R) or document the exact commands in the README.
#   * Likewise the GSE280318 tissue-position files that script 22 looks for under
#     04_processed_data/GSE280318/.
#   * Stage 12 now runs 12b_KUL_from_saved.R, which builds the KUL count matrix
#     itself if it is missing (previously only the unlisted script 12 did this).
#   * GSE226997 (4 archives, ~44 GB) is fetched with 03_original_data/GSE226997/
#     download_GSE226997.bat and then extracted; stage 21 expects the extracted
#     files under 04_processed_data/GSE226997/Ajou_Visium_P1..P4/.
#   * Stages 18 and 19 (scripts 18 and 18b) are a PAIR: 18 writes RCTD weights
#     without dimnames and 18b verifies and restores them. Script 19 refuses to
#     run unless 18b has set labels_verified = TRUE.
#   * Disk: ~150 GB free is needed with all raw archives present. NTFS required
#     (several files exceed the 4 GB FAT32 limit).

#   * Stages 27-39 add the lymphoid (CD8) compartment (protocol v1.8-v1.9). They are
#     additive: frozen_panel_v2.rds is never modified, and stages 38-39 supersede the
#     seven-program figures and submission bundle produced by stages 25 and 24.
# -----------------------------------------------------------------------------
