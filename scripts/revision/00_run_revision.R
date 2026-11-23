###############################################################################
# 00_run_revision.R - runs the revision analyses in order, each in a fresh R
# process, with per-script logs. Put all 5x/6x scripts and the two 50_ files in
# 05_analysis_scripts/revision/ and run from that folder:
#   setwd("F:/CRC_Biomarker_Project 2/05_analysis_scripts/revision"); source("00_run_revision.R")
# The workbook-only scripts (55, 57, 58) run even if pseudobulk files are not set.
###############################################################################
STEPS <- c(
  "51_inspect_inputs.R",              # 1  check inputs                       (all points)
  "50c_save_lymphoid_validation.R",   # 1b optional: save GSE166555 CD8 log-CPM (re-runs script 34)
  "52_reproduce_published.R",         # 2  GATE: reproduce S8/S8b             (prerequisite)
  "55_meta_analysis_sensitivity.R",   # 3  REML/HK, leave-one-cohort-out      (point 6)   workbook only
  "57_multiplicity_fdr.R",            # 4  BH/Holm                            (point 15)  workbook (+52)
  "58_spatial_table_sensitivity.R",   # 5  thresholds, excl. GSE226997        (points 11, 14) workbook only
  "54_leave_one_donor_out.R",         # 6  LODO                               (point 4)
  "56_calibration_transport.R",       # 7  offsets/slopes, recalibration      (point 8)
  "50b_make_feature_map.R",           # 8a discovery feature map (needs rhdf5 or hdf5r)
  "53_null_random_panels.R",          # 8  random-gene nulls                  (point 7)   ~1-2 h
  "61_alternative_panels.R",          # 9  alternative panels                 (point 13)
  "62_mmr_site_strata.R",             # 10 MMR / site                         (point 26)
  "59a_export_spatial_sections.R",     # 11a export spatial sections (reads your RCTD files)
  "60_independent_rctd_reference.R",  # 11 independent RCTD reference         (point 10)  needs spacexr + section export
  "59_spatial_rerun.R",               # 12 block size, negative controls      (points 9, 10)
  "63_reproducibility_clean_room.R",  # 13 manifest (+ optional clean room)   (point 27)
  "64_collect_revision_tables.R")     # 14 assemble tables + numbers
SKIP <- c("50c_save_lymphoid_validation.R")   # remove from SKIP to check the published CD8 GSE166555 result             # e.g. c("60_independent_rctd_reference.R", "59_spatial_rerun.R") until sections are exported
STOP_ON_ERROR <- FALSE  # TRUE: stop at the first failing script

source("50_revision_config.R")
dir.create("logs", showWarnings = FALSE)
rs <- file.path(R.home("bin"), "Rscript")
status <- data.frame(script = STEPS, exit = NA_integer_, minutes = NA_real_)
for (i in seq_along(STEPS)) {
  s <- STEPS[i]
  if (s %in% SKIP) { cat(sprintf("[%2d] %-36s skipped\n", i, s)); next }
  t0 <- Sys.time()
  ex <- system2(rs, s, stdout = file.path("logs", sub("\\.R$", ".out", s)), stderr = file.path("logs", sub("\\.R$", ".err", s)))
  status$exit[i] <- ex; status$minutes[i] <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  cat(sprintf("[%2d] %-36s exit %d  (%.1f min)\n", i, s, ex, status$minutes[i]))
  if (s == "52_reproduce_published.R") {
    g <- tryCatch(readLines(file.path(OUT_DIR, "52_gate_status.txt"))[1], error = function(e) NA)
    cat("     gate:", g, "\n")
  }
  if (ex != 0) { cat("     see logs/", sub("\\.R$", ".err", s), "\n", sep = "")
                 if (STOP_ON_ERROR) break }
}
write.csv(status, file.path("logs", "run_revision_status.csv"), row.names = FALSE)
