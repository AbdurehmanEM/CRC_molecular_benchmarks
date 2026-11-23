# 42_spatial_sensitivity_random_reference.R  (sensitivity analysis, added 20 Oct 2026)
#
# Why: 17b_check_reference_donor_balance.R showed that the RCTD reference built by
# script 17 favours donors stored early in the GSE178341 h5 file (rho -0.47 to -0.73
# for every cell type; 13-27 donors per type absent). This script rebuilds the
# reference with random selection (FIX_SELECTION_ORDER = TRUE) and repeats the whole
# spatial analysis with it.
#
# Nothing published is overwritten. The original scripts 17, 18, 18b, 19, 22, 23 and
# 36, 37, 38 are read as text, every input/output path is redirected to *_randomsel
# or 06_results/sens_randomref/, and the redirected code is run. The script stops if
# any published path would still be written.
#
# Runtime: ~4-5 h (RCTD on 12 Valdeolivas + 8 extra sections). Re-sourcing resumes:
# finished RCTD sections are skipped.
# Run from the project root:  source("05_analysis_scripts/42_spatial_sensitivity_random_reference.R")

stopifnot(dir.exists("05_analysis_scripts"))
S <- "05_analysis_scripts"
dir.create("06_results/sens_randomref", recursive = TRUE, showWarnings = FALSE)

# Order matters: longer, more specific strings first.
REDIRECT <- c(
  "GSE178341_reference_counts.rds"              = "GSE178341_reference_counts_randomsel.rds",
  "04_processed_data/Valdeolivas_Visium/rctd"   = "04_processed_data/Valdeolivas_Visium/rctd_randomsel",
  "04_processed_data/spatial_extra"             = "04_processed_data/spatial_extra_randomsel",
  "spot_program_scores"                         = "spot_program_scores_randomsel",
  "06_results/spatial_RCTD_"                    = "06_results/sens_randomref/spatial_RCTD_",
  "06_results/spatial_support"                  = "06_results/sens_randomref/spatial_support")
# Published locations that must not appear in any redirected script:
FORBIDDEN <- c('reference_counts.rds"', 'Valdeolivas_Visium/rctd"', 'spatial_extra"',
               'spot_program_scores.rds"', 'spot_program_scores_lymphoid.rds"',
               '"06_results/spatial_RCTD_', '"06_results/spatial_support')

find_script <- function(prefix) {
  f <- list.files(S, pattern = paste0("^", prefix, ".*\\.R$"), full.names = TRUE)
  f <- f[!grepl("_before_|_backup", f)]
  if (length(f) != 1) stop("Expected exactly one script for prefix ", prefix, ", found: ", paste(f, collapse = ", "))
  f
}
run_redirected <- function(prefix, extra = list()) {
  f <- find_script(prefix); txt <- readLines(f, warn = FALSE)
  for (k in names(REDIRECT)) txt <- gsub(k, REDIRECT[[k]], txt, fixed = TRUE)
  txt <- gsub("_randomsel_randomsel", "_randomsel", txt, fixed = TRUE)          # guard double substitution
  txt <- gsub("sens_randomref/sens_randomref/", "sens_randomref/", txt, fixed = TRUE)
  # All stages share one session here, so SingleCellExperiment (attached by script 18
  # via SpatialExperiment) masks edgeR::cpm with an S4 generic that rejects matrices.
  # Qualify every bare cpm( call. Same function, same result.
  txt <- gsub("(^|[^:A-Za-z0-9_.])cpm\\(", "\\1edgeR::cpm(", txt, perl = TRUE)
  bad <- FORBIDDEN[vapply(FORBIDDEN, function(p) any(grepl(p, txt, fixed = TRUE)), logical(1))]
  if (length(bad)) stop(basename(f), ": published path still present after redirection: ", paste(bad, collapse = ", "))
  out <- file.path("06_results/sens_randomref", paste0("redirected_", basename(f)))
  writeLines(txt, out)                                                         # kept for inspection
  cat("\n================ running", basename(f), "(redirected) ================\n")
  env <- new.env(parent = globalenv()); for (n in names(extra)) assign(n, extra[[n]], envir = env)
  eval(parse(text = txt), envir = env)
  invisible(TRUE)
}

# 1. reference with random selection (script 17 already writes ..._randomsel.rds when the flag is TRUE)
if (!file.exists("04_processed_data/spatial_reference/GSE178341_reference_counts_randomsel.rds")) {
  env17 <- new.env(parent = globalenv()); assign("FIX_SELECTION_ORDER", TRUE, envir = env17)
  sys.source(find_script("17_spatial_reference"), envir = env17)
} else cat("Random-selection reference already built - reusing\n")

# 2. Valdeolivas: RCTD, label check, spatial support
run_redirected("18_spatial_RCTD")
run_redirected("18b_RCTD_label_check")
if (!all(vapply(list.files("04_processed_data/Valdeolivas_Visium/rctd_randomsel", "_rctd.rds$", full.names = TRUE),
                function(p) isTRUE(readRDS(p)$labels_verified), logical(1))))
  stop("18b did not verify the labels for the random-selection RCTD run - inspect its output before continuing")
run_redirected("19_spatial_support")

# 3. extra cohorts and pooling
run_redirected("22_spatial_extra")
run_redirected("23_spatial_combined")

# 4. lymphoid spatial (descriptive)
run_redirected("36_lymphoid_spatial")
run_redirected("37_lymphoid_spatial")
run_redirected("38_lymphoid_spatial")

# 5. comparison: published vs random-selection reference
library(data.table)
rd <- function(p) fread(p)[analysis == "All spatial cohorts (main)",
                           .(program, n_patients, pooled_effect, lo, hi, A_positive, B_positive, spatial_label)]
pub <- rbind(rd("06_results/spatial_support_combined/spatial_support_summary_combined.csv"),
             rd("06_results/spatial_support_lymphoid_combined/spatial_support_summary_combined.csv"))
sen <- rbind(rd("06_results/sens_randomref/spatial_support_combined/spatial_support_summary_combined.csv"),
             rd("06_results/sens_randomref/spatial_support_lymphoid_combined/spatial_support_summary_combined.csv"))
cmp <- merge(pub, sen, by = "program", suffixes = c("_published", "_randomref"))
cmp[, label_changed := spatial_label_published != spatial_label_randomref]
# Spatial label only affects the evidence category of programs that passed independent
# single-cell validation (END_TNFA_NFKB, EPI_HYPOXIA under the v1.3 rule).
cmp[, affects_category := program %in% c("END_TNFA_NFKB", "EPI_HYPOXIA") & label_changed]
fwrite(cmp, "06_results/sens_randomref/comparison_published_vs_randomref.csv")
cat("\n================ PUBLISHED vs RANDOM-SELECTION REFERENCE (all spatial cohorts, main) ================\n")
print(cmp[, .(program, effect_pub = pooled_effect_published, ci_pub = sprintf("%.3f-%.3f", lo_published, hi_published),
              effect_rand = pooled_effect_randomref, ci_rand = sprintf("%.3f-%.3f", lo_randomref, hi_randomref),
              label_pub = spatial_label_published, label_rand = spatial_label_randomref, label_changed, affects_category)])
cat("\nNo spatial label changes affecting an evidence category:", !any(cmp$affects_category), "\n")

dir.create("06_results/sessionInfo", recursive = TRUE, showWarnings = FALSE); writeLines(capture.output(sessionInfo()), "06_results/sessionInfo/sessionInfo_42.txt")
