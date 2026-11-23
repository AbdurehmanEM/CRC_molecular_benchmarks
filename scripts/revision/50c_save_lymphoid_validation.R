###############################################################################
# 50c_save_lymphoid_validation.R - saves the CD8 (lymphoid) validation log-CPM
# for GSE166555 (and GSE144735), which script 34 computed but did not save.
# It runs your 34_lymphoid_extra_cohorts.R in a separate environment, then looks
# for gene x donor matrices whose donor IDs belong to each cohort, and writes
# <cohort>_Lymphoid_logcpm.rds to EXTRA_LYMPHOID_DIR. Nothing of script 34's own
# outputs is changed (it simply runs again). Optional: only needed so that the
# reproduction gate can check the published CD8 GSE166555 result (0.214).
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision()
dir.create(EXTRA_LYMPHOID_DIR, recursive = TRUE, showWarnings = FALSE)

s34 <- list.files(file.path(PROJECT_ROOT, "05_analysis_scripts"), pattern = "^34_.*\\.R$", full.names = TRUE)
s34 <- s34[!grepl("_before_", s34)][1]
if (is.na(s34)) stop("Script 34 not found in 05_analysis_scripts")
msg("Running", basename(s34), "(this repeats its computation; outputs are identical)")
env <- new.env()
old <- setwd(PROJECT_ROOT); on.exit(setwd(old), add = TRUE)
sys.source(s34, envir = env)
setwd(old)

donor_pat <- c(GSE166555 = "^p0\\d\\d$", GSE144735 = "^KUL\\d+$")
found <- list()
scan <- function(x, path) {
  if ((is.matrix(x) || inherits(x, "Matrix")) && nrow(x) > 5000 && !is.null(colnames(x))) {
    for (co in names(donor_pat)) if (all(grepl(donor_pat[[co]], colnames(x))))
      found[[length(found) + 1]] <<- list(cohort = co, path = path, x = x)
  } else if (is.list(x) && !is.data.frame(x) && length(x) < 50) {
    for (nm in names(x)) scan(x[[nm]], paste0(path, "$", nm))
  }
}
for (nm in ls(env)) scan(get(nm, envir = env), nm)
if (!length(found)) stop("No gene x donor matrix for GSE166555/GSE144735 found in script 34's objects. ",
                         "Send the output of ls() after sourcing script 34.")
for (f in found) {
  x <- as.matrix(f$x)
  is_counts <- all(abs(x[1:min(200, nrow(x)), ] - round(x[1:min(200, nrow(x)), ])) < 1e-8) && max(x) > 50
  lc <- if (is_counts) logcpm(x) else x
  msg(f$cohort, ": object", f$path, "|", nrow(x), "genes x", ncol(x), "donors |", if (is_counts) "counts -> log-CPM" else "log-CPM")
  saveRDS(list(logcpm = lc, ncells = NULL, source = paste(basename(s34), f$path), created = format(Sys.time())),
          file.path(EXTRA_LYMPHOID_DIR, paste0(f$cohort, "_Lymphoid_logcpm.rds")))
}
msg("Saved to", EXTRA_LYMPHOID_DIR, ". If more than one object per cohort was found, the last one was kept;",
    "check that its donors (12 for GSE166555) match S1c.")
