###############################################################################
# 63_reproducibility_clean_room.R - Reviewer point 27 (reproducibility)
#
# Part 1 (always): input manifest with SHA-256, size and date of every file in
#   03_original_data and of the frozen panels; R and package versions.
# Part 2 (RUN_CLEAN_ROOM <- TRUE): copies the scripts, original data and the
#   frozen panels to CLEAN_ROOM_DIR, runs 00_run_all.R from stage
#   CLEAN_ROOM_RUN_FROM (11 keeps the supplied frozen panel; re-running stage 10
#   would rewrite it and change its hash), then compares every CSV in
#   06_results between the original and the clean-room run.
# Needs: digest. The clean-room run takes ~12 h and several GB of disk.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision()
RUN_CLEAN_ROOM <- FALSE
NUM_TOL <- 1e-6

## ---- Part 1: manifest -------------------------------------------------------------
sha <- function(f) digest::digest(file = f, algo = "sha256")
inp <- list.files(file.path(PROJECT_ROOT, "03_original_data"), recursive = TRUE, full.names = TRUE)
inp <- c(inp, unname(FROZEN_FILES[file.exists(FROZEN_FILES)]))
man <- data.table(file = sub(paste0("^", PROJECT_ROOT, "/?"), "", inp), bytes = file.size(inp),
                  modified = format(file.mtime(inp), "%Y-%m-%d %H:%M"))
man[, sha256 := vapply(inp, sha, character(1))]
write_out(man, "63_input_manifest_sha256")
pk <- as.data.table(installed.packages()[, c("Package", "Version")])
pk <- pk[Package %in% c("Matrix", "data.table", "edgeR", "glmnet", "metafor", "Seurat", "SeuratObject", "spacexr",
                        "SpatialExperiment", "rhdf5", "hdf5r", "msigdbr", "readxl", "openxlsx", "digest", "ggplot2")]
pk <- rbind(data.table(Package = "R", Version = R.version.string), pk)
write_out(pk, "63_software_versions_revision")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "63_sessionInfo.txt"))

## ---- Part 2: clean room ----------------------------------------------------------
if (RUN_CLEAN_ROOM) {
  CR <- CLEAN_ROOM_DIR
  if (dir.exists(CR)) stop(CR, " exists; remove it or choose another CLEAN_ROOM_DIR")
  dir.create(CR, recursive = TRUE)
  copy_dir <- function(from, to) { dir.create(to, recursive = TRUE, showWarnings = FALSE)
    file.copy(list.files(from, full.names = TRUE, all.files = FALSE), to, recursive = TRUE) }
  for (d in c("01_protocol", "02_metadata", "03_original_data", "05_analysis_scripts")) copy_dir(file.path(PROJECT_ROOT, d), file.path(CR, d))
  dir.create(file.path(CR, "04_processed_data", "frozen_panel"), recursive = TRUE)
  file.copy(FROZEN_FILES[file.exists(FROZEN_FILES)], file.path(CR, "04_processed_data", "frozen_panel"))
  # Stages before RUN_FROM produce intermediate files that later stages need: copy those folders too.
  # ADAPT if your stage table differs: everything in 04_processed_data except revision folders.
  copy_dir(file.path(PROJECT_ROOT, "04_processed_data"), file.path(CR, "04_processed_data"))
  runner <- readLines(file.path(CR, "05_analysis_scripts", "00_run_all.R"))
  i <- grep("^\\s*RUN_FROM\\s*<-", runner)
  if (!length(i)) stop("RUN_FROM assignment not found in 00_run_all.R")
  runner[i[1]] <- paste0("RUN_FROM <- ", CLEAN_ROOM_RUN_FROM)
  writeLines(runner, file.path(CR, "05_analysis_scripts", "00_run_all_cleanroom.R"))
  msg("Running clean-room pipeline from stage", CLEAN_ROOM_RUN_FROM, "in", CR, "(long)")
  old_wd <- setwd(CR)                                    # pipeline scripts use paths relative to the project root
  st <- system2(file.path(R.home("bin"), "Rscript"), c("05_analysis_scripts/00_run_all_cleanroom.R"),
                stdout = file.path(CR, "cleanroom_stdout.txt"), stderr = file.path(CR, "cleanroom_stderr.txt"),
                wait = TRUE, env = paste0("CRC_PROJECT_ROOT=", CR))
  setwd(old_wd)
  msg("Clean-room exit status", st)

  ## compare CSV outputs
  orig <- list.files(file.path(PROJECT_ROOT, "06_results"), pattern = "\\.csv$", recursive = TRUE)
  orig <- orig[!grepl("revision_", orig)]
  cmp <- rbindlist(lapply(orig, function(f) {
    a <- file.path(PROJECT_ROOT, "06_results", f); b <- file.path(CR, "06_results", f)
    if (!file.exists(b)) return(data.table(file = f, status = "missing in clean room", max_abs_diff = NA_real_))
    x <- fread(a); y <- fread(b)
    if (!identical(dim(x), dim(y)) || !identical(names(x), names(y))) return(data.table(file = f, status = "different shape", max_abs_diff = NA_real_))
    nc <- names(x)[sapply(x, is.numeric)]
    md <- if (length(nc)) max(sapply(nc, function(c) max(abs(x[[c]] - y[[c]]), na.rm = TRUE)), 0) else 0
    oc <- setdiff(names(x), nc); same_txt <- all(sapply(oc, function(c) identical(as.character(x[[c]]), as.character(y[[c]]))))
    data.table(file = f, status = if (md <= NUM_TOL && same_txt) "identical" else if (same_txt) "numeric differences" else "text differences",
               max_abs_diff = md)
  }))
  write_out(cmp, "63_clean_room_comparison")
  print(cmp[, .N, by = status])
}
