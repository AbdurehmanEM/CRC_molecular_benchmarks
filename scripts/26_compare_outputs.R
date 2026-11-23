# 26_compare_outputs.R - compare two runs of the pipeline file by file
# Usage:
#   PROJECT_A <- "F:/.../CRC_Biomarker_Project"          # reference (published) run
#   PROJECT_B <- "F:/.../CRC_Biomarker_Project_rerun"    # clean re-run
#   source("05_analysis_scripts/26_compare_outputs.R")
suppressPackageStartupMessages({ library(data.table); library(digest) })
stopifnot(exists("PROJECT_A"), exists("PROJECT_B"), dir.exists(PROJECT_A), dir.exists(PROJECT_B))
DIRS <- c("06_results", "08_tables", "04_processed_data")
TOL  <- 1e-6          # numeric tolerance for CSV comparison

lsfiles <- function(root) {
  f <- list.files(file.path(root, DIRS), recursive = TRUE, full.names = FALSE)
  f <- f[!grepl("run_logs|sessionInfo|\\.txt$|\\.pdf$|\\.tiff$|\\.png$", f)]
  sort(unique(f))
}
A <- lsfiles(PROJECT_A); B <- lsfiles(PROJECT_B)
cat("Files in reference run:", length(A), "| in re-run:", length(B), "\n")
only_a <- setdiff(A, B); only_b <- setdiff(B, A)
if (length(only_a)) { cat("\nMISSING from the re-run (", length(only_a), "):\n"); print(head(only_a, 40)) }
if (length(only_b)) { cat("\nEXTRA in the re-run (", length(only_b), "):\n"); print(head(only_b, 40)) }

cmp_csv <- function(pa, pb) {
  x <- tryCatch(fread(pa), error = function(e) NULL); y <- tryCatch(fread(pb), error = function(e) NULL)
  if (is.null(x) || is.null(y)) return("unreadable")
  if (!identical(dim(x), dim(y)))      return(sprintf("dim %s vs %s", paste(dim(x), collapse = "x"), paste(dim(y), collapse = "x")))
  if (!identical(names(x), names(y)))  return("different column names")
  num <- names(x)[vapply(x, is.numeric, logical(1))]
  chr <- setdiff(names(x), num)
  bad <- character(0)
  for (k in chr) if (!identical(as.character(x[[k]]), as.character(y[[k]]))) bad <- c(bad, k)
  worst <- 0
  for (k in num) { d <- suppressWarnings(max(abs(x[[k]] - y[[k]]), na.rm = TRUE))
                   if (is.finite(d)) worst <- max(worst, d)
                   if (!identical(is.na(x[[k]]), is.na(y[[k]]))) bad <- c(bad, paste0(k, " (NA pattern)")) }
  if (length(bad)) return(paste("differs in:", paste(unique(bad), collapse = ", ")))
  if (worst > TOL) return(sprintf("numeric diff up to %.3g", worst))
  if (worst > 0)   return(sprintf("identical within tolerance (max %.3g)", worst))
  "identical"
}

common <- intersect(A, B); res <- list()
for (f in common) {
  pa <- file.path(PROJECT_A, f); pb <- file.path(PROJECT_B, f)
  v <- if (grepl("\\.(csv|tsv)$", f)) cmp_csv(pa, pb)
       else if (digest(file = pa) == digest(file = pb)) "identical (hash)"
       else "BINARY DIFFERS"
  res[[f]] <- data.table(file = f, result = v)
}
R <- rbindlist(res)
R[, status := fifelse(grepl("^identical", result), "same", "DIFFERENT")]
cat("\n================ SUMMARY ================\n")
print(R[, .N, by = status])
cat("\n-- files that differ --\n")
print(R[status == "DIFFERENT"], nrows = 100)
fwrite(R, file.path(PROJECT_A, "06_results", "run_logs", "compare_outputs.csv"))
cat("\nFull report: 06_results/run_logs/compare_outputs.csv\n")
cat("\nNOTE: frozen_panel_v2.rds and any .rds storing a creation time will differ by hash\n",
    "even when their contents match. Compare those with readRDS() if flagged.\n")
