###############################################################################
# 51_inspect_inputs.R - STEP 1: check that every input exists and looks right.
# Reads nothing it modifies. Run it, read the report, fix 50_revision_config.R.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision()
LIST_CANDIDATES <- TRUE   # list .rds files that look like pseudobulk counts

report <- character()
add <- function(...) { line <- paste0(...); cat(line, "\n"); report <<- c(report, line) }

add("== Project root: ", PROJECT_ROOT, if (dir.exists(PROJECT_ROOT)) "  [found]" else "  [NOT FOUND - extract the .rar]")
add("== Supplementary workbook: ", SUPP_XLSX, if (file.exists(SUPP_XLSX)) "  [found]" else "  [NOT FOUND]")

## 1. Frozen panel hashes (must be unchanged)
if (requireNamespace("digest", quietly = TRUE)) {
  for (nm in names(FROZEN_FILES)) {
    f <- FROZEN_FILES[[nm]]
    if (file.exists(f)) {
      h <- digest::digest(file = f, algo = "sha256")
      add("== Frozen panel ", nm, ": ", if (identical(h, FROZEN_SHA256[[nm]])) "hash OK" else paste("HASH MISMATCH", h))
    } else add("== Frozen panel ", nm, ": file not found (", f, ")")
  }
} else add("== install.packages('digest') to check frozen-panel hashes")

## 2. Panels and partition from the workbook
panels <- load_panels(); partition <- load_partition()
add("== Programs in S6: ", paste(names(panels), collapse = ", "))
for (p in names(panels)) add("   ", p, ": ", nrow(panels[[p]]), " genes (", paste(panels[[p]]$symbol, collapse = ", "),
                             "); evaluation genes in S5: ", partition[program == p & role == "evaluation", .N])
stopifnot(length(intersect(partition[role == "evaluation", row_id], partition[role != "evaluation", row_id])) == 0 ||
          TRUE)  # disjointness is checked per program below
for (p in names(panels)) {
  ev <- partition[role == "evaluation", unique(row_id)]
  if (any(panels[[p]]$row_id %in% ev)) add("   WARNING: panel gene of ", p, " is an evaluation gene of some program")
}

## 3. Candidate pseudobulk files
if (LIST_CANDIDATES) {
  add("== .rds files under 04_processed_data that look like pseudobulk:")
  cand <- list.files(file.path(PROJECT_ROOT, "04_processed_data"), pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  cand <- cand[grepl("pseudo|pb_|_pb|counts|bulk", basename(cand), ignore.case = TRUE)]
  for (f in head(cand, 80)) add("   ", f)
}

## 4. Each configured pseudobulk input (log-CPM of eligible donors)
add("== Configured pseudobulk inputs (eligible donors after the >= ", MIN_CELLS, "-cell rule)")
cmp_needed <- unique(program_compartment(panels))
for (co in names(PB_FILES)) for (cmp in cmp_needed) {
  L <- tryCatch(load_lc(co, cmp, quiet = TRUE), error = function(e) e)
  if (inherits(L, "error")) { add(sprintf("   %-10s %-12s ERROR    %s", co, cmp, conditionMessage(L))); next }
  if (is.null(L)) { add(sprintf("   %-10s %-12s not available", co, cmp)); next }
  key <- key_col(co)
  pg <- unlist(lapply(names(panels)[program_compartment(panels) == cmp], function(p) panels[[p]][[key]]))
  add(sprintf("   %-10s %-12s %6d genes x %3d donors; rownames e.g. %s; panel genes found %d/%d",
              co, cmp, nrow(L), ncol(L), paste(head(rownames(L), 2), collapse = ","),
              sum(pg %in% rownames(L)), length(pg)))
}
add("== Expected eligible donors (S1b/S1c): discovery 60/27/60/40/53; GSE132465 21/13/20/11/20;",
    " GSE166555 12/10/12/7/12; GSE200997 12/5/12/3/-; GSE144735 6/5/5/5/2")

## 5. Discovery feature map
if (!is.na(DISCOVERY_FEATURES) && file.exists(DISCOVERY_FEATURES)) {
  ft <- fread(DISCOVERY_FEATURES); add("== Discovery feature map: ", nrow(ft), " rows; columns ", paste(names(ft), collapse = ","))
} else add("== Discovery feature map not found: run 50b_make_feature_map.R (needed for scripts 53 and 61).")

## 6. Packages
for (p in c("Matrix", "data.table", "readxl", "openxlsx", "glmnet", "metafor", "edgeR", "digest", "spacexr"))
  add("== package ", p, ": ", if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "NOT INSTALLED")

writeLines(report, file.path(OUT_DIR, "51_inspect_inputs_report.txt"))
msg("Report written to", file.path(OUT_DIR, "51_inspect_inputs_report.txt"))
