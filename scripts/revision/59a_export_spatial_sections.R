###############################################################################
# 59a_export_spatial_sections.R - writes one file per spatial section in the
# format used by 59_spatial_rerun.R and 60_independent_rctd_reference.R:
#   list(section, patient, cohort, grid, counts, coords, weights = list(published, randomsel))
# Sources (your existing outputs; nothing is re-run):
#   04_processed_data/spatial_extra/*_rctd.rds            (GSE226997, GSE280318)
#   04_processed_data/spatial_extra_randomsel/*_rctd.rds  (random-selection reference)
#   04_processed_data/Valdeolivas_Visium/rctd/*_rctd.rds  (+ rctd_randomsel/)
# Each *_rctd.rds is read as either
#   (a) a list with $weights, $cnt, $pos (as in spatial_extra), or
#   (b) a spacexr RCTD object (counts/coords taken from its spatialRNA slot), with
#       array coordinates from the Space Ranger tissue_positions file if present.
# Sections with role != "main" (second section of a patient already included) and
# patient A798015 (non-neoplastic) are skipped, as in the published analysis.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision()
dir.create(SPATIAL_SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

src <- data.table(file = c(list.files(file.path(PB_DIR, "spatial_extra"), "_rctd\\.rds$", full.names = TRUE),
                           list.files(file.path(PB_DIR, "Valdeolivas_Visium", "rctd"), "_rctd\\.rds$", full.names = TRUE)))
src[, section := sub("_rctd\\.rds$", "", basename(file))]
src[, rand_file := sub("/spatial_extra/", "/spatial_extra_randomsel/", sub("/rctd/", "/rctd_randomsel/", file, fixed = TRUE), fixed = TRUE)]
msg(nrow(src), "section files found")

as_weights <- function(w) { w <- as.matrix(w); w / rowSums(w) }
find_positions <- function(section) {
  root <- file.path(PB_DIR, "Valdeolivas_Visium", "unzipped")
  f <- list.files(root, pattern = "tissue_positions.*\\.csv$", recursive = TRUE, full.names = TRUE)
  f <- f[grepl(section, f, fixed = TRUE)]
  if (!length(f)) return(NULL)
  p <- fread(f[1], header = grepl("list", basename(f[1])) == FALSE)
  if (ncol(p) >= 6 && !"barcode" %in% names(p)) setnames(p, 1:6, c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  p
}
grid_from_pixels <- function(xy) {         # fallback: convert pixel coordinates to an integer grid
  d <- median(apply(as.matrix(dist(xy[sample(nrow(xy), min(500, nrow(xy))), ])), 1, function(r) min(r[r > 0])))
  data.frame(array_row = round((xy[, 2] - min(xy[, 2])) / (d * sqrt(3) / 2)),
             array_col = round((xy[, 1] - min(xy[, 1])) / (d / 2)))
}

elog <- list()
for (i in seq_len(nrow(src))) {
  s <- src[i]; o <- readRDS(s$file)
  role <- if (is.list(o) && !is.null(o$role)) o$role else "main"
  patient <- if (is.list(o) && !is.null(o$patient)) o$patient else sub("^SN\\d+_(A\\d+)_.*$", "\\1", s$section)
  cohort <- if (is.list(o) && !is.null(o$cohort)) o$cohort else "Valdeolivas"
  platform <- if (is.list(o) && !is.null(o$platform)) o$platform else "Visium"
  if (patient %in% SPATIAL_EXCLUDE_PATIENTS) { elog[[i]] <- data.table(section = s$section, status = "skipped: excluded patient"); next }
  if (!identical(role, "main")) { elog[[i]] <- data.table(section = s$section, status = paste("skipped: role", role)); next }

  if (is.list(o) && !isS4(o) && !is.null(o$weights) && !is.null(o$cnt)) {          # format (a)
    W <- as_weights(o$weights); cnt <- o$cnt
    pos <- as.data.frame(o$pos)
    crd <- data.frame(spot = pos$barcode, array_row = pos$array_row, array_col = pos$array_col)
  } else {                                                                          # format (b): spacexr RCTD
    obj <- if (isS4(o)) o else (o$rctd %||% o$RCTD %||% o[[which(sapply(o, isS4))[1]]])
    W <- as_weights(obj@results$weights %||% obj@results$weights_doublet)
    cnt <- obj@spatialRNA@counts
    if (is.list(o) && !isS4(o) && !is.null(o$weights)) W <- as_weights(o$weights)
    p <- find_positions(s$section)
    if (!is.null(p)) {
      crd <- data.frame(spot = p$barcode, array_row = p$array_row, array_col = p$array_col)
    } else {
      xy <- as.matrix(obj@spatialRNA@coords[, 1:2])
      crd <- cbind(spot = rownames(obj@spatialRNA@coords), grid_from_pixels(xy))
      msg("  ", s$section, ": no tissue_positions file; grid derived from pixel coordinates")
    }
  }
  spots <- intersect(rownames(W), colnames(cnt)); spots <- intersect(spots, crd$spot)
  weights <- list(published = W[spots, , drop = FALSE])
  if (file.exists(s$rand_file)) {
    r <- readRDS(s$rand_file)
    Wr <- if (is.list(r) && !isS4(r) && !is.null(r$weights)) r$weights else
          (if (isS4(r)) r@results$weights else r[[which(sapply(r, isS4))[1]]]@results$weights)
    Wr <- as_weights(Wr); weights$randomsel <- Wr[intersect(spots, rownames(Wr)), , drop = FALSE]
  }
  out <- list(section = s$section, patient = patient, cohort = cohort,
              grid = if (grepl("HD", platform, ignore.case = TRUE)) "square" else "hex",
              counts = as(cnt[, spots, drop = FALSE], "CsparseMatrix"),
              coords = crd[match(spots, crd$spot), ], weights = weights)
  saveRDS(out, file.path(SPATIAL_SECTION_DIR, paste0(s$section, ".rds")))
  elog[[i]] <- data.table(section = s$section, status = "exported", cohort = cohort, patient = patient,
                         grid = out$grid, spots = length(spots), genes = nrow(cnt),
                         types = ncol(W), randomsel = !is.null(weights$randomsel))
  msg("exported", s$section, "|", cohort, patient, "|", length(spots), "spots")
}
elog <- rbindlist(elog, fill = TRUE)
write_out(elog, "59a_export_log")
print(elog)
cat("\nExpected: 13 patients (Valdeolivas 6, GSE226997 4, GSE280318 3). Patients exported:",
    uniqueN(elog[status == "exported", patient]), "\n")
