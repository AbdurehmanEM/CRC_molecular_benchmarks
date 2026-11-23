###############################################################################
# 60_independent_rctd_reference.R - Reviewer point 10 (reference circularity)
#
# The published RCTD reference was built from GSE178341 (the discovery cohort).
# This script builds an alternative reference from GSE132465 single cells
# (Lee et al. 2020; independent of discovery), with the same nine cell types,
# re-runs RCTD on every exported section (from script 59, section 0), and adds
# the weights as weights$independent to each section file. Then re-run 59.
#
# NOTE: GSE132465 is also a single-cell validation cohort; using it here is
# independent of panel DEVELOPMENT (discovery), which is the circularity the
# reviewer raised. State this in the Methods.
# The spacexr calls below follow the Bioconductor API used by your scripts
# 18/22 (createRctd/runRctd). If your script 18 passes extra arguments
# (UMI_min, max_cores, rctd_mode), copy them into RCTD_ARGS so both runs match.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision()
suppressPackageStartupMessages({ library(spacexr); library(SummarizedExperiment); library(SpatialExperiment) })
RCTD_MODE <- "full"; RCTD_CORES <- 4

## ---- 1. Annotation and cell-type mapping (printed for checking) -----------------
ann <- fread(GSE132465_ANNOT)
setnames(ann, names(ann)[1], "Index")
ann[, ref_type := fcase(
  Cell_type == "Epithelial cells" & Class == "Tumor" & grepl("^CMS", Cell_subtype), "Malignant_epithelial",
  Cell_type == "Epithelial cells" & Class == "Normal",                                "Normal_epithelial",
  grepl("EC|[Ee]ndothel", Cell_subtype),                                              "Endothelial",
  grepl("Pericyte|Smooth muscle", Cell_subtype),                                      "Pericyte_SMC",
  Cell_type == "Stromal cells" & grepl("Stromal|[Ff]ibroblast|Myofibro", Cell_subtype), "Fibroblast",
  Cell_type == "Myeloids",                                                            "Myeloid",
  Cell_type == "T cells",                                                             "T_NK",
  Cell_type == "B cells",                                                             "B_Plasma",
  Cell_type == "Mast cells",                                                          "Mast",
  default = NA_character_)]
map_tab <- ann[, .N, by = .(Cell_type, Cell_subtype, Class, ref_type)][order(ref_type, Cell_type)]
write_out(map_tab, "60_reference_celltype_mapping")
print(map_tab, nrows = 200)
cat("CHECK the mapping above (NA rows are excluded). Edit the fcase() rules if a subtype is misassigned.\n")

## ---- 2. Sample up to REF_CELLS_PER_TYPE cells per type, balanced across donors ---
set.seed(SEED)
pick <- ann[!is.na(ref_type)][, {
  # round-robin over patients so that every donor contributes (fixes the h5-order issue of script 17)
  sp <- split(Index, Patient); sp <- lapply(sp, function(v) v[sample.int(length(v))])
  out <- character(); i <- 1
  while (length(out) < REF_CELLS_PER_TYPE && any(lengths(sp) >= i)) {
    out <- c(out, unlist(lapply(sp, function(v) if (length(v) >= i) v[i])))
    i <- i + 1 }
  .(Index = head(out, REF_CELLS_PER_TYPE))
}, by = ref_type]
write_out(pick[, .N, by = ref_type], "60_reference_cells_per_type")
write_out(merge(pick, ann[, .(Index, Patient)], by = "Index")[, .N, by = .(ref_type, Patient)], "60_reference_donor_balance")

## ---- 3. Read only the sampled cells from the raw UMI matrix -----------------------
hdr <- names(fread(GSE132465_COUNTS, nrows = 0))
cols <- c(hdr[1], intersect(pick$Index, hdr))
msg("Reading", length(cols) - 1, "cells from", basename(GSE132465_COUNTS))
mat <- fread(GSE132465_COUNTS, select = cols)
genes <- mat[[1]]; mat <- as.matrix(mat[, -1]); rownames(mat) <- genes
stopifnot(all(abs(mat - round(mat)) < 1e-8))
if (anyDuplicated(genes)) mat <- rowsum(mat, genes)        # sum duplicated symbols
ref_counts <- as(mat, "dgCMatrix")
qc <- colSums(ref_counts) >= 500 & colSums(ref_counts > 0) >= 200     # same fixed QC as script 17
ref_counts <- ref_counts[, qc]
ref_types <- factor(pick$ref_type[match(colnames(ref_counts), pick$Index)])
ref_se <- SummarizedExperiment(assays = list(counts = ref_counts), colData = DataFrame(cell_type = ref_types))
saveRDS(ref_se, file.path(OUT_DIR, "60_reference_GSE132465.rds"))
msg("Reference:", ncol(ref_counts), "cells;", paste(names(table(ref_types)), table(ref_types), collapse = "; "))

## ---- 4. Run RCTD on each exported section and store the weights --------------------
files <- list.files(SPATIAL_SECTION_DIR, pattern = "\\.rds$", full.names = TRUE)
for (f in files) {
  s <- readRDS(f)
  if (!is.null(s$weights$independent)) { msg("skip (done)", s$section); next }
  xy <- as.matrix(s$coords[match(colnames(s$counts), s$coords$spot), c("array_col", "array_row")])
  rownames(xy) <- colnames(s$counts)
  spe <- SpatialExperiment(assays = list(counts = s$counts), spatialCoords = xy)
  rctd <- createRctd(spe, ref_se, cell_type_col = "cell_type")
  out <- runRctd(rctd, rctd_mode = RCTD_MODE, max_cores = RCTD_CORES)
  w <- t(as.matrix(assay(out, "weights")))                  # spots x types
  w <- w / rowSums(w)                                        # normalised proportions
  s$weights$independent <- w
  saveRDS(s, f)
  msg("RCTD (independent reference) done:", s$section, nrow(w), "spots")
}
cat("\nNext: re-run 59_spatial_rerun.R; rows with weights = 'independent' give criterion A with the\n",
    "independent reference. Report pooled effects, patients positive and any label change.\n")
