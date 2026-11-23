
# 18_spatial_RCTD_deconvolution.R  (protocol v1.2, 10d) - spot QC + RCTD full mode (spacexr 1.4 Bioconductor API)
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(spacexr)
  library(SpatialExperiment); library(SummarizedExperiment) })
set.seed(20261012)
dst <- "04_processed_data/Valdeolivas_Visium/unzipped"; out <- "04_processed_data/Valdeolivas_Visium/rctd"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
R <- readRDS("04_processed_data/spatial_reference/GSE178341_reference_counts.rds")
ref_se <- SummarizedExperiment(assays = list(counts = R$counts),
  colData = S4Vectors::DataFrame(cell_type = factor(R$cell_types[colnames(R$counts)])))
secs <- grep("^SN", list.dirs(dst, recursive = FALSE, full.names = FALSE), value = TRUE)
ncores <- max(1, parallel::detectCores() - 1); cat("Sections:", length(secs), "| cores:", ncores, "\n")
mad_hi <- function(x) { m <- median(x); s <- mad(x); if (!is.finite(s) || s == 0) rep(FALSE, length(x)) else x > m + 3 * s }
for (s in secs) {
  of <- file.path(out, paste0(s, "_rctd.rds")); if (file.exists(of)) { message("Done already: ", s); next }
  t0 <- Sys.time(); message("=== ", s, " started ", format(t0, "%H:%M"))
  cnt <- Seurat::Read10X_h5(file.path(dst, s, "filtered_feature_bc_matrix.h5"))
  pf <- file.path(dst, s, "spatial", "tissue_positions_list.csv")
  pos <- fread(pf, header = FALSE); if (!is.numeric(pos[[2]])) pos <- fread(pf, header = TRUE)
  setnames(pos, 1:6, c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  pos <- pos[in_tissue == 1 & barcode %in% colnames(cnt)]; cnt <- cnt[, pos$barcode]
  mt <- grepl("^MT-", rownames(cnt))
  qc <- data.table(barcode = pos$barcode, nUMI = Matrix::colSums(cnt), nGene = as.numeric(diff(cnt@p)),
                   pct_mt = 100 * Matrix::colSums(cnt[mt, , drop = FALSE]) / pmax(Matrix::colSums(cnt), 1))
  qc[, pass := nUMI >= 500 & nGene >= 250 & !mad_hi(pct_mt)]
  keep <- qc[pass == TRUE, barcode]; pk <- pos[match(keep, barcode)]
  xy <- cbind(x = pk$pxl_col, y = -pk$pxl_row); rownames(xy) <- keep
  spe <- SpatialExperiment(assays = list(counts = cnt[, keep]), spatialCoords = xy)
  rd <- createRctd(spe, ref_se, cell_type_col = "cell_type")
  res <- runRctd(rd, rctd_mode = "full", max_cores = ncores)
  W <- as.matrix(assay(res, "weights")); if (!all(rownames(W) %in% keep)) W <- t(W)
  W <- W / rowSums(W)
  saveRDS(list(section = s, patient = sub("^SN[0-9]+_(A[0-9]+)_.*$", "\\1", s), weights = W, qc = qc, pos = pos), of)
  message("=== ", s, " finished: ", nrow(W), " of ", nrow(qc), " spots deconvolved; ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}
S <- rbindlist(lapply(list.files(out, "_rctd.rds$", full.names = TRUE), function(p) { x <- readRDS(p)
  data.table(section = x$section, patient = x$patient, spots_total = nrow(x$qc), spots_kept = nrow(x$weights),
             as.list(round(colMeans(x$weights), 3))) }), fill = TRUE)
fwrite(S, "06_results/spatial_RCTD_section_summary.csv")
cat("\n=== Mean estimated proportion per section ===\n"); print(S)
 

dir.create("06_results/sessionInfo", recursive = TRUE, showWarnings = FALSE); writeLines(capture.output(sessionInfo()), "06_results/sessionInfo/sessionInfo_18.txt")
