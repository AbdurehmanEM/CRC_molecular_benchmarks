
# 18b_RCTD_label_check.R - restore spot/cell-type labels on RCTD weights and verify them with marker genes
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
dst <- "04_processed_data/Valdeolivas_Visium/unzipped"; out <- "04_processed_data/Valdeolivas_Visium/rctd"
R <- readRDS("04_processed_data/spatial_reference/GSE178341_reference_counts.rds")
types <- levels(factor(R$cell_types)); cat("Assumed column order:", paste(types, collapse = ", "), "\n")
mk <- list(B_Plasma = c("MS4A1","CD79A","JCHAIN","IGKC"), Endothelial = c("PECAM1","VWF","PLVAP"),
  Fibroblast = c("COL1A1","DCN","LUM"), Malignant_epithelial = c("EPCAM","CEACAM5","OLFM4"),
  Mast = c("TPSAB1","CPA3"), Myeloid = c("C1QA","CD68","LYZ"), Normal_epithelial = c("CA1","SLC26A3","AQP8","GUCA2A"),
  Pericyte_SMC = c("RGS5","ACTA2","DES"), T_NK = c("CD3E","NKG7","CD2"))
C <- list(); chk <- list()
for (p in list.files(out, "_rctd.rds$", full.names = TRUE)) {
  x <- readRDS(p); W <- x$weights; keep <- x$qc[pass == TRUE, barcode]
  chk[[x$section]] <- data.table(section = x$section, rows = nrow(W), qc_pass = length(keep), cols = ncol(W))
  if (nrow(W) != length(keep) || ncol(W) != length(types)) next
  cnt <- Seurat::Read10X_h5(file.path(dst, x$section, "filtered_feature_bc_matrix.h5"))[, keep]
  lc <- log1p(t(t(cnt) / Matrix::colSums(cnt)) * 1e4)
  S <- sapply(mk, function(g) Matrix::colMeans(lc[intersect(g, rownames(lc)), , drop = FALSE]))
  C[[x$section]] <- cor(S, W, method = "spearman")
}
chk <- rbindlist(chk); cat("\nDimension check:\n"); print(chk)
M <- Reduce(`+`, C) / length(C); colnames(M) <- paste0("col", seq_along(types), "_", types)
cat("\n=== Mean Spearman: marker score (rows) vs RCTD column (cols), across", length(C), "sections ===\n"); print(round(M, 2))
best <- apply(M, 1, which.max); ok <- all(best == seq_along(types)) && all(chk$rows == chk$qc_pass)
cat("\nEach marker set matches its own column:", ok, "\n")
if (ok) {
  for (p in list.files(out, "_rctd.rds$", full.names = TRUE)) { x <- readRDS(p)
    dimnames(x$weights) <- list(x$qc[pass == TRUE, barcode], types); x$labels_verified <- TRUE; saveRDS(x, p) }
  S <- rbindlist(lapply(list.files(out, "_rctd.rds$", full.names = TRUE), function(p) { x <- readRDS(p)
    data.table(section = x$section, patient = x$patient, spots = nrow(x$weights), t(round(colMeans(x$weights), 3))) }))
  fwrite(S, "06_results/spatial_RCTD_section_summary.csv"); fwrite(data.table(marker_set = rownames(M), round(M, 3)), "06_results/spatial_RCTD_marker_check.csv")
  cat("\n=== Mean estimated proportion per section ===\n"); print(S)
} else cat("Labels NOT saved - send me this output\n")
 

dir.create("06_results/sessionInfo", recursive = TRUE, showWarnings = FALSE); writeLines(capture.output(sessionInfo()), "06_results/sessionInfo/sessionInfo_18b.txt")
