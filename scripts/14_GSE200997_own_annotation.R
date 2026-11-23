
# 14_GSE200997_own_annotation.R  (protocol v1.0, section 10b) - own annotation of GSE200997
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
set.seed(20261009)
in_dir <- "03_original_data/GSE200997"; out_dir <- "04_processed_data/GSE200997"
res_dir <- "06_results/GSE200997_annotation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
cnt_rds <- file.path(out_dir, "GSE200997_counts_sparse.rds")
if (!file.exists(cnt_rds)) {
  con <- gzfile(file.path(in_dir, "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz"), "r")
  cells <- gsub("\"", "", strsplit(readLines(con, 1), ",", fixed = TRUE)[[1]])[-1]
  blocks <- list(); genes <- character(0)
  repeat {
    lines <- readLines(con, n = 1000); if (!length(lines)) break
    dt <- fread(text = lines, header = FALSE, sep = ",")
    genes <- c(genes, dt[[1]]); blocks[[length(blocks) + 1]] <- as(as.matrix(dt[, -1]), "CsparseMatrix")
    message(length(genes), " genes read")
    if (length(lines) < 1000) break
  }
  close(con)
  counts <- do.call(rbind, blocks); rm(blocks); invisible(gc())
  dimnames(counts) <- list(make.unique(genes), cells)
  stopifnot(all(abs(counts@x - round(counts@x)) < 1e-8))
  saveRDS(counts, cnt_rds)
} else counts <- readRDS(cnt_rds)
cat("Counts:", nrow(counts), "genes x", ncol(counts), "cells | duplicated gene names made unique:", sum(duplicated(sub("\\.[0-9]+$", "", rownames(counts)))), "\n")
 

# ---- Part 2: metadata, QC, clustering, pre-specified marker assignment (no program scoring) ----
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
ann <- fread(file.path(in_dir, "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz"), header = TRUE)
setnames(ann, 1, "cell")
cat("Annotation columns:", paste(names(ann), collapse = ", "), "\n")
scol <- grep("sample", names(ann), ignore.case = TRUE, value = TRUE)[1]
ccol <- grep("condition", names(ann), ignore.case = TRUE, value = TRUE)[1]
meta <- ann[match(colnames(counts), cell)]; stopifnot(!anyNA(meta$cell))
meta[, sample_id := get(scol)]; meta[, condition := get(ccol)]
meta[, patient := sub("^[^_]+_", "", sample_id)]
cat("\nCells per patient x condition:\n")
print(dcast(meta[, .N, by = .(patient, condition)], patient ~ condition, value.var = "N", fill = 0))
mt <- grepl("^MT-", rownames(counts)); cat("Mito genes:", sum(mt), "\n")
meta[, nCount := Matrix::colSums(counts)]; meta[, nFeature := as.numeric(diff(counts@p))]
meta[, pct_mt := 100 * Matrix::colSums(counts[mt, , drop = FALSE]) / pmax(nCount, 1)]
mad_flag <- function(x, n = 3, side = c("both", "high")) { side <- match.arg(side); m <- median(x); s <- mad(x)
  if (!is.finite(s) || s == 0) return(rep(FALSE, length(x)))
  if (side == "high") x > m + n * s else x < m - n * s | x > m + n * s }
meta[, qc_flag := mad_flag(log10(nCount + 1)) | mad_flag(log10(nFeature + 1)) | mad_flag(pct_mt, side = "high"), by = sample_id]
cat("QC removed:", sum(meta$qc_flag), "of", nrow(meta), "cells\n")
keep <- !meta$qc_flag
so <- CreateSeuratObject(counts[, keep], meta.data = data.frame(meta[keep], row.names = meta$cell[keep]))
rm(counts); invisible(gc())
so <- NormalizeData(so, verbose = FALSE); so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE); so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
so <- FindClusters(so, resolution = 0.5, random.seed = 20261009, verbose = FALSE)
markers <- list(
  Epithelial = c("EPCAM","KRT8","KRT18","KRT19","CDH1"),
  Stromal = c("COL1A1","COL1A2","DCN","LUM","PDGFRA"),
  Endothelial = c("PECAM1","VWF","CDH5","PLVAP","CLDN5"),
  Myeloid = c("CD14","CD68","LYZ","C1QA","CSF1R"),
  ctx_LymphaticEC = c("PROX1","LYVE1","CCL21"), ctx_Pericyte = c("RGS5","NOTCH3"),
  ctx_DC = c("CD1C","CLEC10A","LAMP3"), ctx_T_NK = c("CD3E","NKG7"),
  ctx_B_Plasma = c("MS4A1","CD79A","JCHAIN"), ctx_Mast = c("TPSAB1","CPA3"))
D <- LayerData(so, assay = "RNA", layer = "data"); cl <- so$seurat_clusters
allm <- unique(unlist(markers)); cat("Marker genes missing:", paste(setdiff(allm, rownames(D)), collapse = ", "), "\n")
ind <- sparseMatrix(i = seq_along(cl), j = as.integer(cl), x = 1, dims = c(length(cl), nlevels(cl)))
ind <- ind %*% Diagonal(x = 1 / as.numeric(table(cl)))
CM <- as.matrix(D[intersect(allm, rownames(D)), ] %*% ind); colnames(CM) <- levels(cl)
S <- sapply(markers, function(g) colMeans(CM[intersect(g, rownames(CM)), , drop = FALSE]))
best <- colnames(S)[max.col(S, ties.method = "first")]
margin <- apply(S, 1, function(v) { s <- sort(v, decreasing = TRUE); s[1] - s[2] })
ctab <- data.table(cluster = levels(cl), n_cells = as.integer(table(cl)), assigned = best, margin = round(margin, 2), round(S, 2))
fwrite(ctab, file.path(res_dir, "cluster_marker_assignment.csv"))
so$assigned <- best[as.integer(cl)]
p <- DotPlot(so, features = allm[allm %in% rownames(D)], group.by = "seurat_clusters") + RotatedAxis()
ggsave(file.path(res_dir, "marker_dotplot.pdf"), p, width = 16, height = 8)
cellmap <- data.table(cell = colnames(so), sample_id = so$sample_id, patient = so$patient,
  condition = so$condition, cluster = as.character(cl), assigned = so$assigned)
saveRDS(cellmap, file.path(out_dir, "GSE200997_cell_assignment.rds"))
cat("\n=== Cluster assignment ===\n"); print(ctab[, .(cluster, n_cells, assigned, margin)], nrows = 100)
cat("\n=== Tumour-sample cells per patient x modelled compartment ===\n")
print(dcast(cellmap[grepl("tum", condition, ignore.case = TRUE) & assigned %in% c("Epithelial","Stromal","Endothelial","Myeloid"), .N, by = .(patient, assigned)],
  patient ~ assigned, value.var = "N", fill = 0))
 

dir.create("06_results/sessionInfo", recursive = TRUE, showWarnings = FALSE); writeLines(capture.output(sessionInfo()), "06_results/sessionInfo/sessionInfo_14.txt")
