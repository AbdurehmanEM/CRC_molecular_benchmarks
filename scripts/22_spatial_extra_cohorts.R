
# 22_spatial_extra_cohorts.R  (protocol v1.5) - extra spatial cohorts: GSE280318 (Visium HD, 56-um pseudo-spots) [+ GSE226997 later]
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(digest)
  library(spacexr); library(SpatialExperiment); library(SummarizedExperiment) })
set.seed(20261013); NB <- 1000; BLK <- 8
frz_file <- "04_processed_data/frozen_panel/frozen_panel_v2.rds"
if (digest(file = frz_file, algo = "sha256") != "a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb") stop("Frozen panel hash mismatch - STOP")
frz <- readRDS(frz_file)
odir <- "04_processed_data/spatial_extra"; dir.create(odir, recursive = TRUE, showWarnings = FALSE)
res_dir <- "06_results/spatial_support_extra"; dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
raw <- "03_original_data/GSE280318"; unz <- "04_processed_data/GSE280318"
SECTIONS <- data.table(
  cohort  = "GSE280318",
  section = c("GSM8594567_P1CRC", "GSM8594568_P2CRC", "GSM8594569_P5CRC", "GSM8594561_P2CRC"),
  patient = c("GSE280318_P1", "GSE280318_P2", "GSE280318_P5", "GSE280318_P2"),
  platform = c("VisiumHD", "VisiumHD", "VisiumHD", "Visium"),
  role = c("main", "main", "main", "sensitivity_same_patient"))
SECTIONS[, h5 := file.path(raw, paste0(section, "_filtered_feature_bc_matrix.h5"))]
SECTIONS[, pos_file := vapply(section, function(x) list.files(unz, paste0("^", x, "_tissue_positions"), full.names = TRUE)[1], "")]
x226 <- "04_processed_data/GSE226997"   # GSE226997 (Ajou, FFPE Visium; primary CRC; treatment status not reported) - protocol v1.6
SECTIONS <- rbind(SECTIONS, data.table(cohort = "GSE226997",
  section = paste0(c("GSM7089855", "GSM7089856", "GSM7089857", "GSM7089858"), "_P", 1:4), patient = paste0("GSE226997_P", 1:4),
  platform = "Visium", role = "main",
  h5 = file.path(x226, paste0("Ajou_Visium_P", 1:4), "filtered_feature_bc_matrix.h5"),
  pos_file = file.path(x226, paste0("Ajou_Visium_P", 1:4), "Spatial", "tissue_positions_list.csv")))
stopifnot(all(file.exists(SECTIONS$h5)), all(file.exists(SECTIONS$pos_file)))
mad_hi <- function(x) { m <- median(x); s <- mad(x); if (!is.finite(s) || s == 0) rep(FALSE, length(x)) else x > m + 3 * s }
load_section <- function(r) {
  cnt <- Seurat::Read10X_h5(r$h5); if (is.list(cnt)) cnt <- cnt[["Gene Expression"]]
  pf <- r$pos_file
  pos <- if (grepl("parquet$", pf)) as.data.table(arrow::read_parquet(pf)) else if (grepl("_list\\.csv$", pf)) fread(pf, header = FALSE) else fread(pf)
  setnames(pos, 1:6, c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  pos <- pos[in_tissue == 1 & barcode %in% colnames(cnt)]; cnt <- cnt[, pos$barcode]
  if (r$platform == "VisiumHD") {
    pos[, ps := paste0("ps_", array_row %/% 7, "_", array_col %/% 7)]
    f <- factor(pos$ps)
    ind <- sparseMatrix(i = seq_len(nrow(pos)), j = as.integer(f), x = 1, dims = c(nrow(pos), nlevels(f)), dimnames = list(NULL, levels(f)))
    cnt <- cnt %*% ind
    pos <- pos[, .(in_tissue = 1L, array_row = array_row[1] %/% 7, array_col = array_col[1] %/% 7,
                   pxl_row = mean(pxl_row), pxl_col = mean(pxl_col), n_bins = .N), by = .(barcode = ps)]
    pos <- pos[match(colnames(cnt), barcode)]
  }
  list(cnt = as(cnt, "CsparseMatrix"), pos = pos)
}
R <- readRDS("04_processed_data/spatial_reference/GSE178341_reference_counts.rds")
ref_se <- SummarizedExperiment(assays = list(counts = R$counts),
  colData = S4Vectors::DataFrame(cell_type = factor(R$cell_types[colnames(R$counts)])))
types <- levels(factor(R$cell_types)); ncores <- max(1, parallel::detectCores() - 1)
for (k in seq_len(nrow(SECTIONS))) {
  r <- SECTIONS[k]; of <- file.path(odir, paste0(r$section, "_rctd.rds"))
  if (file.exists(of)) { message("RCTD done already: ", r$section); next }
  t0 <- Sys.time(); message("=== ", r$section, " (", r$platform, ") started ", format(t0, "%H:%M"))
  S <- load_section(r); cnt <- S$cnt; pos <- S$pos
  mt <- grepl("^MT-", rownames(cnt))
  qc <- data.table(barcode = colnames(cnt), nUMI = Matrix::colSums(cnt), nGene = as.numeric(diff(cnt@p)),
                   pct_mt = if (any(mt)) 100 * Matrix::colSums(cnt[mt, , drop = FALSE]) / pmax(Matrix::colSums(cnt), 1) else 0)
  qc[, pass := nUMI >= 500 & nGene >= 250 & !mad_hi(pct_mt)]
  keep <- qc[pass == TRUE, barcode]; pk <- pos[match(keep, barcode)]
  message("   spots/pseudo-spots: ", nrow(qc), " | pass QC: ", length(keep))
  xy <- cbind(x = pk$pxl_col, y = -pk$pxl_row); rownames(xy) <- keep
  spe <- SpatialExperiment(assays = list(counts = cnt[, keep]), spatialCoords = xy)
  res <- runRctd(createRctd(spe, ref_se, cell_type_col = "cell_type"), rctd_mode = "full", max_cores = ncores)
  W <- as.matrix(assay(res, "weights")); if (!all(rownames(W) %in% keep)) W <- t(W)
  if (is.null(colnames(W)) || !all(colnames(W) %in% types)) colnames(W) <- types
  if (is.null(rownames(W)) || !all(rownames(W) %in% keep)) { stopifnot(nrow(W) == length(keep)); rownames(W) <- keep }
  W <- W / rowSums(W)
  saveRDS(list(section = r$section, patient = r$patient, cohort = r$cohort, platform = r$platform, role = r$role,
               weights = W, qc = qc, pos = pos, cnt = cnt[, keep]), of)
  message("=== ", r$section, " finished in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}
mk <- list(B_Plasma = c("MS4A1","CD79A","JCHAIN","IGKC"), Endothelial = c("PECAM1","VWF","PLVAP"),
  Fibroblast = c("COL1A1","DCN","LUM"), Malignant_epithelial = c("EPCAM","CEACAM5","OLFM4"),
  Mast = c("TPSAB1","CPA3"), Myeloid = c("C1QA","CD68","LYZ"), Normal_epithelial = c("CA1","SLC26A3","AQP8","GUCA2A"),
  Pericyte_SMC = c("RGS5","ACTA2","DES"), T_NK = c("CD3E","NKG7","CD2"))
chk <- rbindlist(lapply(list.files(odir, "_rctd.rds$", full.names = TRUE), function(p) { x <- readRDS(p)
  lc <- log1p(t(t(x$cnt) / Matrix::colSums(x$cnt)) * 1e4)
  Sm <- sapply(mk, function(g) { g <- intersect(g, rownames(lc)); if (length(g)) Matrix::colMeans(lc[g, , drop = FALSE]) else rep(NA_real_, ncol(lc)) })
  C <- cor(Sm, x$weights[, types], method = "spearman", use = "pairwise.complete.obs")
  data.table(section = x$section, type = types, own_rho = round(diag(C), 2), own_is_max = apply(C, 1, which.max) == seq_along(types),
             mean_prop = round(colMeans(x$weights[, types]), 3)) }))
fwrite(chk, file.path(res_dir, "RCTD_marker_check_extra.csv"))
comp2type <- c(Epithelial = "Malignant_epithelial", Stromal = "Fibroblast", Myeloid = "Myeloid", Endothelial = "Endothelial")
zs <- function(M) { Z <- scale(M); Z[!is.finite(Z)] <- 0; Z }
fitb <- function(y, p, Wm) { X <- cbind(1, p, Wm); unname(qr.coef(qr(X), y)[2]) }
A <- list(); Bq <- list(); cover <- list()
for (p in list.files(odir, "_rctd.rds$", full.names = TRUE)) {
  x <- readRDS(p); W <- x$weights; bc <- rownames(W); cnt <- x$cnt[, bc]
  L <- edgeR::cpm(as.matrix(cnt[intersect(unique(unlist(lapply(frz$programs, function(P) c(P$panel$gene_symbol, P$evaluation_genes$gene_symbol)))), rownames(cnt)), , drop = FALSE]),
           lib.size = Matrix::colSums(cnt), log = TRUE, prior.count = 2)
  pos <- x$pos[match(bc, barcode)]; blk <- paste(pos$array_row %/% BLK, pos$array_col %/% BLK)
  Wm <- W[, setdiff(colnames(W), "Mast"), drop = FALSE]
  for (pg in names(frz$programs)) {
    P <- frz$programs[[pg]]; ev <- intersect(P$evaluation_genes$gene_symbol, rownames(L))
    evc <- length(ev) / nrow(P$evaluation_genes); pan_ok <- all(P$panel$gene_symbol %in% rownames(L))
    cover[[paste(x$section, pg)]] <- data.table(section = x$section, program = pg, panel_present = sum(P$panel$gene_symbol %in% rownames(L)),
      panel_genes = nrow(P$panel), missing_panel = paste(setdiff(P$panel$gene_symbol, rownames(L)), collapse = ";"), eval_coverage = round(evc, 3))
    if (!pan_ok || evc < 0.80) next
    y <- rowMeans(zs(t(L[ev, , drop = FALSE]))); pz <- as.numeric(P$intercept + zs(t(L[P$panel$gene_symbol, , drop = FALSE])) %*% P$panel$coefficient)
    ys <- as.numeric(scale(y)); ps <- if (sd(pz) > 0) as.numeric(scale(pz)) else rep(0, length(pz))
    est <- fitb(ys, ps, Wm); ub <- unique(blk); bi <- split(seq_along(blk), blk)
    bs <- replicate(NB, { i <- unlist(bi[sample(ub, length(ub), TRUE)], use.names = FALSE); fitb(ys[i], ps[i], Wm[i, , drop = FALSE]) })
    cp <- W[, comp2type[[P$compartment]]]; hi <- cp >= quantile(cp, 0.75); lo <- cp <= quantile(cp, 0.25)
    A[[paste(x$section, pg)]] <- data.table(cohort = x$cohort, section = x$section, patient = x$patient, platform = x$platform, role = x$role,
      program = pg, n_spots = length(ys), effect = est, se = sd(bs, na.rm = TRUE), lo = quantile(bs, 0.025, na.rm = TRUE), hi = quantile(bs, 0.975, na.rm = TRUE))
    Bq[[paste(x$section, pg)]] <- data.table(cohort = x$cohort, section = x$section, patient = x$patient, role = x$role, program = pg,
      rule = "composition_quartile", n_hi = sum(hi), n_lo = sum(lo), diff = mean(y[hi]) - mean(y[lo]))
  }
  message("scored ", x$section)
}
A <- rbindlist(A); Bq <- rbindlist(Bq); cover <- rbindlist(cover)
rnd <- function(X) { nn <- names(X)[vapply(X, is.numeric, logical(1))]; X[, (nn) := lapply(.SD, round, 3), .SDcols = nn][] }
fwrite(rnd(A), file.path(res_dir, "criterionA_by_section_extra.csv")); fwrite(rnd(Bq), file.path(res_dir, "criterionB_by_section_extra.csv"))
fwrite(cover, file.path(res_dir, "spatial_gene_coverage_extra.csv"))
cat("\n=== RCTD marker check (own-marker rho) ===\n"); print(dcast(chk, type ~ section, value.var = "own_rho"))
cat("\n=== Mean estimated proportions ===\n"); print(dcast(chk, type ~ section, value.var = "mean_prop"))
cat("\n=== Coverage failures ===\n"); print(cover[panel_present < panel_genes | eval_coverage < 0.8])
cat("\n=== Criterion A by section (standardised panel coefficient [95% CI]) ===\n")
print(dcast(A[, .(program, section, v = sprintf("%.3f [%.3f, %.3f]", effect, lo, hi))], program ~ section, value.var = "v"))
cat("\n=== Criterion B by section (composition quartiles) ===\n"); print(dcast(rnd(Bq), program ~ section, value.var = "diff"))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_22.txt"))
 
