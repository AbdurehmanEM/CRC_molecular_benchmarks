# 13_validation_GSE166555.R  (protocol v0.9; frozen panel v2) - Uhlitz 2021, 12 patients, tumour samples
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(digest) })
set.seed(20261008)
MIN_CELLS <- 20; NBOOT <- 2000; MIN_DONORS <- 10; MIN_EVAL_COVER <- 0.80
frz_file <- "04_processed_data/frozen_panel/frozen_panel_v2.rds"
if (digest(file = frz_file, algo = "sha256") != "a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb") stop("Frozen panel hash mismatch - STOP")
frz <- readRDS(frz_file); cat("Frozen panel verified\n")
# 19 Oct 2026: archive members are extracted to 04_processed_data (03_original_data
# must contain only the downloaded files), and cell barcodes are checked for
# uniqueness after the per-sample matrices are column-bound - duplicated barcodes
# would make match() assign the first sample's metadata to cells of later samples.
in_dir <- "03_original_data/GSE166555"; raw_dir <- "04_processed_data/GSE166555/raw"; res_dir <- "06_results/validation_GSE166555"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE); dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
meta <- fread(file.path(in_dir, "GSE166555_meta_data.tsv.gz"))
meta[, val_compartment := fcase(
  minor_cell_type %in% c("TC1", "TC2", "TC3", "TC4"), "Epithelial",
  minor_cell_type %in% c("CAFs", "CB FBs", "UC FBs", "CCL FBs", "MyoFBs 1", "MyoFBs 2"), "Stromal",
  minor_cell_type %in% c("Monocyte 1", "Monocyte 2"), "Myeloid",
  minor_cell_type == "Endothelial", "Endothelial", default = NA_character_)]
tum_files <- grep("t[0-9]*\\.tsv\\.gz$", untar(file.path(in_dir, "GSE166555_RAW.tar"), list = TRUE), value = TRUE)
if (!all(file.exists(file.path(raw_dir, tum_files)))) untar(file.path(in_dir, "GSE166555_RAW.tar"), files = tum_files, exdir = raw_dir)
cat("Tumour sample files:", length(tum_files), "\n")
mats <- list()
for (fl in tum_files) {
  x <- fread(file.path(raw_dir, fl)); g <- x[[1]]
  m <- as(as.matrix(x[, -1]), "CsparseMatrix"); rownames(m) <- g
  mats[[fl]] <- m; message("read ", fl, ": ", ncol(m), " cells")
}
genes <- Reduce(intersect, lapply(mats, rownames))
counts <- do.call(cbind, lapply(mats, function(m) m[genes, , drop = FALSE])); rm(mats)
stopifnot(all(abs(counts@x - round(counts@x)) < 1e-8))
if (anyDuplicated(colnames(counts))) stop(sum(duplicated(colnames(counts))), " duplicated cell barcodes across sample files - metadata matching would be wrong")
cm <- meta[match(colnames(counts), cell)]
cat("Cells:", ncol(counts), "| genes:", nrow(counts), "| matched:", sum(!is.na(cm$case_id)), "\n")
mt <- grepl("^MT-", genes)
cm[, nCount := Matrix::colSums(counts)]; cm[, nFeature := as.numeric(diff(counts@p))]
cm[, pct_mt := 100 * Matrix::colSums(counts[mt, , drop = FALSE]) / pmax(nCount, 1)]
mad_flag <- function(x, n = 3, side = c("both", "high")) { side <- match.arg(side); m <- median(x); s <- mad(x)
  if (!is.finite(s) || s == 0) return(rep(FALSE, length(x)))
  if (side == "high") x > m + n * s else x < m - n * s | x > m + n * s }
cm[, qc_flag := mad_flag(log10(nCount + 1)) | mad_flag(log10(nFeature + 1)) | mad_flag(pct_mt, side = "high"), by = sample_id]
use <- !is.na(cm$val_compartment) & !cm$qc_flag & cm$sample_origin == "Tumor"
cat("\nTumour cells per patient x compartment (post-QC; t1/t2 pooled per patient):\n")
print(dcast(cm[use, .N, by = .(case_id, val_compartment)], case_id ~ val_compartment, value.var = "N", fill = 0))
val <- list()
for (comp in unique(na.omit(cm$val_compartment))) {
  ix <- which(use & cm$val_compartment == comp); f2 <- factor(cm$case_id[ix])
  ind <- sparseMatrix(i = seq_along(ix), j = as.integer(f2), x = 1, dims = c(length(ix), nlevels(f2)), dimnames = list(NULL, levels(f2)))
  pb <- as.matrix(counts[, ix] %*% ind); nc <- table(f2); keep <- names(nc)[nc >= MIN_CELLS]
  if (length(keep)) val[[comp]] <- cpm(pb[, keep, drop = FALSE], log = TRUE, prior.count = 2)
}
saveRDS(val, "04_processed_data/GSE166555_validation_pseudobulk_logcpm.rds")
 
score_cohort <- function(L, P) {
  ev <- P$evaluation_genes[gene_symbol %in% rownames(L)]; pn <- P$panel
  E <- t(L[ev$gene_symbol, , drop = FALSE]); X <- t(L[pn$gene_symbol, , drop = FALSE])
  data.table(donor = rownames(E),
    y_frozen = rowMeans(sweep(sweep(E, 2, ev$eval_mean, "-"), 2, ev$eval_sd, "/")),
    p_frozen = as.numeric(P$intercept + sweep(sweep(X, 2, pn$pred_mean, "-"), 2, pn$pred_sd, "/") %*% pn$coefficient),
    y_within = if (nrow(E) > 2) rowMeans(scale(E), na.rm = TRUE) else NA_real_,
    p_within = if (nrow(E) > 2) { Z <- scale(X); Z[!is.finite(Z)] <- 0; as.numeric(P$intercept + Z %*% pn$coefficient) } else NA_real_)
}
res <- list(); cover <- list(); preds <- list()
for (pg in names(frz$programs)) {
  P <- frz$programs[[pg]]; L <- val[[P$compartment]]
  nd <- if (is.null(L)) 0 else ncol(L)
  pan <- if (is.null(L)) 0 else sum(P$panel$gene_symbol %in% rownames(L))
  evc <- if (is.null(L)) NA else mean(P$evaluation_genes$gene_symbol %in% rownames(L))
  cover[[pg]] <- data.table(program = pg, n_donors = nd, panel_present = pan, panel_genes = nrow(P$panel),
    missing_panel = if (is.null(L)) "" else paste(setdiff(P$panel$gene_symbol, rownames(L)), collapse = ";"),
    eval_coverage = round(evc, 3))
  if (is.null(L) || pan < nrow(P$panel) || evc < MIN_EVAL_COVER || nd < MIN_DONORS) {
    res[[pg]] <- data.table(program = pg, scaling = c("frozen", "within"), n_donors = nd, status = "Not evaluable"); next }
  d <- score_cohort(L, P); preds[[pg]] <- cbind(program = pg, d)
  for (sc in c("frozen", "within")) {
    y <- d[[paste0("y_", sc)]]; p <- d[[paste0("p_", sc)]]; n <- length(y)
    ratio <- mean(abs(y - p)) / mean(abs(y))
    bs <- replicate(NBOOT, { i <- sample.int(n, n, TRUE)
      c(mean(abs(y[i] - p[i])) / mean(abs(y[i])), suppressWarnings(cor(y[i], p[i], method = "spearman"))) })
    ci <- apply(bs, 1, quantile, c(0.025, 0.975), na.rm = TRUE)
    status <- if (!is.finite(ratio)) "Not evaluable" else if (ratio >= 1) "Not supported" else if (ci[2, 1] < 1) "Supported" else "Partially supported"
    res[[paste(pg, sc)]] <- data.table(program = pg, scaling = sc, n_donors = n, MAE_ratio = ratio,
      ratio_lo = ci[1, 1], ratio_hi = ci[2, 1], spearman = suppressWarnings(cor(y, p, method = "spearman")),
      rho_lo = ci[1, 2], rho_hi = ci[2, 2], cal_slope = unname(coef(lm(y ~ p))[2]), status = status)
  }
}
res <- rbindlist(res, fill = TRUE); cover <- rbindlist(cover); preds <- rbindlist(preds, fill = TRUE)
num <- names(res)[vapply(res, is.numeric, logical(1))]; res[, (num) := lapply(.SD, round, 3), .SDcols = num]
fwrite(res, file.path(res_dir, "validation_results.csv")); fwrite(cover, file.path(res_dir, "validation_gene_coverage.csv"))
fwrite(preds, file.path(res_dir, "validation_predictions.csv"))
cat("\n=== Coverage ===\n"); print(cover)
cat("\n=== GSE166555 PRIMARY (frozen) ===\n"); print(res[scaling == "frozen", .(program, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, status)])
cat("\n=== GSE166555 secondary (within) ===\n"); print(res[scaling == "within", .(program, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, status)])
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_13.txt"))
 
