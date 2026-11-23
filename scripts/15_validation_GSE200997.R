
# 15_validation_GSE200997.R  (protocol v1.0; frozen panel v2; own annotation per section 10b)
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(digest) })
set.seed(20261009)
MIN_CELLS <- 20; NBOOT <- 2000; MIN_DONORS <- 10; MIN_EVAL_COVER <- 0.80
frz_file <- "04_processed_data/frozen_panel/frozen_panel_v2.rds"
if (digest(file = frz_file, algo = "sha256") != "a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb") stop("Frozen panel hash mismatch - STOP")
frz <- readRDS(frz_file); cat("Frozen panel verified\n")
out_dir <- "04_processed_data/GSE200997"; res_dir <- "06_results/validation_GSE200997"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
counts <- readRDS(file.path(out_dir, "GSE200997_counts_sparse.rds"))
cellmap <- readRDS(file.path(out_dir, "GSE200997_cell_assignment.rds"))   # QC-passing cells only
cellmap <- cellmap[condition == "Tumor" & assigned %in% c("Epithelial", "Stromal", "Endothelial", "Myeloid")]
val <- list()
for (comp in unique(cellmap$assigned)) {
  cc <- cellmap[assigned == comp]; ix <- match(cc$cell, colnames(counts)); stopifnot(!anyNA(ix))
  f2 <- factor(cc$patient)
  ind <- sparseMatrix(i = seq_along(ix), j = as.integer(f2), x = 1, dims = c(length(ix), nlevels(f2)), dimnames = list(NULL, levels(f2)))
  pb <- as.matrix(counts[, ix] %*% ind); nc <- table(f2); keep <- names(nc)[nc >= MIN_CELLS]
  if (length(keep)) val[[comp]] <- cpm(pb[, keep, drop = FALSE], log = TRUE, prior.count = 2)
}
saveRDS(val, file.path(out_dir, "GSE200997_validation_pseudobulk_logcpm.rds")); rm(counts); invisible(gc())
score_cohort <- function(L, P) {
  ev <- P$evaluation_genes[gene_symbol %in% rownames(L)]; pn <- P$panel
  E <- t(L[ev$gene_symbol, , drop = FALSE]); X <- t(L[pn$gene_symbol, , drop = FALSE])
  zs <- function(M) { Z <- scale(M); Z[!is.finite(Z)] <- 0; Z }
  data.table(donor = rownames(E),
    y_frozen = rowMeans(sweep(sweep(E, 2, ev$eval_mean, "-"), 2, ev$eval_sd, "/")),
    p_frozen = as.numeric(P$intercept + sweep(sweep(X, 2, pn$pred_mean, "-"), 2, pn$pred_sd, "/") %*% pn$coefficient),
    y_within = rowMeans(zs(E)), p_within = as.numeric(P$intercept + zs(X) %*% pn$coefficient))
}
res <- list(); cover <- list(); preds <- list()
for (pg in names(frz$programs)) {
  P <- frz$programs[[pg]]; L <- val[[P$compartment]]
  nd <- if (is.null(L)) 0 else ncol(L)
  pan <- if (is.null(L)) 0 else sum(P$panel$gene_symbol %in% rownames(L))
  evc <- if (is.null(L)) NA_real_ else mean(P$evaluation_genes$gene_symbol %in% rownames(L))
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
    status <- if (!is.finite(ratio) || ratio >= 1) "Not supported" else if (ci[2, 1] < 1) "Supported" else "Partially supported"
    res[[paste(pg, sc)]] <- data.table(program = pg, scaling = sc, n_donors = n, MAE_ratio = ratio,
      ratio_lo = ci[1, 1], ratio_hi = ci[2, 1], spearman = suppressWarnings(cor(y, p, method = "spearman")),
      rho_lo = ci[1, 2], rho_hi = ci[2, 2], cal_slope = if (sd(p) > 0) unname(coef(lm(y ~ p))[2]) else NA_real_, status = status)
  }
}
res <- rbindlist(res, fill = TRUE); cover <- rbindlist(cover); preds <- rbindlist(preds, fill = TRUE)
num <- names(res)[vapply(res, is.numeric, logical(1))]; res[, (num) := lapply(.SD, round, 3), .SDcols = num]
fwrite(res, file.path(res_dir, "validation_results.csv")); fwrite(cover, file.path(res_dir, "validation_gene_coverage.csv"))
fwrite(preds, file.path(res_dir, "validation_predictions.csv"))
cat("\n=== Coverage ===\n"); print(cover)
cat("\n=== GSE200997 PRIMARY (frozen) ===\n"); print(res[scaling == "frozen", .(program, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, status)])
cat("\n=== GSE200997 secondary (within) ===\n"); print(res[scaling == "within", .(program, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, status)])
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_15.txt"))
 
