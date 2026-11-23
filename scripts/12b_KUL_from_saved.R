# 12b_KUL_from_saved.R
# Protocol v0.9: frozen panel v2 applied to GSE144735 (KUL3, Belgium; tumour core only).
# 19 Oct 2026: builds GSE144735_counts_sparse.rds itself when absent (loader copied
# verbatim from 12_validation_GSE144735_pooled.R), so stage 12 of 00_run_all.R works
# on a clean machine. Results are unchanged: the loader only reads the raw matrix.
# KUL has < 10 donors -> reported descriptively; pooled SMC + KUL = pre-specified sensitivity analysis.
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(digest) })
set.seed(20261008)
MIN_CELLS <- 20; NBOOT <- 2000; MIN_DONORS <- 10; MIN_EVAL_COVER <- 0.80
frz_file <- "04_processed_data/frozen_panel/frozen_panel_v2.rds"
EXPECTED_SHA <- "a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb"
if (digest(file = frz_file, algo = "sha256") != EXPECTED_SHA) stop("Frozen panel hash mismatch - STOP")
frz <- readRDS(frz_file); cat("Frozen panel verified\n")
in_dir <- "03_original_data/GSE144735"; res_dir <- "06_results/validation_GSE144735"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("04_processed_data/GSE144735", recursive = TRUE, showWarnings = FALSE)
annot <- fread(file.path(in_dir, "GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz"))
cnt_rds <- "04_processed_data/GSE144735/GSE144735_counts_sparse.rds"
if (!file.exists(cnt_rds)) {
  con <- gzfile(file.path(in_dir, "GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz"), "r")
  header <- gsub("\"", "", strsplit(readLines(con, 1), "\t", fixed = TRUE)[[1]])
  first <- readLines(con, 1); nf <- length(strsplit(first, "\t", fixed = TRUE)[[1]])
  cells <- if (length(header) == nf) header[-1] else header
  blocks <- list(); genes <- character(0); pending <- first
  repeat {
    lines <- c(pending, readLines(con, n = 2000 - length(pending))); pending <- character(0)
    if (!length(lines)) break
    dt <- fread(text = lines, header = FALSE, sep = "\t")
    genes <- c(genes, dt[[1]]); blocks[[length(blocks) + 1]] <- as(as.matrix(dt[, -1]), "CsparseMatrix"); message("KUL rows read: ", length(genes))
    if (length(lines) < 2000) break
  }
  close(con)
  counts <- do.call(rbind, blocks); dimnames(counts) <- list(genes, cells); rm(blocks)
  stopifnot(all(abs(counts@x - round(counts@x)) < 1e-8))
  saveRDS(counts, cnt_rds)
}
counts <- readRDS(cnt_rds)
genes <- rownames(counts); cells <- colnames(counts)
stopifnot(!anyDuplicated(cells))
meta <- annot[match(cells, Index)]
cat("Cells:", ncol(counts), "| genes:", nrow(counts), "| matched to annotation:", sum(!is.na(meta$Patient)), "\n")
mt <- grepl("^MT-", genes)
meta[, nCount := Matrix::colSums(counts)]; meta[, nFeature := Matrix::colSums(counts > 0)]
meta[, pct_mt := 100 * Matrix::colSums(counts[mt, , drop = FALSE]) / nCount]
mad_flag <- function(x, n = 3, side = c("both", "high")) { side <- match.arg(side); m <- median(x); s <- mad(x)
  if (!is.finite(s) || s == 0) return(rep(FALSE, length(x)))
  if (side == "high") x > m + n * s else x < m - n * s | x > m + n * s }
meta[, qc_flag := mad_flag(log10(nCount + 1)) | mad_flag(log10(nFeature + 1)) | mad_flag(pct_mt, side = "high"), by = Sample]
meta[, val_compartment := fcase(
  Cell_subtype %in% c("Tip-like ECs", "Stalk-like ECs", "Proliferative ECs"), "Endothelial",
  Cell_subtype %in% c("Stromal 1", "Stromal 2", "Stromal 3", "Myofibroblasts"), "Stromal",
  Cell_type == "Myeloids" & !Cell_subtype %in% c("cDC", "Unknown"), "Myeloid",
  Cell_type == "Epithelial cells" & grepl("CMS", Cell_subtype), "Epithelial",
  default = NA_character_)]
use <- meta$Class == "Tumor" & !meta$qc_flag & !is.na(meta$val_compartment)
 
cat("\nKUL tumour-core cells per patient x compartment (post-QC):\n")
print(dcast(meta[use, .N, by = .(Patient, val_compartment)], Patient ~ val_compartment, value.var = "N", fill = 0))
kul <- list()
for (comp in unique(na.omit(meta$val_compartment))) {
  ix <- which(use & meta$val_compartment == comp); if (!length(ix)) next
  f2 <- factor(meta$Patient[ix])
  ind <- sparseMatrix(i = seq_along(ix), j = as.integer(f2), x = 1, dims = c(length(ix), nlevels(f2)), dimnames = list(NULL, levels(f2)))
  pb <- as.matrix(counts[, ix] %*% ind); nc <- table(f2); keep <- names(nc)[nc >= MIN_CELLS]
  if (length(keep)) kul[[comp]] <- cpm(pb[, keep, drop = FALSE], log = TRUE, prior.count = 2)
}
saveRDS(kul, "04_processed_data/GSE144735/validation_pseudobulk_logcpm.rds")
smc <- readRDS("04_processed_data/GSE132465/validation_pseudobulk_logcpm.rds")$logcpm
score_cohort <- function(L, P) {
  ev <- P$evaluation_genes[gene_symbol %in% rownames(L)]; pn <- P$panel
  E <- t(L[ev$gene_symbol, , drop = FALSE]); X <- t(L[pn$gene_symbol, , drop = FALSE])
  data.table(donor = rownames(E),
    y_frozen = rowMeans(sweep(sweep(E, 2, ev$eval_mean, "-"), 2, ev$eval_sd, "/")),
    p_frozen = as.numeric(P$intercept + sweep(sweep(X, 2, pn$pred_mean, "-"), 2, pn$pred_sd, "/") %*% pn$coefficient),
    y_within = if (nrow(E) > 2) rowMeans(scale(E)) else NA_real_,
    p_within = if (nrow(E) > 2) as.numeric(P$intercept + scale(X) %*% pn$coefficient) else NA_real_)
}
evaluate <- function(d, label, pg) {
  out <- list()
  for (sc in c("frozen", "within")) {
    y <- d[[paste0("y_", sc)]]; p <- d[[paste0("p_", sc)]]; ok <- is.finite(y) & is.finite(p)
    y <- y[ok]; p <- p[ok]; coh <- d$cohort[ok]; n <- length(y)
    if (n < 3) { out[[sc]] <- data.table(analysis = label, program = pg, scaling = sc, n_donors = n, status = "Not evaluable"); next }
    ratio <- mean(abs(y - p)) / mean(abs(y))
    bs <- replicate(NBOOT, { i <- unlist(lapply(split(seq_len(n), coh), function(z) z[sample.int(length(z), length(z), TRUE)]))
      c(mean(abs(y[i] - p[i])) / mean(abs(y[i])), suppressWarnings(cor(y[i], p[i], method = "spearman"))) })
    ci <- apply(bs, 1, quantile, c(0.025, 0.975), na.rm = TRUE)
    status <- if (n < MIN_DONORS) "Descriptive only (<10 donors)" else if (ratio >= 1) "Not supported" else if (ci[2, 1] < 1) "Supported" else "Partially supported"
    out[[sc]] <- data.table(analysis = label, program = pg, scaling = sc, n_donors = n,
      MAE_ratio = ratio, ratio_lo = ci[1, 1], ratio_hi = ci[2, 1],
      spearman = suppressWarnings(cor(y, p, method = "spearman")), rho_lo = ci[1, 2], rho_hi = ci[2, 2],
      cal_slope = if (sd(p) > 0) unname(coef(lm(y ~ p))[2]) else NA_real_, status = status)
  }
  rbindlist(out, fill = TRUE)
}
 
res <- list(); cover <- list()
for (pg in names(frz$programs)) {
  P <- frz$programs[[pg]]; comp <- P$compartment
  Lk <- kul[[comp]]; Ls <- smc[[comp]]
  ok_k <- !is.null(Lk) && all(P$panel$gene_symbol %in% rownames(Lk)) &&
          mean(P$evaluation_genes$gene_symbol %in% rownames(Lk)) >= MIN_EVAL_COVER
  cover[[pg]] <- data.table(program = pg, kul_donors = if (is.null(Lk)) 0 else ncol(Lk),
    kul_panel_present = if (is.null(Lk)) 0 else sum(P$panel$gene_symbol %in% rownames(Lk)),
    kul_eval_coverage = if (is.null(Lk)) NA else round(mean(P$evaluation_genes$gene_symbol %in% rownames(Lk)), 3))
  dk <- if (ok_k) cbind(cohort = "KUL", score_cohort(Lk, P)) else NULL
  ds <- cbind(cohort = "SMC", score_cohort(Ls, P))
  if (!is.null(dk)) res[[paste(pg, "KUL")]] <- evaluate(dk, "KUL alone", pg)
  res[[paste(pg, "pooled")]] <- evaluate(rbind(ds, dk), "SMC + KUL pooled", pg)
}
res <- rbindlist(res, fill = TRUE); cover <- rbindlist(cover)
num <- names(res)[vapply(res, is.numeric, logical(1))]; res[, (num) := lapply(.SD, round, 3), .SDcols = num]
fwrite(res, file.path(res_dir, "validation_results_KUL_and_pooled.csv")); fwrite(cover, file.path(res_dir, "KUL_gene_coverage.csv"))
cat("\n=== KUL coverage ===\n"); print(cover)
cat("\n=== KUL alone (descriptive) ===\n"); print(res[analysis == "KUL alone" & scaling == "frozen", .(program, n_donors, MAE_ratio, spearman)])
cat("\n=== POOLED SMC + KUL (frozen scaling) ===\n"); print(res[analysis == "SMC + KUL pooled" & scaling == "frozen", .(program, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, status)])
cat("\n=== POOLED SMC + KUL (within-cohort scaling) ===\n"); print(res[analysis == "SMC + KUL pooled" & scaling == "within", .(program, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, status)])
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_12.txt"))
 
