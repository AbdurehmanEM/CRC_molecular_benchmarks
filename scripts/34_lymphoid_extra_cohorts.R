
# 34_lymphoid_extra_cohorts.R  (protocol v1.9, section 14.10) - frozen lymphoid panel in GSE166555 and GSE144735.
# CD8 definitions fixed before running:
#   GSE166555: minor_cell_type matching ^CD8\+ T   (Tconv/Treg/other excluded)
#   GSE144735: Cell_type == "T cells" & Cell_subtype == "CD8+ T cells"
# GSE200997 EXCLUDED: its annotation resolves T and NK jointly (ctx_T_NK); CD8 is not separable.
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(digest) })
set.seed(20261008)
MIN_CELLS <- 20; NBOOT <- 2000; MIN_DONORS <- 10; MIN_EVAL_COVER <- 0.80
frz_file <- "04_processed_data/frozen_panel/frozen_panel_LYM_v1.rds"
if (digest(file = frz_file, algo = "sha256") != "7dc3256d97477a796fcacd1df71021740888f7978a01ea36caf70787da7f513f")
  stop("Lymphoid frozen panel hash mismatch - STOP")
frz <- readRDS(frz_file); cat("Lymphoid frozen panel verified\n")
mad_flag <- function(x, n = 3, side = c("both", "high")) { side <- match.arg(side); m <- median(x); s <- mad(x)
  if (!is.finite(s) || s == 0) return(rep(FALSE, length(x)))
  if (side == "high") x > m + n * s else x < m - n * s | x > m + n * s }
pseudo <- function(counts, donor, ix) {
  f2 <- factor(donor[ix])
  ind <- sparseMatrix(i = seq_along(ix), j = as.integer(f2), x = 1, dims = c(length(ix), nlevels(f2)),
                      dimnames = list(NULL, levels(f2)))
  pb <- as.matrix(counts[, ix] %*% ind); nc <- table(f2); keep <- names(nc)[nc >= MIN_CELLS]
  if (!length(keep)) return(NULL)
  edgeR::cpm(pb[, keep, drop = FALSE], log = TRUE, prior.count = 2) }
score_cohort <- function(L, P) {
  ev <- P$evaluation_genes[gene_symbol %in% rownames(L)]; pn <- P$panel
  E <- t(L[ev$gene_symbol, , drop = FALSE]); X <- t(L[pn$gene_symbol, , drop = FALSE])
  zs <- function(M) { Z <- scale(M); Z[!is.finite(Z)] <- 0; Z }
  data.table(donor = rownames(E),
    y_frozen = rowMeans(sweep(sweep(E, 2, ev$eval_mean, "-"), 2, ev$eval_sd, "/")),
    p_frozen = as.numeric(P$intercept + sweep(sweep(X, 2, pn$pred_mean, "-"), 2, pn$pred_sd, "/") %*% pn$coefficient),
    y_within = rowMeans(zs(E)), p_within = as.numeric(P$intercept + zs(X) %*% pn$coefficient)) }
evaluate <- function(L, cohort) {
  out <- list(); P <- frz$programs[["LYM_IL2_STAT5"]]
  nd <- if (is.null(L)) 0 else ncol(L)
  pan <- if (is.null(L)) 0 else sum(P$panel$gene_symbol %in% rownames(L))
  evc <- if (is.null(L)) NA_real_ else mean(P$evaluation_genes$gene_symbol %in% rownames(L))
  cov <- data.table(cohort = cohort, program = "LYM_IL2_STAT5", n_donors = nd, panel_genes = nrow(P$panel),
    panel_present = pan, missing_panel = if (is.null(L)) "" else paste(setdiff(P$panel$gene_symbol, rownames(L)), collapse = ";"),
    eval_coverage = round(evc, 3))
  if (is.null(L) || pan < nrow(P$panel) || evc < MIN_EVAL_COVER || nd < MIN_DONORS) {
    return(list(res = data.table(cohort = cohort, program = "LYM_IL2_STAT5", scaling = c("frozen", "within"),
                                 n_donors = nd, status = "Not evaluable"), cover = cov, preds = NULL)) }
  d <- score_cohort(L, P)
  for (sc in c("frozen", "within")) {
    y <- d[[paste0("y_", sc)]]; p <- d[[paste0("p_", sc)]]; n <- length(y)
    ratio <- mean(abs(y - p)) / mean(abs(y))
    bs <- replicate(NBOOT, { i <- sample.int(n, n, TRUE)
      c(mean(abs(y[i] - p[i])) / mean(abs(y[i])), suppressWarnings(cor(y[i], p[i], method = "spearman"))) })
    ci <- apply(bs, 1, quantile, c(0.025, 0.975), na.rm = TRUE)
    status <- if (!is.finite(ratio) || ratio >= 1) "Not supported" else if (ci[2, 1] < 1) "Supported" else "Partially supported"
    out[[sc]] <- data.table(cohort = cohort, program = "LYM_IL2_STAT5", scaling = sc, n_donors = n,
      MAE_ratio = ratio, ratio_lo = ci[1, 1], ratio_hi = ci[2, 1],
      spearman = suppressWarnings(cor(y, p, method = "spearman")), rho_lo = ci[1, 2], rho_hi = ci[2, 2],
      cal_slope = if (sd(p) > 0) unname(coef(lm(y ~ p))[2]) else NA_real_, status = status) }
  list(res = rbindlist(out), cover = cov, preds = cbind(cohort = cohort, d)) }

## ---------------- GSE166555 ----------------
in1 <- "03_original_data/GSE166555"; raw1 <- "04_processed_data/GSE166555/raw"   # same extraction folder as script 13
dir.create(raw1, recursive = TRUE, showWarnings = FALSE)
meta1 <- fread(file.path(in1, "GSE166555_meta_data.tsv.gz"))
meta1[, is_cd8 := grepl("^CD8\\+ T", minor_cell_type)]
cat("GSE166555 CD8 cells in metadata:", sum(meta1$is_cd8), "\n")
tum <- grep("t[0-9]*\\.tsv\\.gz$", untar(file.path(in1, "GSE166555_RAW.tar"), list = TRUE), value = TRUE)
if (!all(file.exists(file.path(raw1, tum)))) untar(file.path(in1, "GSE166555_RAW.tar"), files = tum, exdir = raw1)
mats <- list()
for (fl in tum) { x <- fread(file.path(raw1, fl)); g <- x[[1]]
  m <- as(as.matrix(x[, -1]), "CsparseMatrix"); rownames(m) <- g; mats[[fl]] <- m; message("read ", fl) }
gg <- Reduce(intersect, lapply(mats, rownames))
c1 <- do.call(cbind, lapply(mats, function(m) m[gg, , drop = FALSE])); rm(mats); invisible(gc())
if (anyDuplicated(colnames(c1))) stop("duplicated cell barcodes across GSE166555 sample files")
cm <- meta1[match(colnames(c1), cell)]
mt <- grepl("^MT-", gg)
cm[, nCount := Matrix::colSums(c1)]; cm[, nFeature := as.numeric(diff(c1@p))]
cm[, pct_mt := 100 * Matrix::colSums(c1[mt, , drop = FALSE]) / pmax(nCount, 1)]
cm[, qc_flag := mad_flag(log10(nCount + 1)) | mad_flag(log10(nFeature + 1)) | mad_flag(pct_mt, side = "high"), by = sample_id]
u1 <- which(cm$is_cd8 %in% TRUE & !cm$qc_flag & cm$sample_origin == "Tumor")
cat("GSE166555 CD8 cells used:", length(u1), "| patients:", uniqueN(cm$case_id[u1]), "\n")
L1 <- pseudo(c1, cm$case_id, u1); rm(c1); invisible(gc())
r1 <- evaluate(L1, "GSE166555")

## ---------------- GSE144735 (KUL) ----------------
c2 <- readRDS("04_processed_data/GSE144735/GSE144735_counts_sparse.rds")
a2 <- fread("03_original_data/GSE144735/GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz")
m2 <- a2[match(colnames(c2), Index)]
mt2 <- grepl("^MT-", rownames(c2))
m2[, nCount := Matrix::colSums(c2)]; m2[, nFeature := Matrix::colSums(c2 > 0)]
m2[, pct_mt := 100 * Matrix::colSums(c2[mt2, , drop = FALSE]) / nCount]
m2[, qc_flag := mad_flag(log10(nCount + 1)) | mad_flag(log10(nFeature + 1)) | mad_flag(pct_mt, side = "high"), by = Sample]
u2 <- which(m2$Class == "Tumor" & !m2$qc_flag & m2$Cell_type == "T cells" & m2$Cell_subtype == "CD8+ T cells")
cat("KUL CD8 cells used:", length(u2), "| patients:", uniqueN(m2$Patient[u2]), "\n")
L2 <- pseudo(c2, m2$Patient, u2); rm(c2); invisible(gc())
r2 <- evaluate(L2, "GSE144735_KUL")

res <- rbindlist(list(r1$res, r2$res), fill = TRUE)
cover <- rbindlist(list(r1$cover, r2$cover))
preds <- rbindlist(Filter(Negate(is.null), list(r1$preds, r2$preds)), fill = TRUE)
res <- rbind(res, data.table(cohort = "GSE200997", program = "LYM_IL2_STAT5", scaling = c("frozen", "within"),
                             n_donors = 0, status = "Excluded: annotation resolves T and NK jointly (ctx_T_NK)"), fill = TRUE)
num <- names(res)[vapply(res, is.numeric, logical(1))]; res[, (num) := lapply(.SD, round, 3), .SDcols = num]
dir.create("06_results/validation_lymphoid", recursive = TRUE, showWarnings = FALSE)
fwrite(res, "06_results/validation_lymphoid/validation_results_extra_cohorts.csv")
fwrite(cover, "06_results/validation_lymphoid/validation_coverage_extra_cohorts.csv")
if (nrow(preds)) fwrite(preds, "06_results/validation_lymphoid/validation_predictions_extra_cohorts.csv")
cat("\n=== Coverage ===\n"); print(cover)
cat("\n=== Extra cohorts, frozen scaling (primary) ===\n")
print(res[scaling == "frozen", .(cohort, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, cal_slope, status)])
cat("\n=== Extra cohorts, within-cohort scaling ===\n")
print(res[scaling == "within", .(cohort, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, status)])
writeLines(capture.output(sessionInfo()), "06_results/validation_lymphoid/sessionInfo_34.txt")
 
