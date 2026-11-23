# 11_validation_GSE132465.R
# Protocol v0.9: apply the FROZEN panel (frozen_panel_v2) to independent validation cohort GSE132465.
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(digest); library(ggplot2) })
set.seed(20261008)
MIN_CELLS <- 20; NBOOT <- 2000; MIN_DONORS <- 10; MIN_EVAL_COVER <- 0.80
frz_file <- "04_processed_data/frozen_panel/frozen_panel_v2.rds"
EXPECTED_SHA <- "a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb"
res_dir <- "06_results/validation_GSE132465"; fig_dir <- "07_figures/validation"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
sha <- digest(file = frz_file, algo = "sha256")
if (sha != EXPECTED_SHA) stop("Frozen panel hash does not match - STOP. Got: ", sha)
frz <- readRDS(frz_file)
cat("Frozen panel verified:", sha, "\n")
counts <- readRDS("04_processed_data/GSE132465/GSE132465_full_counts_sparse.rds")
meta   <- fread("06_results/GSE132465/cell_meta_with_qc_and_compartment.csv.gz")
meta   <- meta[match(colnames(counts), Index)]; stopifnot(identical(meta$Index, colnames(counts)))
meta[, val_compartment := fcase(
  compartment == "Epithelial" & epi_status == "putative_malignant",            "Epithelial",
  compartment == "Stromal",                                                    "Stromal",
  compartment == "Endothelial",                                                "Endothelial",
  Cell_type == "Myeloids" & !Cell_subtype %in% c("cDC", "Unknown"),            "Myeloid",
  default = NA_character_)]
use <- meta$Class == "Tumor" & !meta$qc_flag & !is.na(meta$val_compartment)
cat("Validation cells used:", sum(use), "\n")
logcpm <- list(); ncells <- list()
for (comp in unique(na.omit(meta$val_compartment))) {
  ix <- which(use & meta$val_compartment == comp)
  f2 <- factor(meta$Patient[ix])
  ind <- sparseMatrix(i = seq_along(ix), j = as.integer(f2), x = 1, dims = c(length(ix), nlevels(f2)),
                      dimnames = list(NULL, levels(f2)))
  pb <- as.matrix(counts[, ix] %*% ind)
  nc <- table(f2)
  keep <- names(nc)[nc >= MIN_CELLS]
  logcpm[[comp]] <- cpm(pb[, keep, drop = FALSE], log = TRUE, prior.count = 2)
  ncells[[comp]] <- nc
  cat(comp, ":", length(keep), "donors with >=", MIN_CELLS, "cells\n")
}
saveRDS(list(logcpm = logcpm, ncells = ncells), "04_processed_data/GSE132465/validation_pseudobulk_logcpm.rds")
 
boot_ci <- function(obs, pred, base) {
  n <- length(obs); r <- replicate(NBOOT, { i <- sample(n, n, TRUE)
    c(mean(abs(obs[i] - pred[i])) / mean(abs(obs[i] - base)),
      suppressWarnings(cor(obs[i], pred[i], method = "spearman"))) })
  apply(r, 1, quantile, c(0.025, 0.975), na.rm = TRUE)
}
res <- list(); cover <- list(); preds_all <- list()
for (pg in names(frz$programs)) {
  P <- frz$programs[[pg]]; comp <- P$compartment
  L <- logcpm[[comp]]
  panel_ok <- P$panel$gene_symbol %in% rownames(L)
  ev_ok    <- P$evaluation_genes$gene_symbol %in% rownames(L)
  cover[[pg]] <- data.table(program = pg, compartment = comp, n_donors = ncol(L),
                            panel_genes = nrow(P$panel), panel_present = sum(panel_ok),
                            missing_panel = paste(P$panel$gene_symbol[!panel_ok], collapse = ";"),
                            eval_genes = nrow(P$evaluation_genes), eval_present = sum(ev_ok),
                            eval_coverage = round(mean(ev_ok), 3))
  evaluable <- all(panel_ok) && mean(ev_ok) >= MIN_EVAL_COVER && ncol(L) >= MIN_DONORS
  if (!evaluable) { res[[pg]] <- data.table(program = pg, scaling = c("frozen", "within_cohort"),
                                            status = "Not evaluable"); next }
  ev <- P$evaluation_genes[ev_ok]; pn <- P$panel
  E <- t(L[ev$gene_symbol, , drop = FALSE]); X <- t(L[pn$gene_symbol, , drop = FALSE])
  for (sc in c("frozen", "within_cohort")) {
    if (sc == "frozen") {
      y  <- rowMeans(sweep(sweep(E, 2, ev$eval_mean, "-"), 2, ev$eval_sd, "/"))
      Z  <- sweep(sweep(X, 2, pn$pred_mean, "-"), 2, pn$pred_sd, "/")
    } else {
      y  <- rowMeans(scale(E)); Z <- scale(X)
    }
    pred <- as.numeric(P$intercept + Z %*% pn$coefficient)
    base <- 0
    mae  <- mean(abs(y - pred)); mae_b <- mean(abs(y - base))
    ci   <- boot_ci(y, pred, base)
    cal  <- coef(lm(y ~ pred))
    ratio <- mae / mae_b
    status <- if (ratio >= 1) "Not supported" else if (ci[2, 1] < 1) "Supported" else "Partially supported"
    res[[paste(pg, sc)]] <- data.table(program = pg, scaling = sc, n_donors = length(y),
      MAE_panel = mae, MAE_baseline = mae_b, MAE_ratio = ratio, ratio_lo = ci[1, 1], ratio_hi = ci[2, 1],
      spearman = suppressWarnings(cor(y, pred, method = "spearman")), rho_lo = ci[1, 2], rho_hi = ci[2, 2],
      cal_intercept = cal[1], cal_slope = cal[2], status = status)
    preds_all[[paste(pg, sc)]] <- data.table(program = pg, scaling = sc, Patient = rownames(E), observed = y, predicted = pred)
  }
}
 
res <- rbindlist(res, fill = TRUE); cover <- rbindlist(cover); preds_all <- rbindlist(preds_all)
num <- names(res)[vapply(res, is.numeric, logical(1))]
res[, (num) := lapply(.SD, round, 3), .SDcols = num]
fwrite(res, file.path(res_dir, "validation_results.csv"))
fwrite(cover, file.path(res_dir, "validation_gene_coverage.csv"))
fwrite(preds_all, file.path(res_dir, "validation_predictions.csv"))
cat("\n=== Gene coverage in GSE132465 ===\n"); print(cover)
cat("\n=== Validation results (primary = frozen scaling) ===\n")
print(res[scaling == "frozen", .(program, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, rho_lo, rho_hi, cal_slope, status)])
cat("\n=== Secondary: within-cohort standardisation ===\n")
print(res[scaling == "within_cohort", .(program, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, cal_slope, status)])
if (nrow(preds_all)) {
  p <- ggplot(preds_all[scaling == "frozen"], aes(predicted, observed)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_point(size = 1.6, alpha = 0.8) + geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
    facet_wrap(~program, scales = "free") + theme_bw(10) +
    labs(x = "Predicted program score (frozen panel)", y = "Observed program score",
         title = "Independent validation: GSE132465 (frozen panel v2)")
  ggsave(file.path(fig_dir, "validation_GSE132465_obs_vs_pred.pdf"), p, width = 10, height = 7)
}
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_11.txt"))
 
