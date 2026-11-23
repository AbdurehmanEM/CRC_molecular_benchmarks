
# 16_cross_cohort_meta.R  (protocol v1.0, section 10c) - specified after per-cohort results (disclosed)
suppressPackageStartupMessages({ library(data.table); library(digest) })
frz_file <- "04_processed_data/frozen_panel/frozen_panel_v2.rds"
if (digest(file = frz_file, algo = "sha256") != "a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb") stop("Frozen panel hash mismatch - STOP")
frz <- readRDS(frz_file)
res_dir <- "06_results/cross_cohort_meta"; dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
rd <- function(path, cohort, filt = NULL) { x <- fread(path); if (!is.null(filt)) x <- x[analysis == filt]
  x[scaling == "within_cohort", scaling := "within"]
  x[, .(cohort, program, scaling, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, status)] }
R <- rbind(
  rd("06_results/validation_GSE132465/validation_results.csv", "GSE132465"),
  rd("06_results/validation_GSE144735/validation_results_KUL_and_pooled.csv", "GSE144735", "KUL alone"),
  rd("06_results/validation_GSE166555/validation_results.csv", "GSE166555"),
  rd("06_results/validation_GSE200997/validation_results.csv", "GSE200997"), fill = TRUE)
fwrite(R, file.path(res_dir, "per_cohort_results_long.csv"))
dl <- function(est, se) {   # DerSimonian-Laird random effects
  k <- length(est); w <- 1 / se^2; fe <- sum(w * est) / sum(w)
  Q <- sum(w * (est - fe)^2); tau2 <- if (k > 1) max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))) else 0
  wr <- 1 / (se^2 + tau2); m <- sum(wr * est) / sum(wr); s <- sqrt(1 / sum(wr))
  c(k = k, est = m, lo = m - 1.96 * s, hi = m + 1.96 * s, tau2 = tau2,
    I2 = if (k > 1 && Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0) }
meta_one <- function(d) {
  lr <- log(d$MAE_ratio); se_r <- (log(d$ratio_hi) - log(d$ratio_lo)) / (2 * 1.96)
  z <- atanh(pmin(pmax(d$spearman, -0.999), 0.999)); se_z <- 1 / sqrt(d$n_donors - 3)
  a <- dl(lr, se_r); b <- dl(z, se_z)
  data.table(k_cohorts = a[["k"]], cohorts = paste(d$cohort, collapse = ";"), total_donors = sum(d$n_donors),
    n_supported = sum(d$status == "Supported"),
    pooled_MAE_ratio = exp(a[["est"]]), ratio_lo = exp(a[["lo"]]), ratio_hi = exp(a[["hi"]]), I2_ratio = a[["I2"]],
    pooled_rho = tanh(b[["est"]]), rho_lo = tanh(b[["lo"]]), rho_hi = tanh(b[["hi"]]), I2_rho = b[["I2"]]) }
M <- list()
for (an in c("main", "with_GSE144735")) for (sc in c("frozen", "within")) for (pg in names(frz$programs)) {
  d <- R[program == pg & scaling == sc & is.finite(MAE_ratio) & is.finite(ratio_lo) & ratio_lo > 0 & n_donors > 3]
  if (an == "main") d <- d[n_donors >= 10]
  if (nrow(d) < 2) { M[[paste(an, sc, pg)]] <- data.table(analysis = an, scaling = sc, program = pg, k_cohorts = nrow(d)); next }
  M[[paste(an, sc, pg)]] <- cbind(data.table(analysis = an, scaling = sc, program = pg), meta_one(d))
}
M <- rbindlist(M, fill = TRUE)
M[, meta_status := fifelse(is.na(pooled_MAE_ratio), "Fewer than 2 evaluable cohorts",
  fifelse(ratio_hi < 1, "Supported", fifelse(pooled_MAE_ratio < 1, "Partially supported", "Not supported")))]
num <- names(M)[vapply(M, is.numeric, logical(1))]; M[, (num) := lapply(.SD, round, 3), .SDcols = num]
fwrite(M, file.path(res_dir, "meta_analysis_results.csv"))
# ---- Offset diagnostic (descriptive only; no panel change) ----
pbs <- list(GSE132465 = readRDS("04_processed_data/GSE132465/validation_pseudobulk_logcpm.rds")$logcpm,
  GSE144735 = readRDS("04_processed_data/GSE144735/validation_pseudobulk_logcpm.rds"),
  GSE166555 = readRDS("04_processed_data/GSE166555_validation_pseudobulk_logcpm.rds"),
  GSE200997 = readRDS("04_processed_data/GSE200997/GSE200997_validation_pseudobulk_logcpm.rds"))
G <- list(); O <- list()
for (co in names(pbs)) for (pg in names(frz$programs)) {
  P <- frz$programs[[pg]]; L <- pbs[[co]][[P$compartment]]
  if (is.null(L) || !all(P$panel$gene_symbol %in% rownames(L)) || ncol(L) < 3) next
  pn <- P$panel; ev <- P$evaluation_genes[gene_symbol %in% rownames(L)]
  Zx <- sweep(sweep(t(L[pn$gene_symbol, , drop = FALSE]), 2, pn$pred_mean, "-"), 2, pn$pred_sd, "/")
  y <- rowMeans(sweep(sweep(t(L[ev$gene_symbol, , drop = FALSE]), 2, ev$eval_mean, "-"), 2, ev$eval_sd, "/"))
  p <- as.numeric(P$intercept + Zx %*% pn$coefficient)
  G[[paste(co, pg)]] <- data.table(cohort = co, program = pg, gene = pn$gene_symbol, coefficient = pn$coefficient,
    mean_z_vs_discovery = colMeans(Zx), contribution_to_mean_pred = colMeans(Zx) * pn$coefficient)
  O[[paste(co, pg)]] <- data.table(cohort = co, program = pg, n_donors = ncol(L), mean_observed = mean(y),
    mean_predicted = mean(p), mean_offset_pred_minus_obs = mean(p) - mean(y), sd_observed = sd(y))
}
G <- rbindlist(G); O <- rbindlist(O)
rnd <- function(X) { nn <- names(X)[vapply(X, is.numeric, logical(1))]; X[, (nn) := lapply(.SD, round, 3), .SDcols = nn][] }
G <- rnd(G); O <- rnd(O)
fwrite(G, file.path(res_dir, "offset_diagnostic_panel_genes.csv")); fwrite(O, file.path(res_dir, "offset_diagnostic_targets.csv"))
cat("\n=== MAIN meta-analysis (cohorts with >=10 donors), FROZEN scaling ===\n")
print(M[analysis == "main" & scaling == "frozen", .(program, k_cohorts, total_donors, n_supported, pooled_MAE_ratio, ratio_lo, ratio_hi, I2_ratio, pooled_rho, rho_lo, rho_hi, meta_status)])
cat("\n=== MAIN meta-analysis, WITHIN-cohort scaling ===\n")
print(M[analysis == "main" & scaling == "within", .(program, k_cohorts, pooled_MAE_ratio, ratio_lo, ratio_hi, I2_ratio, pooled_rho, meta_status)])
cat("\n=== Sensitivity incl. GSE144735, FROZEN ===\n")
print(M[analysis == "with_GSE144735" & scaling == "frozen", .(program, k_cohorts, pooled_MAE_ratio, ratio_lo, ratio_hi, pooled_rho, meta_status)])
cat("\n=== Offset diagnostic: mean predicted - mean observed (frozen scale) ===\n")
print(dcast(O, program ~ cohort, value.var = "mean_offset_pred_minus_obs"))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_16.txt"))
 
