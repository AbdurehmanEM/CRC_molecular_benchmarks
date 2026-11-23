
# 35_lymphoid_meta.R  (protocol v1.9, section 14.10) - pooled lymphoid validation.
# dl() and meta_one() copied verbatim from 16_cross_cohort_meta.R.
suppressPackageStartupMessages({ library(data.table); library(digest) })
if (digest(file = "04_processed_data/frozen_panel/frozen_panel_LYM_v1.rds", algo = "sha256") !=
    "7dc3256d97477a796fcacd1df71021740888f7978a01ea36caf70787da7f513f") stop("Lymphoid panel hash mismatch - STOP")
res_dir <- "06_results/validation_lymphoid"

a <- fread("06_results/validation_GSE132465/validation_results_lymphoid.csv")
a[scaling == "within_cohort", scaling := "within"]; a[, cohort := "GSE132465"]
b <- fread(file.path(res_dir, "validation_results_extra_cohorts.csv"))
R <- rbind(a[, .(cohort, program, scaling, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, status)],
           b[, .(cohort, program, scaling, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, status)], fill = TRUE)
fwrite(R, file.path(res_dir, "per_cohort_results_long_lymphoid.csv"))
cat("=== All cohorts attempted (reported in full, protocol v1.9 s14.10.3) ===\n")
print(R[scaling == "frozen", .(cohort, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, status)])

dl <- function(est, se) { k <- length(est); w <- 1 / se^2; fe <- sum(w * est) / sum(w)
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
    pooled_MAE_ratio = exp(a[["est"]]), ratio_lo = exp(a[["lo"]]), ratio_hi = exp(a[["hi"]]),
    tau2_ratio = a[["tau2"]], I2_ratio = a[["I2"]],
    pooled_rho = tanh(b[["est"]]), rho_lo = tanh(b[["lo"]]), rho_hi = tanh(b[["hi"]]), I2_rho = b[["I2"]]) }

M <- list()
for (sc in c("frozen", "within")) {
  d <- R[scaling == sc & is.finite(MAE_ratio) & is.finite(ratio_lo) & ratio_lo > 0 & n_donors >= 10]
  M[[sc]] <- if (nrow(d) < 2) data.table(scaling = sc, program = "LYM_IL2_STAT5", k_cohorts = nrow(d))
             else cbind(data.table(scaling = sc, program = "LYM_IL2_STAT5"), meta_one(d)) }
M <- rbindlist(M, fill = TRUE)
M[, meta_status := fifelse(is.na(pooled_MAE_ratio), "Fewer than 2 evaluable cohorts",
    fifelse(ratio_hi < 1, "Supported", fifelse(pooled_MAE_ratio < 1, "Partially supported", "Not supported")))]
num <- names(M)[vapply(M, is.numeric, logical(1))]; M[, (num) := lapply(.SD, round, 3), .SDcols = num]
fwrite(M, file.path(res_dir, "meta_analysis_results_lymphoid.csv"))

# offset diagnostic from saved per-donor predictions
p1 <- fread("06_results/validation_GSE132465/validation_predictions_lymphoid.csv")[scaling == "frozen",
        .(cohort = "GSE132465", obs = observed, pred = predicted)]
p2 <- fread(file.path(res_dir, "validation_predictions_extra_cohorts.csv"))[, .(cohort, obs = y_frozen, pred = p_frozen)]
O <- rbind(p1, p2)[, .(n = .N, mean_observed = mean(obs), mean_predicted = mean(pred),
                       mean_offset = mean(pred) - mean(obs), sd_observed = sd(obs)), by = cohort]
nn <- names(O)[vapply(O, is.numeric, logical(1))]; O[, (nn) := lapply(.SD, round, 3), .SDcols = nn]
fwrite(O, file.path(res_dir, "offset_diagnostic_lymphoid.csv"))

cat("\n=== POOLED lymphoid validation ===\n")
print(M[, .(scaling, k_cohorts, total_donors, n_supported, pooled_MAE_ratio, ratio_lo, ratio_hi,
            tau2_ratio, I2_ratio, pooled_rho, rho_lo, rho_hi, meta_status)])
cat("\n=== Offset diagnostic (frozen scale) ===\n"); print(O)
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_35.txt"))
 
