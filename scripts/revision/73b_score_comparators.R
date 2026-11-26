###############################################################################
# 73b_score_comparators.R - Protocol amendment v2.3, PART B (external scoring)
#
# Scores the frozen panels and the frozen comparators (73a) in every validation
# cohort with identical targets, frozen discovery scaling, donor eligibility and
# missing-gene rule (same code path as the published analysis, script 52).
#
# Outputs (06_results/revision_2026-11):
#   73b_external_by_cohort.csv  MAE ratio (95% donor bootstrap CI), Spearman, status,
#                               and the PAIRED difference log(MAE panel) - log(MAE comparator)
#                               on the same bootstrap resamples (< 0 = panel better)
#   73b_external_pooled.csv     DerSimonian-Laird pooling over main cohorts (>= 10 donors),
#                               as in script 52/16; single-cohort programs reported as k = 1
# Descriptive: comparison is reported in full; it changes no panel, category or headline.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision(); require_gate()
suppressPackageStartupMessages(library(digest))

FROZEN_COMP_FILE <- file.path(PROJECT_ROOT, "04_processed_data", "frozen_panel", "frozen_comparators_v2.3.rds")
stopifnot(file.exists(FROZEN_COMP_FILE),
          digest(file = FROZEN_COMP_FILE, algo = "sha256") == readLines(sub("\\.rds$", ".sha256", FROZEN_COMP_FILE))[1])
for (k in names(FROZEN_FILES))
  stopifnot(digest(file = FROZEN_FILES[[k]], algo = "sha256") == FROZEN_SHA256[[k]])
FC <- readRDS(FROZEN_COMP_FILE)

panels <- load_panels(); partition <- load_partition()
pc <- program_compartment(panels); cmps <- unique(pc)
cohorts <- c(VALIDATION, DESCRIPTIVE)
lc <- build_lc_cache(c(DISCOVERY, cohorts), cmps)
dstat <- all_disc_stats(panels, partition, lc)
pub <- fread(file.path(OUT_DIR, "52_validation_recomputed.csv"))[scaling == "frozen"]

status_of <- function(n, r, hi) if (n < MIN_DONORS) "Descriptive only (<10 donors)" else
  if (r < 1 && hi < 1) "Supported" else if (r < 1) "Partially supported" else "Not supported"

rows <- list()
for (p in names(panels)) for (co in cohorts) {
  e0 <- eval_program_cohort(p, co, panels, partition, dstat, lc)
  if (is.null(e0) || is.null(e0$yhat) || e0$coverage < MIN_EVAL_COVER) next
  methods <- list(panel = panels[[p]], C1_DE = FC$comparators[[p]]$C1_DE, C2_CORR = FC$comparators[[p]]$C2_CORR)
  preds <- lapply(methods, function(m) if (is.null(m)) NULL else
    eval_program_cohort(p, co, panels, partition, dstat, lc, panel = m)$yhat)
  y <- unname(e0$y); n <- length(y)
  # Sanity: panel must reproduce the published (script 52) frozen ratio
  pr <- pub[program == p & cohort == co, ratio]
  if (length(pr) == 1 && is.finite(pr)) stopifnot(abs(mae_ratio(y, preds$panel) - pr) < 1e-6)
  set.seed(SEED)
  idx <- replicate(N_BOOT, sample.int(n, n, replace = TRUE), simplify = FALSE)
  rb <- function(yh) vapply(idx, function(i) mae_ratio(y[i], yh[i]), numeric(1))
  bp <- rb(preds$panel)
  for (m in names(methods)) {
    yh <- preds[[m]]
    if (is.null(yh)) {
      rows[[length(rows) + 1]] <- data.table(program = p, compartment = pc[[p]], cohort = co, method = m, donors = n,
        genes = paste(methods[[m]]$symbol, collapse = ", "), note = "not evaluable (gene missing or none selected)")
      next
    }
    r <- mae_ratio(y, yh); b <- rb(yh); ok <- is.finite(b) & is.finite(bp)
    ci <- quantile(b[ok], c(0.025, 0.975), names = FALSE)
    d  <- log(mae_ratio(y, preds$panel)) - log(r); db <- log(bp[ok]) - log(b[ok])
    rows[[length(rows) + 1]] <- data.table(
      program = p, compartment = pc[[p]], cohort = co, method = m, donors = n,
      genes = paste(methods[[m]]$symbol, collapse = ", "),
      ratio = r, lo = ci[1], hi = ci[2], rho = cor(y, yh, method = "spearman"),
      status = status_of(n, r, ci[2]),
      d_log_panel_minus_method = if (m == "panel") NA_real_ else d,
      d_lo = if (m == "panel") NA_real_ else quantile(db, 0.025, names = FALSE),
      d_hi = if (m == "panel") NA_real_ else quantile(db, 0.975, names = FALSE),
      d_se = if (m == "panel") NA_real_ else sd(db),
      note = if (co %in% DESCRIPTIVE) "descriptive cohort, not pooled" else "")
  }
}
res <- rbindlist(rows, fill = TRUE)
write_out(res, "73b_external_by_cohort")

## ---- Pooling over main cohorts (>= 10 donors), as script 52 -----------------
main <- res[cohort %in% VALIDATION & donors >= MIN_DONORS & is.finite(ratio) & lo > 0]
pool_ratio <- main[, { d <- dl_pool(log(ratio), se_from_ci(lo, hi))
  .(k = .N, donors = sum(donors), pooled_ratio = exp(d$est), plo = exp(d$lo), phi = exp(d$hi), I2 = d$I2) },
  by = .(program, compartment, method)]
pool_diff <- main[method != "panel" & is.finite(d_se) & d_se > 0, { d <- dl_pool(d_log_panel_minus_method, d_se)
  .(pooled_ratio_of_MAE = exp(d$est), rlo = exp(d$lo), rhi = exp(d$hi), I2_diff = d$I2) },
  by = .(program, method)]
pooled <- merge(pool_ratio, pool_diff, by = c("program", "method"), all.x = TRUE)
## k = 1 (e.g. endothelial programs, GSE132465 only): report the cohort's own bootstrap values
for (i in which(pooled$k == 1)) {
  s <- main[program == pooled$program[i] & method == pooled$method[i]]
  set(pooled, i, c("pooled_ratio", "plo", "phi"), list(s$ratio, s$lo, s$hi))
  if (pooled$method[i] != "panel")
    set(pooled, i, c("pooled_ratio_of_MAE", "rlo", "rhi"),
        list(exp(s$d_log_panel_minus_method), exp(s$d_lo), exp(s$d_hi)))
}
pooled[, verdict := fifelse(method == "panel", "",
          fifelse(rhi < 1, "panel lower error",
          fifelse(rlo > 1, "comparator lower error", "no clear difference")))]
write_out(pooled, "73b_external_pooled")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "73b_sessionInfo.txt"))

cat("\npooled_ratio_of_MAE = MAE(panel) / MAE(comparator); < 1 favours the frozen panel.\n")
print(pooled[order(program, method), .(program, method, k, donors, pooled = round(pooled_ratio, 3),
                                       ci = sprintf("%.3f-%.3f", plo, phi),
                                       panel_vs = round(pooled_ratio_of_MAE, 3),
                                       vs_ci = sprintf("%.3f-%.3f", rlo, rhi), verdict)])
cat("\nReport every row, including programs where a comparator has the lower error.\n")
