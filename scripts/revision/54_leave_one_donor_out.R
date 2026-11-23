###############################################################################
# 54_leave_one_donor_out.R - Reviewer point 4 (precision with n = 11)
#
# For every evaluable program x cohort (frozen scaling), recompute the MAE
# ratio leaving out one donor at a time. Reports the range, the most
# influential donor, and a jackknife 95% CI on the log scale, next to the
# published bootstrap CI. The endothelial programs (GSE132465, n = 11) are the
# focus; all others are reported for completeness. The frozen panel is never
# refitted, so this measures donor influence only, not between-cohort variation.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision(); require_gate()

pr <- readRDS(file.path(OUT_DIR, "52_donor_predictions.rds"))$preds
rec <- fread(file.path(OUT_DIR, "52_validation_recomputed.csv"))[scaling == "frozen"]

out <- pr[, {
  n <- .N; full <- mae_ratio(y, yhat)
  loo <- vapply(seq_len(n), function(i) mae_ratio(y[-i], yhat[-i]), numeric(1))
  lj <- log(loo); jk_se <- sqrt((n - 1) / n * sum((lj - mean(lj))^2))
  infl <- donor[which.max(abs(loo - full))]
  .(n = n, ratio = full, loo_min = min(loo), loo_max = max(loo),
    most_influential = infl, ratio_without_it = loo[donor == infl],
    jk_lo = exp(log(full) - 1.959964 * jk_se), jk_hi = exp(log(full) + 1.959964 * jk_se),
    all_loo_below_1 = all(loo < 1))
}, by = .(program, cohort)]
out <- merge(out, rec[, .(program, cohort, boot_lo = lo, boot_hi = hi, status)], by = c("program", "cohort"))
num <- names(out)[sapply(out, is.numeric)]; out[, (num) := lapply(.SD, round, 3), .SDcols = num]
write_out(out[order(program, cohort)], "54_leave_one_donor_out")

## Per-donor detail for the endothelial programs (supplementary)
endo <- pr[grepl("^END_", program) & cohort == "GSE132465"]
if (nrow(endo)) {
  det <- endo[, .(donor, ratio_without = vapply(seq_len(.N), function(i) mae_ratio(y[-i], yhat[-i]), numeric(1)),
                  abs_error = abs(yhat - y), target = y, prediction = yhat), by = program]
  write_out(det, "54_endothelial_donor_detail")
}
print(out[n >= MIN_DONORS])
cat("\nSentence for Results (point 4), fill from the END_TNFA_NFKB row:\n",
    "Leaving out each of the 11 GSE132465 patients in turn gave MAE ratios of [loo_min]-[loo_max] for\n",
    "endothelial TNF/NF-kB (jackknife 95% CI [jk_lo]-[jk_hi]); removing the most influential patient\n",
    "([most_influential]) gave [ratio_without_it].\n")
