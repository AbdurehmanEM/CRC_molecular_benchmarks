###############################################################################
# 56_calibration_transport.R - Reviewer point 8 (external platform transport)
#
# Quantifies how much of the frozen-scale error in each validation cohort is
# a cohort-level offset or slope (calibration) rather than loss of ranking:
#   offset       = mean(prediction - target)
#   intercept/slope of lm(target ~ prediction)  (ideal: 0 / 1)
#   MAE ratio, frozen (published primary; baseline = discovery mean = 0)
#   MAE ratio after leave-one-out OFFSET recalibration
#   MAE ratio after leave-one-out LINEAR recalibration
# Recalibrated ratios use a leave-one-out cohort-mean baseline, so they ask
# "does the panel beat the cohort's own mean once calibrated locally?" - the
# situation of a chip study with within-experiment reference conditions.
# No information from donor i is used to recalibrate donor i.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision(); require_gate()

pr <- readRDS(file.path(OUT_DIR, "52_donor_predictions.rds"))$preds
out <- pr[, {
  n <- .N
  fit <- lm(y ~ yhat); ci <- suppressWarnings(confint(fit))
  off_loo <- vapply(seq_len(n), function(i) yhat[i] - mean(yhat[-i] - y[-i]), numeric(1))
  lin_loo <- vapply(seq_len(n), function(i) { f <- lm(y[-i] ~ yhat[-i]); sum(coef(f) * c(1, yhat[i])) }, numeric(1))
  base_loo <- vapply(seq_len(n), function(i) mean(y[-i]), numeric(1))
  mb <- mean(abs(y - base_loo))
  .(n = n, offset = mean(yhat - y),
    cal_intercept = coef(fit)[1], cal_slope = coef(fit)[2], slope_lo = ci[2, 1], slope_hi = ci[2, 2],
    rho = cor(y, yhat, method = "spearman"),
    ratio_frozen = mae_ratio(y, yhat),
    ratio_offset_loo = mean(abs(off_loo - y)) / mb,
    ratio_linear_loo = mean(abs(lin_loo - y)) / mb,
    share_error_offset = 1 - mean(abs(off_loo - y)) / mean(abs(yhat - y)))
}, by = .(program, cohort)]
num <- names(out)[sapply(out, is.numeric)]; out[, (num) := lapply(.SD, round, 3), .SDcols = num]
write_out(out[order(program, cohort)], "56_calibration_transport")

## Calibration plots (one page per program)
pdf(file.path(FIG_DIR, "56_calibration_plots.pdf"), width = 8, height = 3)
for (p in unique(pr$program)) {
  d <- pr[program == p]; cs <- unique(d$cohort)
  par(mfrow = c(1, length(cs)), mar = c(4, 4, 2, 1), cex = 0.7)
  for (co in cs) {
    x <- d[cohort == co]; lim <- range(c(x$y, x$yhat))
    plot(x$yhat, x$y, xlim = lim, ylim = lim, pch = 19, col = "grey30",
         xlab = "Frozen prediction", ylab = "Observed target", main = paste(p, co, "n =", nrow(x)))
    abline(0, 1, lty = 2); if (nrow(x) > 2) abline(lm(y ~ yhat, data = x), col = "firebrick")
  }
}
invisible(dev.off())
print(out[n >= MIN_DONORS])
cat("\nReport (point 8): offsets, slopes and the frozen vs locally recalibrated MAE ratios; state that\n",
    "recalibration uses each validation cohort's own data (leave-one-out) and is therefore not an\n",
    "independent test of the frozen scale but an estimate of performance after on-chip reference calibration.\n")
