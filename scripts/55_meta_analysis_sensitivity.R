###############################################################################
# 55_meta_analysis_sensitivity.R - Reviewer point 6 (heterogeneity, k = 2-3)
#
# Re-pools the published per-cohort MAE ratios (Supplementary Table S8) with
#   DL           DerSimonian-Laird (published method; reproduced as a check)
#   REML         restricted maximum likelihood
#   REML + HK    REML with Hartung-Knapp(-Sidik-Jonkman) t-based CI
#   PM + HK      Paule-Mandel tau2 with Hartung-Knapp
#   HK ad hoc    Hartung-Knapp with the "ad hoc" correction (never narrower than DL)
# plus leave-one-cohort-out and the 95% prediction interval (k >= 3).
# Works from the workbook alone (no pseudobulk needed). Uses log ratios with
# SE = (log U - log L) / (2 x 1.96), as in the published pooling.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision()
suppressPackageStartupMessages(library(metafor))

s8 <- read_supp("S8_validation_by_cohort")
setnames(s8, c("Program code", "Cohort", "Scaling", "Donors", "MAE ratio", "95% CI, lower", "95% CI, upper"),
         c("program", "cohort", "scaling", "donors", "ratio", "lo", "hi"))
s8[, scaling := ifelse(grepl("^Frozen", scaling), "frozen", "within")]
num <- c("donors", "ratio", "lo", "hi"); s8[, (num) := lapply(.SD, as.numeric), .SDcols = num]
s8 <- s8[is.finite(ratio) & is.finite(lo) & lo > 0]
s8[, `:=`(yi = log(ratio), sei = se_from_ci(lo, hi))]

fit1 <- function(d, method, test = "z") {
  f <- tryCatch(rma(yi = d$yi, sei = d$sei, method = method, test = test), error = function(e) NULL)
  if (is.null(f)) return(data.table(est = NA_real_, lo = NA_real_, hi = NA_real_, tau2 = NA_real_, I2 = NA_real_, pi_lo = NA_real_, pi_hi = NA_real_))
  pr <- if (d[, .N] >= 3) predict(f) else NULL
  data.table(est = exp(f$b[1]), lo = exp(f$ci.lb), hi = exp(f$ci.ub), tau2 = f$tau2, I2 = f$I2,
             pi_lo = if (!is.null(pr)) exp(pr$pi.lb) else NA_real_, pi_hi = if (!is.null(pr)) exp(pr$pi.ub) else NA_real_)
}
methods <- list(DL = c("DL", "z"), REML = c("REML", "z"), `REML + HK` = c("REML", "knha"),
                `PM + HK` = c("PM", "knha"), `DL + HK ad hoc` = c("DL", "adhoc"))

res <- list()
for (sc in c("frozen", "within")) for (set in c("main", "incl_GSE144735")) {
  dd <- if (set == "main") s8[scaling == sc & cohort %in% VALIDATION & donors >= MIN_DONORS]
        else s8[scaling == sc & cohort %in% c(VALIDATION, DESCRIPTIVE) & (donors >= MIN_DONORS | cohort %in% DESCRIPTIVE)]
  for (p in unique(dd$program)) {
    d <- dd[program == p]; if (nrow(d) < 2) next
    for (m in names(methods)) {
      r <- fit1(d, methods[[m]][1], methods[[m]][2])
      res[[length(res) + 1]] <- cbind(data.table(program = p, scaling = sc, cohorts = set, k = nrow(d),
                                                 analysis = "all cohorts", method = m), r)
    }
    if (nrow(d) >= 3) for (drop in d$cohort) {          # leave-one-cohort-out
      dl1 <- d[cohort != drop]
      for (m in c("DL", "REML + HK")) {
        r <- fit1(dl1, methods[[m]][1], methods[[m]][2])
        res[[length(res) + 1]] <- cbind(data.table(program = p, scaling = sc, cohorts = set, k = nrow(dl1),
                                                   analysis = paste("without", drop), method = m), r)
      }
    }
  }
}
res <- rbindlist(res)
res[, supported := is.finite(hi) & hi < 1 & est < 1]
num <- c("est", "lo", "hi", "tau2", "I2", "pi_lo", "pi_hi"); res[, (num) := lapply(.SD, round, 3), .SDcols = num]
write_out(res, "55_meta_analysis_sensitivity")

## compact table for the manuscript: frozen, main cohorts
wide <- dcast(res[scaling == "frozen" & cohorts == "main" & analysis == "all cohorts"],
              program + k ~ method, value.var = c("est", "lo", "hi"))
write_out(wide, "55_meta_frozen_main_wide")
print(res[scaling == "frozen" & cohorts == "main" & analysis == "all cohorts", .(program, k, method, est, lo, hi, I2, supported)])
cat("\nCheck: the DL rows must equal S8b (e.g. hypoxia 0.808, 0.727-0.897). With k = 2 the Hartung-Knapp\n",
    "interval uses a t distribution with 1 df and is expected to be very wide; report it as such.\n")
