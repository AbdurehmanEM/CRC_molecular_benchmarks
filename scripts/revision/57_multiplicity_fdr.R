###############################################################################
# 57_multiplicity_fdr.R - Reviewer point 15 (no multiplicity correction)
#
# Family 1: every evaluable program x validation cohort test of H0: MAE ratio >= 1
#           (frozen scaling, primary endpoint; cohorts with >= 10 donors).
# Family 2: the pooled (DerSimonian-Laird) program-level tests.
# p-values: one-sided donor-bootstrap p from script 52 if available
#           (p = (1 + #boot ratio >= 1)/(1 + B)); otherwise a normal approximation
#           on the log scale from the published CI. Adjustments: Benjamini-Hochberg
#           (FDR) and Holm (FWER). Reported as a sensitivity analysis; the primary
#           endpoint remains the prespecified per-cohort MAE ratio.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision()

f52 <- file.path(OUT_DIR, "52_validation_recomputed.csv")
if (file.exists(f52)) {
  d <- fread(f52)[scaling == "frozen" & cohort %in% VALIDATION & donors >= MIN_DONORS & !is.na(ratio)]
  d[, p := p_boot]; src <- "bootstrap (script 52)"
} else {
  s8 <- read_supp("S8_validation_by_cohort")
  setnames(s8, c("Program code", "Cohort", "Scaling", "Donors", "MAE ratio", "95% CI, lower", "95% CI, upper"),
           c("program", "cohort", "scaling", "donors", "ratio", "lo", "hi"))
  num <- c("donors", "ratio", "lo", "hi"); s8[, (num) := lapply(.SD, as.numeric), .SDcols = num]
  d <- s8[grepl("^Frozen", scaling) & cohort %in% VALIDATION & donors >= MIN_DONORS & is.finite(ratio)]
  d[, p := pnorm(log(ratio) / se_from_ci(lo, hi))]; src <- "normal approximation from published CI"
}
d[, `:=`(q_BH = p.adjust(p, "BH"), p_Holm = p.adjust(p, "holm"))]
d[, `:=`(supported_unadjusted = ratio < 1 & hi < 1, supported_FDR05 = ratio < 1 & q_BH < 0.05)]
fam1 <- d[, .(family = "per cohort", program, cohort, donors, ratio = round(ratio, 3), lo = round(lo, 3), hi = round(hi, 3),
              p = signif(p, 3), q_BH = signif(q_BH, 3), p_Holm = signif(p_Holm, 3),
              supported_unadjusted, supported_FDR05, p_source = src)]

s8b <- read_supp("S8b_meta_analysis")
setnames(s8b, c("Program code", "Analysis", "Scaling", "Cohorts (k)", "Pooled MAE ratio", "95% CI, lower", "95% CI, upper"),
         c("program", "analysis", "scaling", "k", "pooled", "plo", "phi"))
num <- c("k", "pooled", "plo", "phi"); s8b[, (num) := lapply(.SD, as.numeric), .SDcols = num]
pb <- s8b[grepl("^Main", analysis) & scaling == "Frozen" & is.finite(pooled)]
pb[, p := pnorm(log(pooled) / se_from_ci(plo, phi))]
pb[, `:=`(q_BH = p.adjust(p, "BH"), p_Holm = p.adjust(p, "holm"))]
fam2 <- pb[, .(family = "pooled", program, cohort = paste0("k=", k), donors = NA_integer_, ratio = pooled, lo = plo, hi = phi,
               p = signif(p, 3), q_BH = signif(q_BH, 3), p_Holm = signif(p_Holm, 3),
               supported_unadjusted = pooled < 1 & phi < 1, supported_FDR05 = pooled < 1 & q_BH < 0.05,
               p_source = "normal approximation from pooled CI")]
out <- rbind(fam1, fam2, fill = TRUE)
write_out(out, "57_multiplicity_fdr")
print(out)
cat("\nSentence (point 15): Of [n] per-cohort tests supported without adjustment, [m] remained supported at\n",
    "a Benjamini-Hochberg FDR of 5%; [list any change in support].\n")
