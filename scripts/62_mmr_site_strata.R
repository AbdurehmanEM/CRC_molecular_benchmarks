###############################################################################
# 62_mmr_site_strata.R - Reviewer point 26 (clinical generalisation)
#
# Discovery: repeats the published donor-grouped CV (relaxed elastic net with
# selection inside folds, 5-fold x 10, chemistry-balanced) and splits the
# out-of-fold errors by MMR status and tumour site (Supplementary Table S2b).
# For each subgroup: CV MAE ratio (out-of-fold panel error / training-mean
# error), donor-bootstrap 95% CI, and the ratio of subgroup ratios with CI.
# Validation (optional): frozen-panel MAE ratio by MMR status if
# VALIDATION_MMR_FILE is supplied. Subgroups with < 8 donors are not reported.
# Hypothesis-generating; no multiplicity correction.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision(); require_gate()
suppressPackageStartupMessages(library(glmnet))

panels <- load_panels(); partition <- load_partition()
pc <- program_compartment(panels); cmps <- unique(pc)
lc <- build_lc_cache(DISCOVERY, cmps)
s2b <- read_supp("S2b_discovery_donors", "Patient")
clin <- s2b[, .(donor = Patient, MMR = `MMR status`, site = `Tumour site`, chem = Chemistry)]
MIN_SUB <- 8

grp_ratio <- function(err, berr) mean(err) / mean(berr)
out <- list()
for (p in names(panels)) {
  L <- lc[[DISCOVERY]][[pc[[p]]]]
  ev <- intersect(partition[program == p & role == "evaluation", row_id], rownames(L))
  pool <- intersect(partition[program == p & role == "predictor", row_id], rownames(L))
  set.seed(SEED)
  cv <- cv_panel(L, ev, pool, select = TRUE, strata = clin$chem[match(colnames(L), clin$donor)])
  err <- rowMeans(abs(cv$oof - cv$ytrue)); berr <- rowMeans(abs(cv$base - cv$ytrue))
  d <- data.table(donor = colnames(L), err = err, berr = berr)[clin, on = "donor", nomatch = 0]
  msg(p, ": overall CV MAE ratio", round(grp_ratio(d$err, d$berr), 3), "(published S7 for comparison)")
  for (v in c("MMR", "site")) {
    lv <- d[, .N, by = v][N >= MIN_SUB][[v]]
    for (l in lv) {
      s <- d[get(v) == l]; r <- grp_ratio(s$err, s$berr)
      bt <- replicate(N_BOOT, { i <- sample.int(nrow(s), replace = TRUE); grp_ratio(s$err[i], s$berr[i]) })
      out[[length(out) + 1]] <- data.table(program = p, cohort = DISCOVERY, variable = v, level = l, n = nrow(s),
                                           ratio = r, lo = quantile(bt, .025), hi = quantile(bt, .975))
    }
    if (length(lv) == 2) {     # ratio of ratios (level1 / level2) with bootstrap CI
      a <- d[get(v) == lv[1]]; b <- d[get(v) == lv[2]]
      rr <- replicate(N_BOOT, { ia <- sample.int(nrow(a), replace = TRUE); ib <- sample.int(nrow(b), replace = TRUE)
                                grp_ratio(a$err[ia], a$berr[ia]) / grp_ratio(b$err[ib], b$berr[ib]) })
      out[[length(out) + 1]] <- data.table(program = p, cohort = DISCOVERY, variable = v, level = paste(lv[1], "/", lv[2]),
                                           n = nrow(d), ratio = grp_ratio(a$err, a$berr) / grp_ratio(b$err, b$berr),
                                           lo = quantile(rr, .025), hi = quantile(rr, .975))
    }
  }
}
## optional: validation cohorts by MMR
if (!is.na(VALIDATION_MMR_FILE) && file.exists(VALIDATION_MMR_FILE)) {
  mm <- fread(VALIDATION_MMR_FILE)
  pr <- readRDS(file.path(OUT_DIR, "52_donor_predictions.rds"))$preds[mm, on = .(cohort, donor), nomatch = 0]
  vv <- pr[, { r <- mae_ratio(y, yhat); b <- boot_ratio(y, yhat); .(n = .N, ratio = r, lo = b$ci[1], hi = b$ci[2]) },
           by = .(program, cohort, level = mmr)][n >= MIN_SUB]
  out[[length(out) + 1]] <- vv[, .(program, cohort, variable = "MMR", level, n, ratio, lo, hi)]
}
out <- rbindlist(out, fill = TRUE)
num <- c("ratio", "lo", "hi"); out[, (num) := lapply(.SD, round, 3), .SDcols = num]
write_out(out, "62_mmr_site_strata")
print(out)
