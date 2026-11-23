###############################################################################
# 53_null_random_panels.R - Reviewer point 7 (predictor-target dependence)
#
# Question: is the frozen panels' accuracy greater than what random genes,
# matched for expression level and variability, achieve against the same
# targets? And does the frozen panel track arbitrary (negative-control)
# gene-set scores, which would indicate shared technical/batch variation?
#
# Three nulls, all on the frozen discovery scale and the same validation rules:
#  A. Fixed random panels: same size as the frozen panel, genes drawn from
#     non-program genes matched gene-by-gene on discovery mean log-CPM decile
#     x SD tercile; OLS on discovery donors; scored in validation cohorts.
#  B. Selection null: a random matched POOL the size of the program's
#     predictor pool, then the same relaxed elastic net (1-4 genes, OLS refit).
#     This is the fair null for the whole selection procedure.
#  C. Negative-control targets: random gene sets matched to the evaluation
#     genes; Spearman correlation of the FROZEN panel prediction with them.
# Empirical p = (1 + #null <= observed) / (1 + B)   (one-sided, lower is better
# for MAE ratios; for C, p = (1 + #null >= observed rho) / (1 + B)).
# Post hoc (added in revision); label as such.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision(); require_gate()
suppressPackageStartupMessages(library(glmnet))

panels <- load_panels(); partition <- load_partition()
pc <- program_compartment(panels); cmps <- unique(pc)
lc <- build_lc_cache(c(DISCOVERY, VALIDATION), cmps)
dstat <- all_disc_stats(panels, partition, lc)
symmap <- load_symbol_map(partition, lc[[DISCOVERY]][[cmps[1]]])

s2b <- tryCatch(read_supp("S2b_discovery_donors", "Patient"), error = function(e) NULL)
chem <- if (!is.null(s2b)) setNames(s2b$Chemistry, s2b$Patient) else NULL

all_program_rows <- unique(partition$row_id)
all_program_syms <- unique(partition$symbol)
B_A <- N_NULL_FIXED; B_B <- N_NULL_SELECT; B_C <- N_NULL_TARGET

pub <- fread(file.path(OUT_DIR, "52_validation_recomputed.csv"))[scaling == "frozen"]
summ <- list(); draws <- list()

for (p in names(panels)) {
  cmp <- pc[[p]]; L <- lc[[DISCOVERY]][[cmp]]
  msg("Program", p, "(", cmp, ", ", ncol(L), "discovery donors )")
  et <- eval_table(p, partition, dstat)
  y_disc <- target_score(L, et, "row_id")$y
  real <- panels[[p]]
  ## evaluable validation cohorts for this program (as published)
  vco <- pub[program == p & donors >= MIN_DONORS & !is.na(ratio), cohort]
  if (length(vco) == 0) { msg("  no evaluable validation cohort; nulls computed for discovery CV only") }
  ## universe: expressed (log2CPM >= 0 in >= 25% donors), not in any program, symbol present in every evaluable cohort
  expressed <- rownames(L)[rowMeans(L >= 0) >= 0.25]
  uni <- setdiff(expressed, all_program_rows)
  uni <- uni[!is.na(symmap[uni]) & !(symmap[uni] %in% all_program_syms)]
  for (co in vco) uni <- uni[symmap[uni] %in% rownames(lc[[co]][[cmp]])]
  if (length(uni) < 200) { msg("  WARNING: only", length(uni), "universe genes (feature map missing?); skipping"); next }
  bins <- expr_bins(L, unique(c(uni, real$row_id, et$row_id, partition[program == p, row_id])))
  strata <- if (!is.null(chem)) chem[colnames(L)] else NULL

  score_val <- function(panel) sapply(vco, function(co) {
    e <- eval_program_cohort(p, co, panels, partition, dstat, lc, panel = panel)
    if (is.null(e$yhat)) NA_real_ else mae_ratio(e$y, e$yhat)
  })

  ## observed
  set.seed(SEED)
  obs_cv  <- cv_panel(L, et$row_id, NULL, select = FALSE, fixed_genes = real$row_id, repeats = 2, strata = strata)$ratio
  obs_val <- score_val(real)

  ## ---- Null A: fixed matched random panels
  set.seed(SEED + 1)
  nullA <- t(replicate(B_A, {
    g <- unique(draw_matched(real$row_id, bins, uni))
    pn <- fit_fixed_panel(L, g, y_disc, symmap)
    c(cv = cv_panel(L, et$row_id, NULL, select = FALSE, fixed_genes = g, repeats = 2, strata = strata)$ratio,
      score_val(pn))
  }))
  ## ---- Null B: selection from a matched random pool
  pool_real <- partition[program == p & role == "predictor", row_id]
  pool_real <- intersect(pool_real, rownames(L))
  set.seed(SEED + 2)
  nullB <- t(replicate(B_B, {
    pool <- unique(draw_matched(pool_real, bins, uni))
    s <- zstats(L[pool, , drop = FALSE]); ok <- s$sd > 1e-8; pool <- pool[ok]
    X <- t(zapply(L[pool, , drop = FALSE], lapply(s, `[`, ok)))
    sel <- select_relaxed_en(X, y_disc)
    score_val(fit_fixed_panel(L, sel, y_disc, symmap))
  }))
  if (length(vco) == 1) nullB <- matrix(nullB, ncol = 1, dimnames = list(NULL, vco))
  ## ---- Null C: negative-control targets (frozen panel vs random target sets)
  set.seed(SEED + 3)
  nullC <- list()
  obsC <- c()
  for (co in vco) {
    e <- eval_program_cohort(p, co, panels, partition, dstat, lc)
    obsC[co] <- cor(e$yhat, e$y, method = "spearman")
    Lv <- lc[[co]][[cmp]]
    nullC[[co]] <- replicate(B_C, {
      g <- unique(draw_matched(et$row_id[et$row_id %in% names(bins)], bins, uni))
      st <- data.table(row_id = g, symbol = unname(symmap[g]), mu = rowMeans(L[g, , drop = FALSE]),
                       sd = apply(L[g, , drop = FALSE], 1, sd))
      cor(e$yhat, target_score(Lv, st, "symbol")$y, method = "spearman")
    })
  }

  ## ---- summaries
  emp_low  <- function(obs, nul) (1 + sum(nul <= obs, na.rm = TRUE)) / (1 + sum(is.finite(nul)))
  emp_high <- function(obs, nul) (1 + sum(nul >= obs, na.rm = TRUE)) / (1 + sum(is.finite(nul)))
  q <- function(x) quantile(x, c(0.5, 0.025, 0.975), na.rm = TRUE)
  add <- function(null, setting, obs, nul, p_emp) summ[[length(summ) + 1]] <<- data.table(
    program = p, null = null, setting = setting, observed = obs, null_median = q(nul)[1],
    null_q025 = q(nul)[2], null_q975 = q(nul)[3], n_null = sum(is.finite(nul)), p_empirical = p_emp)
  add("A fixed random panels", "discovery CV (fixed genes, OLS)", obs_cv, nullA[, "cv"], emp_low(obs_cv, nullA[, "cv"]))
  for (co in vco) {
    add("A fixed random panels", paste("validation", co), obs_val[[co]], nullA[, co], emp_low(obs_val[[co]], nullA[, co]))
    add("B selection from random pool", paste("validation", co), obs_val[[co]], nullB[, co], emp_low(obs_val[[co]], nullB[, co]))
    add("C negative-control targets (Spearman rho)", paste("validation", co), obsC[[co]], nullC[[co]], emp_high(obsC[[co]], nullC[[co]]))
  }
  if (length(vco) >= 2) {
    gm <- function(M) apply(log(M[, vco, drop = FALSE]), 1, mean)
    o <- mean(log(unlist(obs_val[vco])))
    add("A fixed random panels", "mean log ratio over cohorts", exp(o), exp(gm(nullA)), emp_low(o, gm(nullA)))
    add("B selection from random pool", "mean log ratio over cohorts", exp(o), exp(gm(nullB)), emp_low(o, gm(nullB)))
  }
  draws[[p]] <- list(A = nullA, B = nullB, C = nullC, universe_size = length(uni))
}
summ <- rbindlist(summ)
num <- names(summ)[sapply(summ, is.numeric)]; summ[, (num) := lapply(.SD, signif, 4), .SDcols = num]
write_out(summ, "53_null_random_panels_summary")
saveRDS(draws, file.path(OUT_DIR, "53_null_random_panels_draws.rds"))
print(summ)
cat("\nInterpretation: a frozen panel adds information beyond expression-matched random genes if its\n",
    "observed MAE ratio lies in the lower tail of nulls A/B (small p). For null C, a small p means the\n",
    "panel correlates with its own program target more than with random matched gene sets (specificity);\n",
    "a high null-C median rho means part of the apparent accuracy reflects shared technical/global\n",
    "variation (library complexity, batch), which should then be reported as a limitation.\n")
