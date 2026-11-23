###############################################################################
# 52_reproduce_published.R - STEP 2: REPRODUCTION GATE
# Recomputes every published frozen-panel validation result (Supplementary
# Table S8, both scalings) from the pseudobulk files with the helper functions,
# and compares with the published values. Later scripts reuse these helpers,
# so they are trustworthy only if this gate passes.
# Outputs: 52_validation_recomputed.csv, 52_donor_predictions.rds,
#          52_gate_comparison.csv, 52_gate_status.txt
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision()

panels <- load_panels(); partition <- load_partition()
cmps <- unique(program_compartment(panels))
cohorts <- c(DISCOVERY, VALIDATION, DESCRIPTIVE)
msg("Loading pseudobulk and computing log-CPM ...")
lc <- build_lc_cache(cohorts, cmps)
for (co in cohorts) for (cm in cmps) if (!is.null(lc[[co]][[cm]])) msg(sprintf("%-10s %-12s %3d donors", co, cm, ncol(lc[[co]][[cm]])))
dstat <- all_disc_stats(panels, partition, lc)

## Within-cohort standardisation (secondary analysis): z with the cohort's own statistics
within_cohort <- function(prog, co) {
  cmp <- program_compartment(panels)[[prog]]; L <- lc[[co]][[cmp]]; key <- key_col(co)
  et <- eval_table(prog, partition, dstat); g <- intersect(et[[key]], rownames(L))
  s <- zstats(L[g, , drop = FALSE]); ok <- s$sd > 1e-8
  y <- colMeans(zapply(L[g[ok], , drop = FALSE], lapply(s, `[`, ok)))
  pn <- panels[[prog]]; if (!all(pn[[key]] %in% rownames(L))) return(NULL)
  ps <- zstats(L[pn[[key]], , drop = FALSE])
  yhat <- colSums(zapply(L[pn[[key]], , drop = FALSE], ps) * pn$coef)
  list(y = y, yhat = yhat)
}

res <- list(); preds <- list()
for (p in names(panels)) for (co in c(VALIDATION, DESCRIPTIVE)) {
  e <- eval_program_cohort(p, co, panels, partition, dstat, lc)
  if (is.null(e)) next
  row <- data.table(program = p, cohort = co, donors = e$n, coverage = round(e$coverage, 3))
  if (is.null(e$yhat) || e$coverage < MIN_EVAL_COVER || e$n < 3) {
    res[[length(res) + 1]] <- cbind(row, scaling = "frozen", ratio = NA, lo = NA, hi = NA, rho = NA, p_boot = NA,
                                    status = "Not evaluable"); next
  }
  preds[[paste(p, co)]] <- data.table(program = p, cohort = co, donor = names(e$y), y = e$y, yhat = e$yhat)
  set.seed(SEED)
  b <- boot_ratio(e$y, e$yhat)
  r <- mae_ratio(e$y, e$yhat)
  st <- if (e$n < MIN_DONORS) "Descriptive only (<10 donors)" else if (r < 1 && b$ci[2] < 1) "Supported" else if (r < 1) "Partially supported" else "Not supported"
  res[[length(res) + 1]] <- cbind(row, scaling = "frozen", ratio = r, lo = b$ci[1], hi = b$ci[2],
                                  rho = cor(e$y, e$yhat, method = "spearman"), p_boot = b$p_one_sided, status = st)
  w <- within_cohort(p, co)
  if (!is.null(w)) {
    set.seed(SEED); bw <- boot_ratio(w$y, w$yhat); rw <- mae_ratio(w$y, w$yhat)
    stw <- if (e$n < MIN_DONORS) "Descriptive only (<10 donors)" else if (rw < 1 && bw$ci[2] < 1) "Supported" else if (rw < 1) "Partially supported" else "Not supported"
    res[[length(res) + 1]] <- cbind(row, scaling = "within", ratio = rw, lo = bw$ci[1], hi = bw$ci[2],
                                    rho = cor(w$y, w$yhat, method = "spearman"), p_boot = bw$p_one_sided, status = stw)
  }
}
res <- rbindlist(res, fill = TRUE)
write_out(res, "52_validation_recomputed")
saveRDS(list(preds = rbindlist(preds), lc_note = "frozen-scale targets and predictions per donor"),
        file.path(OUT_DIR, "52_donor_predictions.rds"))

## ---- Compare with published S8 ------------------------------------------------
s8 <- read_supp("S8_validation_by_cohort")
setnames_re(s8, c("^Program code$", "^Program$", "^Cohort$", "^Scaling$", "^Donors$", "^MAE ratio$",
                  "lower$", "upper$", "^Spearman", "^Status$"),
            c("program", "program_name", "cohort", "scaling_lab", "donors_pub", "ratio_pub", "lo_pub", "hi_pub", "rho_pub", "status_pub"))
s8[, scaling := ifelse(grepl("^Frozen", scaling_lab), "frozen", "within")]
num <- c("donors_pub", "ratio_pub", "lo_pub", "hi_pub", "rho_pub"); s8[, (num) := lapply(.SD, as.numeric), .SDcols = num]
cmpt <- merge(s8[, .(program, cohort, scaling, donors_pub, ratio_pub, lo_pub, hi_pub, rho_pub, status_pub)],
              res[, .(program, cohort, scaling, donors, ratio, lo, hi, rho, status)],
              by = c("program", "cohort", "scaling"), all = TRUE)
cmpt[, `:=`(d_ratio = round(ratio - ratio_pub, 4), d_lo = round(lo - lo_pub, 3), d_hi = round(hi - hi_pub, 3),
            donors_match = donors == donors_pub)]
cmpt[, pass := (is.na(ratio_pub) & (is.na(ratio) | status == "Not evaluable")) |
              (!is.na(ratio_pub) & !is.na(ratio) & abs(d_ratio) <= GATE_TOL_RATIO &
                 abs(d_lo) <= GATE_TOL_CI & abs(d_hi) <= GATE_TOL_CI & donors_match %in% TRUE)]
write_out(cmpt, "52_gate_comparison")

## ---- Pooled (main analysis, frozen) vs S8b ---------------------------------------
s8b <- read_supp("S8b_meta_analysis")
setnames_re(s8b, c("^Program code$", "^Analysis$", "^Scaling$", "^Pooled MAE ratio$", "^95% CI, lower$", "^95% CI, upper$", "^I.*ratio$"),
            c("program", "analysis", "scaling", "pooled_pub", "plo_pub", "phi_pub", "I2_pub"))
pool <- res[cohort %in% VALIDATION & donors >= MIN_DONORS & !is.na(ratio) & is.finite(lo) & lo > 0,
            { d <- dl_pool(log(ratio), se_from_ci(lo, hi))
              .(k = .N, pooled = exp(d$est), plo = exp(d$lo), phi = exp(d$hi), I2 = d$I2) },
            by = .(program, scaling)][k >= 2]
pool[, scaling_lab := ifelse(scaling == "frozen", "Frozen", "Within-cohort")]
pc <- merge(pool, s8b[grepl("^Main", analysis), .(program, scaling_lab = scaling, pooled_pub = as.numeric(pooled_pub),
                                                  plo_pub = as.numeric(plo_pub), phi_pub = as.numeric(phi_pub))],
            by = c("program", "scaling_lab"), all.x = TRUE)
pc[, d_pooled := round(pooled - pooled_pub, 3)]
write_out(pc, "52_pooled_comparison")

## ---- Verdict ---------------------------------------------------------------------
fails <- cmpt[pass == FALSE | is.na(pass)]
cat("\n================ REPRODUCTION GATE ================\n")
print(cmpt[, .(program, cohort, scaling, donors_pub, donors, ratio_pub, ratio = round(ratio, 3), d_ratio, d_lo, d_hi, pass)], nrows = 200)
cat("\nPooled (frozen + within) vs S8b:\n"); print(pc[, .(program, scaling_lab, k, pooled = round(pooled, 3), pooled_pub, d_pooled)])
verdict <- if (nrow(fails) == 0) "PASS" else "FAIL"
writeLines(c(verdict, paste("checked", nrow(cmpt), "rows;", nrow(fails), "failures"), format(Sys.time())), gate_file())
cat("\nGATE:", verdict, "-", nrow(fails), "failing rows (see 52_gate_comparison.csv)\n")
if (verdict == "FAIL") cat(
"Typical causes: (1) wrong pseudobulk file or donors not filtered to >=20 cells (check donors_match);\n",
"(2) discovery rownames are symbols instead of matrix row IDs; (3) duplicated symbols handled differently\n",
"in a validation cohort (sum vs make.unique); (4) CI-only differences larger than tolerance usually reflect\n",
"the bootstrap seed; if every ratio matches and only CIs differ slightly, raise GATE_TOL_CI and note it.\n")
