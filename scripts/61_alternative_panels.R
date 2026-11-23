###############################################################################
# 61_alternative_panels.R - Reviewer point 13 (panel gene stability)
#
# Builds post hoc alternative panels from the cross-validation selection
# frequencies (Supplementary Table S7b) and asks whether program-level
# performance survives replacing the frozen genes:
#   next-best : the k most frequently selected NON-frozen genes (k = panel size)
#   swap-<g>  : frozen panel with gene g replaced by the most frequent non-panel gene
# Each alternative is refitted by OLS on discovery donors (fixed genes, frozen
# scale), evaluated by discovery CV and in every evaluable validation cohort.
# EXPLORATORY by construction (genes chosen after validation results were known);
# it does not replace the frozen panels.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision(); require_gate()

panels <- load_panels(); partition <- load_partition()
pc <- program_compartment(panels); cmps <- unique(pc)
lc <- build_lc_cache(c(DISCOVERY, VALIDATION), cmps)
dstat <- all_disc_stats(panels, partition, lc)
symmap <- load_symbol_map(partition, lc[[DISCOVERY]][[cmps[1]]])
s7b <- read_supp("S7b_selection_stability")
setnames(s7b, c("Program code", "Matrix row ID", "Gene symbol", "Selection frequency", "In frozen panel"),
         c("program", "row_id", "symbol", "freq", "frozen"))
s7b[, freq := as.numeric(freq)]
symmap <- c(symmap, setNames(s7b$symbol, s7b$row_id)); symmap <- symmap[!duplicated(names(symmap))]
s2b <- tryCatch(read_supp("S2b_discovery_donors", "Patient"), error = function(e) NULL)
chem <- if (!is.null(s2b)) setNames(s2b$Chemistry, s2b$Patient) else NULL
pub <- fread(file.path(OUT_DIR, "52_validation_recomputed.csv"))[scaling == "frozen"]

out <- list()
for (p in names(panels)) {
  cmp <- pc[[p]]; L <- lc[[DISCOVERY]][[cmp]]; real <- panels[[p]]; k <- nrow(real)
  et <- eval_table(p, partition, dstat); y <- target_score(L, et, "row_id")$y
  vco <- pub[program == p & donors >= MIN_DONORS & !is.na(ratio), cohort]
  cand <- s7b[program == p & frozen != "Yes" & row_id %in% rownames(L)][order(-freq)]
  for (co in vco) cand <- cand[symbol %in% rownames(lc[[co]][[cmp]])]
  if (nrow(cand) == 0) next
  alts <- list(frozen = real$row_id, `next-best` = head(cand$row_id, k))
  for (g in real$row_id) alts[[paste0("swap ", symmap[[g]], " -> ", cand$symbol[1])]] <- c(setdiff(real$row_id, g), cand$row_id[1])
  strata <- if (!is.null(chem)) chem[colnames(L)] else NULL
  for (nm in names(alts)) {
    genes <- alts[[nm]]
    pn <- if (nm == "frozen") real else fit_fixed_panel(L, genes, y, symmap)
    set.seed(SEED)
    cvr <- cv_panel(L, et$row_id, NULL, select = FALSE, fixed_genes = genes, strata = strata)$ratio
    sc_disc <- panel_predict(L, pn, "row_id"); sc_frozen <- panel_predict(L, real, "row_id")
    row <- data.table(program = p, panel = nm, genes = paste(symmap[genes], collapse = ", "),
                      mean_selection_freq = round(mean(s7b[program == p][match(genes, row_id), freq], na.rm = TRUE), 2),
                      cv_ratio = round(cvr, 3), r_with_frozen_discovery = round(cor(sc_disc, sc_frozen), 3))
    for (co in vco) {
      e <- eval_program_cohort(p, co, panels, partition, dstat, lc, panel = pn)
      row[[paste0("ratio_", co)]] <- if (is.null(e$yhat)) NA else round(mae_ratio(e$y, e$yhat), 3)
      row[[paste0("rho_", co)]] <- if (is.null(e$yhat)) NA else round(cor(e$y, e$yhat, method = "spearman"), 3)
    }
    out[[length(out) + 1]] <- row
  }
}
out <- rbindlist(out, fill = TRUE)
write_out(out, "61_alternative_panels")
print(out)
cat("\nReport (point 13): for each program, whether next-best and single-swap panels keep the MAE ratio\n",
    "below 1 in discovery CV and in the validation cohorts where the frozen panel was supported.\n")
