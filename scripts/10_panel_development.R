# 10_panel_development.R
# Protocol v0.7, step 11: compact panel per retained program, DISCOVERY (GSE178341) only.
# - Normalisation: sample-wise log2-CPM (prior count 2) -> transferable to new samples
# - Target: mean z-score of EVALUATION genes (z parameters from training donors only)
# - Panel: elastic net (alpha 0.5) on PREDICTOR genes, at most 4 genes
# - Outer CV: 5-fold x 10 repeats, folds balanced by 10x chemistry; all scaling,
#   selection and tuning inside training folds
# - Baselines: training mean; equal-weight top-4 score; conventional markers (predictor
#   pool only); larger uncapped elastic net
# - Final model refit on all eligible donors and FROZEN (saved with hash)
# No validation-cohort data are read here.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(glmnet); library(digest) })
SEED <- 20261008; MIN_CELLS <- 20; MAX_GENES <- 4; ALPHA <- 0.5; K <- 5; REPS <- 10
pb_dir   <- "04_processed_data/GSE178341"
prog_dir <- "04_processed_data/programs"
res_dir  <- "06_results/panel"
frz_dir  <- "04_processed_data/frozen_panel"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(frz_dir, recursive = TRUE, showWarnings = FALSE)

programs <- fread(file.path(prog_dir, "program_definitions.csv"))[retained == TRUE]
part     <- fread(file.path(prog_dir, "gene_partition.tsv"))
dm       <- fread(file.path(pb_dir, "donor_metadata_tumour.csv"))
cat("Programs modelled:", paste(programs$program, collapse = ", "), "\n")

# Conventional markers (pre-specified); only those in the program's PREDICTOR pool are used
conventional <- list(
  EPI_HYPOXIA   = c("VEGFA", "CA9", "SLC2A1", "PGK1", "LDHA", "BNIP3", "NDRG1", "ADM"),
  EPI_APOPTOSIS = c("BAX", "BCL2", "CASP3", "FAS", "BAK1", "CASP8", "BCL2L1", "BID"),
  STR_TGFB      = c("SERPINE1", "TGFBI", "SMAD7", "SKIL", "LTBP2", "TGFB1", "SMAD3"),
  MYE_INFLAM    = c("IL1B", "TNF", "CXCL8", "IL6", "CCL2", "NFKBIA", "PTGS2"),
  MYE_IFNG      = c("CXCL9", "CXCL10", "STAT1", "IRF1", "GBP1", "HLA-DRA", "IDO1"),
  END_ANGIO     = c("VEGFA", "KDR", "ANGPT2", "POSTN", "PDGFA", "VCAN"),
  END_TNFA_NFKB = c("SELE", "VCAM1", "ICAM1", "NFKBIA", "CXCL2", "IER3", "SELP"))

# ---------------------------------------------------------------- helpers
zfit   <- function(M) { mu <- colMeans(M); s <- apply(M, 2, sd); s[!is.finite(s) | s == 0] <- 1; list(mu = mu, sd = s) }
zapply <- function(M, z) sweep(sweep(M, 2, z$mu, "-"), 2, z$sd, "/")
metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  o <- obs[ok]; p <- pred[ok]
  if (length(o) < 3) return(data.table(MAE = NA_real_, RMSE = NA_real_, spearman = NA_real_,
                                       cal_intercept = NA_real_, cal_slope = NA_real_, n = length(o)))
  cal <- if (sd(p) > 1e-8) coef(lm(o ~ p)) else c(NA_real_, NA_real_)
  data.table(MAE = mean(abs(o - p)), RMSE = sqrt(mean((o - p)^2)),
             spearman = suppressWarnings(cor(o, p, method = "spearman")),
             cal_intercept = unname(cal[1]), cal_slope = unname(cal[2]), n = length(o))
}
calib_score <- function(s_tr, y_tr, s_te) { f <- lm(y_tr ~ s_tr); unname(coef(f)[1] + coef(f)[2] * s_te) }

fit_all <- function(Etr, Ptr, Ete, Pte, sym_pr, conv) {
  ze <- zfit(Etr); y_tr <- rowMeans(zapply(Etr, ze)); y_te <- rowMeans(zapply(Ete, ze))
  keep <- apply(Ptr, 2, sd) > 0
  Ptr <- Ptr[, keep, drop = FALSE]; Pte <- Pte[, keep, drop = FALSE]
  zp <- zfit(Ptr); Xtr <- zapply(Ptr, zp); Xte <- zapply(Pte, zp)
  nf <- min(5, nrow(Xtr))
  cv <- suppressWarnings(cv.glmnet(Xtr, y_tr, alpha = ALPHA, nfolds = nf, standardize = FALSE, grouped = FALSE))
  ok <- which(cv$nzero >= 1 & cv$nzero <= MAX_GENES)
  if (!length(ok)) ok <- which(cv$nzero >= 1)[1]
  lam <- cv$lambda[ok][which.min(cv$cvm[ok])]
  b <- as.matrix(coef(cv, s = lam)); sel <- setdiff(rownames(b)[b[, 1] != 0], "(Intercept)")
  refit <- lm(y_tr ~ ., data = data.frame(y_tr = y_tr, Xtr[, sel, drop = FALSE], check.names = FALSE))
  pred <- list(
    panel_refit = as.numeric(predict(refit, newdata = data.frame(Xte[, sel, drop = FALSE], check.names = FALSE))),
    panel_shrunk = as.numeric(predict(cv, Xte, s = lam)),
    larger_model = as.numeric(predict(cv, Xte, s = "lambda.min")),
    mean = rep(mean(y_tr), nrow(Xte)))
  r <- suppressWarnings(cor(Xtr, y_tr)[, 1]); r[!is.finite(r)] <- 0
  top <- names(sort(abs(r), decreasing = TRUE))[1:min(MAX_GENES, length(r))]
  pred$equal_weight <- calib_score(as.numeric(Xtr[, top, drop = FALSE] %*% sign(r[top])) / length(top), y_tr,
                                   as.numeric(Xte[, top, drop = FALSE] %*% sign(r[top])) / length(top))
  cg <- intersect(names(sym_pr)[sym_pr %in% conv], colnames(Xtr))
  pred$conventional <- if (length(cg)) calib_score(rowMeans(Xtr[, cg, drop = FALSE]), y_tr,
                                                   rowMeans(Xte[, cg, drop = FALSE])) else rep(NA_real_, nrow(Xte))
  list(y_te = y_te, pred = pred, sel = sel, b = b, refit = refit, lam = lam, ze = ze, zp = zp,
       n_larger = cv$nzero[cv$lambda == cv$lambda.min], conv_genes = cg)
}

# ---------------------------------------------------------------- per program
perf <- list(); stab <- list(); frozen <- list(); folds_used <- list()
for (k in seq_len(nrow(programs))) {
  pg <- programs$program[k]; comp <- programs$compartment[k]
  pb <- readRDS(file.path(pb_dir, paste0("pseudobulk_", comp, ".rds")))
  donors <- names(pb$n_cells)[pb$n_cells >= MIN_CELLS]
  L <- cpm(as.matrix(pb$counts[, donors]), log = TRUE, prior.count = 2)
  ev <- part[program == pg & role == "evaluation"]; pr <- part[program == pg & role == "predictor"]
  E <- t(L[ev$row_id, , drop = FALSE]); P <- t(L[pr$row_id, , drop = FALSE])
  sym_pr <- setNames(pr$gene_symbol, pr$row_id)
  chem <- dm$chemistry[match(donors, dm$PID)]
  n <- length(donors)
  cat("\n==", pg, "(", comp, "):", n, "donors,", ncol(E), "evaluation,", ncol(P), "predictor genes ==\n")

  set.seed(SEED + k)
  preds <- list(); sel_count <- setNames(integer(ncol(P)), colnames(P))
  for (r in seq_len(REPS)) {
    fold <- integer(n)
    for (lv in unique(chem)) { ix <- which(chem == lv); fold[ix] <- sample(rep(seq_len(K), length.out = length(ix))) }
    folds_used[[paste(pg, r)]] <- data.table(program = pg, rep = r, PID = donors, fold = fold)
    for (f in seq_len(K)) {
      te <- which(fold == f); tr <- which(fold != f)
      out <- fit_all(E[tr, , drop = FALSE], P[tr, , drop = FALSE], E[te, , drop = FALSE], P[te, , drop = FALSE],
                     sym_pr, conventional[[pg]])
      sel_count[out$sel] <- sel_count[out$sel] + 1L
      for (m in names(out$pred))
        preds[[length(preds) + 1]] <- data.table(rep = r, PID = donors[te], method = m,
                                                 obs = out$y_te, pred = out$pred[[m]])
    }
  }
  preds <- rbindlist(preds)
  per_rep <- preds[, metrics(obs, pred), by = .(method, rep)]
  per_rep <- per_rep[is.finite(MAE)]
  summ <- per_rep[, .(MAE = mean(MAE), MAE_rep_min = min(MAE), MAE_rep_max = max(MAE),
                      RMSE = mean(RMSE), spearman = mean(spearman, na.rm = TRUE),
                      cal_slope = mean(cal_slope, na.rm = TRUE), cal_intercept = mean(cal_intercept, na.rm = TRUE)),
                  by = method]
  summ[, MAE_vs_mean := round(MAE / summ[method == "mean", MAE], 3)]
  summ[method == "mean", `:=`(cal_slope = NA_real_, cal_intercept = NA_real_, spearman = NA_real_)]
  perf[[pg]] <- cbind(program = pg, n_donors = n, summ)
  print(summ)

  stab[[pg]] <- data.table(program = pg, row_id = names(sel_count), gene_symbol = sym_pr[names(sel_count)],
                           selection_freq = round(sel_count / (K * REPS), 3))[order(-selection_freq)]

  # ---------------- final fit on all donors (frozen)
  fin <- fit_all(E, P, E, P, sym_pr, conventional[[pg]])
  rc <- coef(fin$refit); names(rc) <- gsub("`", "", names(rc)); nz <- setdiff(names(rc), "(Intercept)")
  frozen[[pg]] <- list(
    program = pg, compartment = comp, n_donors = n, donors = donors,
    panel = data.table(row_id = nz, gene_symbol = sym_pr[nz], coefficient = unname(rc[nz]), shrunk_coefficient = fin$b[nz, 1],
                       pred_mean = fin$zp$mu[nz], pred_sd = fin$zp$sd[nz]),
    intercept = unname(rc["(Intercept)"]), lambda = fin$lam,
    evaluation_genes = data.table(row_id = ev$row_id, gene_symbol = ev$gene_symbol,
                                  eval_mean = fin$ze$mu[ev$row_id], eval_sd = fin$ze$sd[ev$row_id]),
    conventional_genes_used = fin$conv_genes)
  cat("Final panel:", paste(sym_pr[nz], collapse = ", "), "\n")
}

# ---------------------------------------------------------------- save
perf <- rbindlist(perf); stab <- rbindlist(stab)
fwrite(perf, file.path(res_dir, "cv_performance_discovery_v2.csv"))
fwrite(stab, file.path(res_dir, "selection_stability_v2.csv"))
fwrite(rbindlist(folds_used), file.path(res_dir, "cv_folds_v2.csv"))
final_tab <- rbindlist(lapply(frozen, function(x) cbind(program = x$program, x$panel, intercept = x$intercept)))
fwrite(final_tab, file.path(res_dir, "final_panel_coefficients_v2.csv"))

frozen_obj <- list(
  programs = frozen,
  normalisation = "log2-CPM per sample, edgeR::cpm(log = TRUE, prior.count = 2), library-size only",
  target_definition = "mean over evaluation genes of (log2-CPM - eval_mean) / eval_sd, discovery parameters",
  prediction = "intercept + sum(coefficient * (log2-CPM - pred_mean) / pred_sd)",
  gene_id_rule = "match by Ensembl ID (version stripped) where available, otherwise by gene symbol",
  missing_gene_rule = "if any panel gene is absent in a validation cohort, that program is reported as not evaluable (no imputation)",
  metrics = "primary MAE; secondary RMSE, Spearman, calibration slope/intercept; compared with discovery-mean baseline",
  min_cells = MIN_CELLS, max_genes = MAX_GENES, alpha = ALPHA,
  partition_hash = readLines(file.path(prog_dir, "gene_partition.sha"))[2],
  created = as.character(Sys.time()), protocol = "v0.8", supersedes = "frozen_panel_v1 (over-shrunk coefficients)")
f <- file.path(frz_dir, "frozen_panel_v2.rds")
saveRDS(frozen_obj, f)
sha <- digest(file = f, algo = "sha256")
writeLines(c(format(Sys.time()), sha), file.path(frz_dir, "frozen_panel_v2.sha256"))
cat("\n=== Final panels ===\n"); print(final_tab[, .(program, gene_symbol, coefficient = round(coefficient, 3))])
cat("\nFrozen panel sha256:", sha, "\n")
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_10.txt"))
cat("\n=== CV performance (discovery, v2) ===\n"); print(perf[, .(program, n_donors, method, MAE = round(MAE, 3), MAE_vs_mean, spearman = round(spearman, 2), cal_slope = round(cal_slope, 2))], nrows = 60)

