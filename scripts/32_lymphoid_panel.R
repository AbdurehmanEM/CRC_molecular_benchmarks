
# 32_lymphoid_panel.R  (protocol v1.8, section 14.4) - compact panel for LYM_IL2_STAT5, DISCOVERY only.
# Helpers and settings copied verbatim from 10_panel_development.R so results are comparable.
# Writes frozen_panel_LYM_v1.rds; frozen_panel_v2.rds is NOT touched. No validation data read here.
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(glmnet); library(digest) })
SEED <- 20261008; MIN_CELLS <- 20; MAX_GENES <- 4; ALPHA <- 0.5; K <- 5; REPS <- 10
pb_dir <- "04_processed_data/GSE178341"; prog_dir <- "04_processed_data/programs"
res_dir <- "06_results/panel"; frz_dir <- "04_processed_data/frozen_panel"

programs <- fread(file.path(prog_dir, "program_definitions_lymphoid.csv"))[retained == TRUE]
part <- fread(file.path(prog_dir, "gene_partition_lymphoid.tsv"))
dm   <- fread(file.path(pb_dir, "donor_metadata_tumour.csv"))
conventional <- list(LYM_IL2_STAT5 = fread(file.path(prog_dir, "conventional_markers_lymphoid.csv"))$gene_symbol)
cat("Programs modelled:", paste(programs$program, collapse = ", "), "\n")

zfit   <- function(M) { mu <- colMeans(M); s <- apply(M, 2, sd); s[!is.finite(s) | s == 0] <- 1; list(mu = mu, sd = s) }
zapply <- function(M, z) sweep(sweep(M, 2, z$mu, "-"), 2, z$sd, "/")
metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred); o <- obs[ok]; p <- pred[ok]
  if (length(o) < 3) return(data.table(MAE = NA_real_, RMSE = NA_real_, spearman = NA_real_,
                                       cal_intercept = NA_real_, cal_slope = NA_real_, n = length(o)))
  cal <- if (sd(p) > 1e-8) coef(lm(o ~ p)) else c(NA_real_, NA_real_)
  data.table(MAE = mean(abs(o - p)), RMSE = sqrt(mean((o - p)^2)),
             spearman = suppressWarnings(cor(o, p, method = "spearman")),
             cal_intercept = unname(cal[1]), cal_slope = unname(cal[2]), n = length(o)) }
calib_score <- function(s_tr, y_tr, s_te) { f <- lm(y_tr ~ s_tr); unname(coef(f)[1] + coef(f)[2] * s_te) }
fit_all <- function(Etr, Ptr, Ete, Pte, sym_pr, conv) {
  ze <- zfit(Etr); y_tr <- rowMeans(zapply(Etr, ze)); y_te <- rowMeans(zapply(Ete, ze))
  keep <- apply(Ptr, 2, sd) > 0; Ptr <- Ptr[, keep, drop = FALSE]; Pte <- Pte[, keep, drop = FALSE]
  zp <- zfit(Ptr); Xtr <- zapply(Ptr, zp); Xte <- zapply(Pte, zp)
  nf <- min(5, nrow(Xtr))
  cv <- suppressWarnings(cv.glmnet(Xtr, y_tr, alpha = ALPHA, nfolds = nf, standardize = FALSE, grouped = FALSE))
  ok <- which(cv$nzero >= 1 & cv$nzero <= MAX_GENES); if (!length(ok)) ok <- which(cv$nzero >= 1)[1]
  lam <- cv$lambda[ok][which.min(cv$cvm[ok])]
  b <- as.matrix(coef(cv, s = lam)); sel <- setdiff(rownames(b)[b[, 1] != 0], "(Intercept)")
  refit <- lm(y_tr ~ ., data = data.frame(y_tr = y_tr, Xtr[, sel, drop = FALSE], check.names = FALSE))
  pred <- list(panel_refit = as.numeric(predict(refit, newdata = data.frame(Xte[, sel, drop = FALSE], check.names = FALSE))),
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
       n_larger = cv$nzero[cv$lambda == cv$lambda.min], conv_genes = cg) }

perf <- list(); stab <- list(); frozen <- list(); folds_used <- list()
for (k in seq_len(nrow(programs))) {
  pg <- programs$program[k]; comp <- programs$compartment[k]
  pb <- readRDS(file.path(pb_dir, "pseudobulk_Lymphoid_CD8.rds"))
  donors <- names(pb$n_cells)[pb$n_cells >= MIN_CELLS]
  L <- edgeR::cpm(as.matrix(pb$counts[, donors]), log = TRUE, prior.count = 2)
  ev <- part[program == pg & role == "evaluation"]; pr <- part[program == pg & role == "predictor"]
  E <- t(L[ev$row_id, , drop = FALSE]); P <- t(L[pr$row_id, , drop = FALSE])
  sym_pr <- setNames(pr$gene_symbol, pr$row_id)
  chem <- dm$chemistry[match(donors, dm$PID)]; n <- length(donors)
  cat("\n==", pg, "(", comp, "):", n, "donors,", ncol(E), "evaluation,", ncol(P), "predictor genes ==\n")
  set.seed(SEED + 7 + k)      # existing programs used SEED+1 .. SEED+7
  preds <- list(); sel_count <- setNames(integer(ncol(P)), colnames(P))
  for (r in seq_len(REPS)) {
    fold <- integer(n)
    for (lv in unique(chem)) { ix <- which(chem == lv); fold[ix] <- sample(rep(seq_len(K), length.out = length(ix))) }
    folds_used[[paste(pg, r)]] <- data.table(program = pg, rep = r, PID = donors, fold = fold)
    for (fo in seq_len(K)) {
      te <- which(fold == fo); tr <- which(fold != fo)
      out <- fit_all(E[tr, , drop = FALSE], P[tr, , drop = FALSE], E[te, , drop = FALSE], P[te, , drop = FALSE],
                     sym_pr, conventional[[pg]])
      sel_count[out$sel] <- sel_count[out$sel] + 1L
      for (m in names(out$pred))
        preds[[length(preds) + 1]] <- data.table(rep = r, PID = donors[te], method = m,
                                                 obs = out$y_te, pred = out$pred[[m]]) } }
  preds <- rbindlist(preds)
  per_rep <- preds[, metrics(obs, pred), by = .(method, rep)][is.finite(MAE)]
  summ <- per_rep[, .(MAE = mean(MAE), MAE_rep_min = min(MAE), MAE_rep_max = max(MAE), RMSE = mean(RMSE),
                      spearman = mean(spearman, na.rm = TRUE), cal_slope = mean(cal_slope, na.rm = TRUE),
                      cal_intercept = mean(cal_intercept, na.rm = TRUE)), by = method]
  summ[, MAE_vs_mean := round(MAE / summ[method == "mean", MAE], 3)]
  summ[method == "mean", `:=`(cal_slope = NA_real_, cal_intercept = NA_real_, spearman = NA_real_)]
  perf[[pg]] <- cbind(program = pg, n_donors = n, summ); print(summ)
  stab[[pg]] <- data.table(program = pg, row_id = names(sel_count), gene_symbol = sym_pr[names(sel_count)],
                           selection_freq = round(sel_count / (K * REPS), 3))[order(-selection_freq)]
  fin <- fit_all(E, P, E, P, sym_pr, conventional[[pg]])
  rc <- coef(fin$refit); names(rc) <- gsub("`", "", names(rc)); nz <- setdiff(names(rc), "(Intercept)")
  frozen[[pg]] <- list(program = pg, compartment = comp, n_donors = n, donors = donors,
    panel = data.table(row_id = nz, gene_symbol = sym_pr[nz], coefficient = unname(rc[nz]),
                       shrunk_coefficient = fin$b[nz, 1], pred_mean = fin$zp$mu[nz], pred_sd = fin$zp$sd[nz]),
    intercept = unname(rc["(Intercept)"]), lambda = fin$lam,
    evaluation_genes = data.table(row_id = ev$row_id, gene_symbol = ev$gene_symbol,
                                  eval_mean = fin$ze$mu[ev$row_id], eval_sd = fin$ze$sd[ev$row_id]),
    conventional_genes_used = fin$conv_genes)
  cat("Final panel:", paste(sym_pr[nz], collapse = ", "), "\n") }

perf <- rbindlist(perf); stab <- rbindlist(stab)
fwrite(perf, file.path(res_dir, "cv_performance_discovery_lymphoid.csv"))
fwrite(stab, file.path(res_dir, "selection_stability_lymphoid.csv"))
fwrite(rbindlist(folds_used), file.path(res_dir, "cv_folds_lymphoid.csv"))
final_tab <- rbindlist(lapply(frozen, function(x) cbind(program = x$program, x$panel, intercept = x$intercept)))
fwrite(final_tab, file.path(res_dir, "final_panel_coefficients_lymphoid.csv"))
frozen_obj <- list(programs = frozen,
  normalisation = "log2-CPM per sample, edgeR::cpm(log = TRUE, prior.count = 2), library-size only",
  target_definition = "mean over evaluation genes of (log2-CPM - eval_mean) / eval_sd, discovery parameters",
  prediction = "intercept + sum(coefficient * (log2-CPM - pred_mean) / pred_sd)",
  gene_id_rule = "match by Ensembl ID (version stripped) where available, otherwise by gene symbol",
  missing_gene_rule = "if any panel gene is absent in a validation cohort, that program is reported as not evaluable (no imputation)",
  metrics = "primary MAE; secondary RMSE, Spearman, calibration slope/intercept; compared with discovery-mean baseline",
  min_cells = MIN_CELLS, max_genes = MAX_GENES, alpha = ALPHA,
  partition_hash = readLines(file.path(prog_dir, "gene_partition_lymphoid.sha"))[2],
  created = as.character(Sys.time()), protocol = "v1.8",
  note = "Additive lymphoid compartment. Independent of frozen_panel_v2.rds, which is unchanged.")
fp <- file.path(frz_dir, "frozen_panel_LYM_v1.rds"); saveRDS(frozen_obj, fp)
sha <- digest(file = fp, algo = "sha256")
writeLines(c(format(Sys.time()), sha), file.path(frz_dir, "frozen_panel_LYM_v1.sha256"))
stopifnot(digest(file = file.path(frz_dir, "frozen_panel_v2.rds"), algo = "sha256") ==
          "a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb")
cat("\nfrozen_panel_v2.rds verified unchanged\n")
cat("\n=== Final lymphoid panel ===\n"); print(final_tab[, .(program, gene_symbol, coefficient = round(coefficient, 3))])
cat("\nLymphoid frozen panel sha256:", sha, "\n")
cat("\n=== CV performance (discovery, lymphoid) ===\n")
print(perf[, .(program, n_donors, method, MAE = round(MAE, 3), MAE_vs_mean,
               spearman = round(spearman, 2), cal_slope = round(cal_slope, 2))])
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_32.txt"))
 
