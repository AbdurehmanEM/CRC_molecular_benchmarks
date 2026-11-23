
# 24_discovery_sensitivity.R  (protocol sections 6, 8, 9 - pre-specified sensitivity analyses; DISCOVERY only, read-only)
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(glmnet) })
SEED <- 20261008; MAX_GENES <- 4; ALPHA <- 0.5; K <- 5; REPS <- 10
pb_dir <- "04_processed_data/GSE178341"; prog_dir <- "04_processed_data/programs"
res_dir <- "06_results/sensitivity_discovery"; dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
programs <- fread(file.path(prog_dir, "program_definitions.csv"))[retained == TRUE]
part0 <- fread(file.path(prog_dir, "gene_partition.tsv"))
dm <- fread(file.path(pb_dir, "donor_metadata_tumour.csv"))
frz <- readRDS("04_processed_data/frozen_panel/frozen_panel_v2.rds")
zfit <- function(M) { mu <- colMeans(M); s <- apply(M, 2, sd); s[!is.finite(s) | s == 0] <- 1; list(mu = mu, sd = s) }
zapply <- function(M, z) sweep(sweep(M, 2, z$mu, "-"), 2, z$sd, "/")
fit_sel <- function(Etr, Ptr, Ete, Pte) {
  ze <- zfit(Etr); y_tr <- rowMeans(zapply(Etr, ze)); y_te <- rowMeans(zapply(Ete, ze))
  keep <- apply(Ptr, 2, sd) > 0; Ptr <- Ptr[, keep, drop = FALSE]; Pte <- Pte[, keep, drop = FALSE]
  zp <- zfit(Ptr); Xtr <- zapply(Ptr, zp); Xte <- zapply(Pte, zp)
  cv <- suppressWarnings(cv.glmnet(Xtr, y_tr, alpha = ALPHA, nfolds = min(5, nrow(Xtr)), standardize = FALSE, grouped = FALSE))
  ok <- which(cv$nzero >= 1 & cv$nzero <= MAX_GENES); if (!length(ok)) ok <- which(cv$nzero >= 1)[1]
  lam <- cv$lambda[ok][which.min(cv$cvm[ok])]
  b <- as.matrix(coef(cv, s = lam)); sel <- setdiff(rownames(b)[b[, 1] != 0], "(Intercept)")
  rf <- lm(y_tr ~ ., data = data.frame(y_tr = y_tr, Xtr[, sel, drop = FALSE], check.names = FALSE))
  list(y_te = y_te, pred = as.numeric(predict(rf, newdata = data.frame(Xte[, sel, drop = FALSE], check.names = FALSE))),
       base = rep(mean(y_tr), nrow(Xte)), sel = sel)
}
make_partition <- function(seed) {
  g <- unique(part0[, .(program, row_id, gene_symbol)]); set.seed(seed)
  g[, role := { r <- rep(c("evaluation", "predictor"), length.out = .N); sample(r) }, by = program]
  ev_all <- g[role == "evaluation", unique(row_id)]
  g[!(role == "predictor" & row_id %in% ev_all)] }
load_L <- function(comp, suffix, min_cells, tmm = FALSE, donors_keep = NULL) {
  f <- file.path(pb_dir, paste0("pseudobulk_", comp, suffix, ".rds")); if (!file.exists(f)) return(NULL)
  pb <- readRDS(f); d <- names(pb$n_cells)[pb$n_cells >= min_cells]
  if (!is.null(donors_keep)) d <- intersect(d, donors_keep)
  C <- as.matrix(pb$counts[, d, drop = FALSE])
  if (tmm) { y <- calcNormFactors(DGEList(C)); edgeR::cpm(y, log = TRUE, prior.count = 2) } else edgeR::cpm(C, log = TRUE, prior.count = 2)
}
run_cv <- function(L, ev, pr, chem) {
  E <- t(L[intersect(ev, rownames(L)), , drop = FALSE]); P <- t(L[intersect(pr, rownames(L)), , drop = FALSE]); n <- nrow(E)
  if (n < 10) return(NULL)
  out <- list()
  for (r in seq_len(REPS)) { set.seed(SEED + r); fold <- integer(n)
    for (lv in unique(chem)) { ix <- which(chem == lv); fold[ix] <- sample(rep(seq_len(K), length.out = length(ix))) }
    for (fk in seq_len(K)) { te <- which(fold == fk); tr <- which(fold != fk); if (length(te) < 1 || length(tr) < 5) next
      o <- fit_sel(E[tr, , drop = FALSE], P[tr, , drop = FALSE], E[te, , drop = FALSE], P[te, , drop = FALSE])
      out[[length(out) + 1]] <- data.table(rep = r, y = o$y_te, p = o$pred, b = o$base) } }
  o <- rbindlist(out); m <- o[, .(r = mean(abs(y - p)) / mean(abs(y - b)), rho = suppressWarnings(cor(y, p, method = "spearman"))), by = rep]
  fin <- fit_sel(E, P, E, P)
  list(n = n, ratio = mean(m$r), ratio_min = min(m$r), ratio_max = max(m$r), rho = mean(m$rho, na.rm = TRUE), sel = fin$sel)
}
SC <- rbind(data.table(scenario = paste0("min_cells_", c(10, 30, 50)), min_cells = c(10, 30, 50), seed = SEED, suffix = "", tmm = FALSE, v2 = FALSE),
  data.table(scenario = paste0("partition_seed_", 1:5), min_cells = 20, seed = 1:5, suffix = "", tmm = FALSE, v2 = FALSE),
  data.table(scenario = "TMM", min_cells = 20, seed = SEED, suffix = "", tmm = TRUE, v2 = FALSE),
  data.table(scenario = "10x_v2_only", min_cells = 20, seed = SEED, suffix = "", tmm = FALSE, v2 = TRUE),
  data.table(scenario = "all_sorting_fractions", min_cells = 20, seed = SEED, suffix = "_allfrac", tmm = FALSE, v2 = FALSE),
  data.table(scenario = "all_epithelial_states", min_cells = 20, seed = SEED, suffix = "_all", tmm = FALSE, v2 = FALSE))
chem_col <- grep("chem|SINGLECELL", names(dm), value = TRUE, ignore.case = TRUE)[1]; id_col <- grep("^PID$|donor", names(dm), value = TRUE)[1]
cat("Donor metadata columns used: id =", id_col, "| chemistry =", chem_col, "\n")
v2_donors <- dm[grepl("v2", get(chem_col), ignore.case = TRUE), get(id_col)]
RES <- list()
for (s in seq_len(nrow(SC))) { sc <- SC[s]
  part <- if (sc$seed == SEED) part0 else make_partition(sc$seed)
  for (k in seq_len(nrow(programs))) { pg <- programs$program[k]; comp <- programs$compartment[k]
    if (sc$scenario == "all_epithelial_states" && comp != "Epithelial") next
    L <- load_L(comp, sc$suffix, sc$min_cells, sc$tmm, if (sc$v2) v2_donors else NULL)
    if (is.null(L)) { RES[[paste(s, pg)]] <- data.table(scenario = sc$scenario, program = pg, note = "profile file not found"); next }
    chem <- dm[[chem_col]][match(colnames(L), dm[[id_col]])]; chem[is.na(chem)] <- "NA"
    r <- run_cv(L, part[program == pg & role == "evaluation", row_id], part[program == pg & role == "predictor", row_id], chem)
    if (is.null(r)) { RES[[paste(s, pg)]] <- data.table(scenario = sc$scenario, program = pg, note = "<10 donors"); next }
    fz <- frz$programs[[pg]]$panel$row_id
    RES[[paste(s, pg)]] <- data.table(scenario = sc$scenario, program = pg, n_donors = r$n, MAE_ratio = round(r$ratio, 3),
      ratio_range = sprintf("%.2f-%.2f", r$ratio_min, r$ratio_max), spearman = round(r$rho, 3),
      frozen_genes_reselected = paste0(sum(fz %in% r$sel), "/", length(fz)),
      selected = paste(part[match(r$sel, row_id), gene_symbol], collapse = ";"))
    message(sc$scenario, " | ", pg, " | n = ", r$n, " | MAE ratio ", round(r$ratio, 3)) } }
RES <- rbindlist(RES, fill = TRUE)
base <- fread("06_results/panel/cv_performance_discovery_v2.csv")[method == "panel_refit", .(program, primary = MAE_vs_mean)]
RES <- merge(RES, base, by = "program", all.x = TRUE)
fwrite(RES, file.path(res_dir, "discovery_sensitivity_results.csv"))
cat("\n=== MAE ratio vs mean (primary analysis in 'primary' column) ===\n")
print(dcast(RES[!is.na(MAE_ratio)], program + primary ~ scenario, value.var = "MAE_ratio"))
cat("\n=== Frozen panel genes re-selected (final fit, all donors) ===\n")
print(dcast(RES[!is.na(MAE_ratio)], program ~ scenario, value.var = "frozen_genes_reselected"))
cat("\n=== Donors per scenario ===\n"); print(dcast(RES[!is.na(MAE_ratio)], program ~ scenario, value.var = "n_donors"))
if (any(!is.na(RES$note))) { cat("\nNotes:\n"); print(RES[!is.na(note), .(scenario, program, note)]) }
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_24.txt"))
 
