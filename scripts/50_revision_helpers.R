###############################################################################
# 50_revision_helpers.R - shared functions for revision scripts 51-64.
# Re-implements the published single-cell pipeline exactly as described in
# protocol v0.6-v0.9 and Supplementary Tables S5/S6, so that new analyses use
# the same targets, scaling and frozen coefficients. Script 52 checks this.
###############################################################################
suppressPackageStartupMessages({
  library(Matrix); library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

init_revision <- function() {
  for (d in c(OUT_DIR, TAB_DIR, FIG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  set.seed(SEED)
  invisible(TRUE)
}

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

## ---------------------------------------------------------------------------
## Supplementary workbook readers (header row located automatically)
## ---------------------------------------------------------------------------
read_supp <- function(sheet, first_col = c("Program code", "Dataset", "Patient", "Table")) {
  raw <- readxl::read_excel(SUPP_XLSX, sheet = sheet, col_names = FALSE, .name_repair = "minimal")
  hdr <- which(as.character(raw[[1]]) %in% first_col)[1]
  if (is.na(hdr)) stop("Header row not found in sheet ", sheet)
  dt <- readxl::read_excel(SUPP_XLSX, sheet = sheet, skip = hdr - 1, .name_repair = "unique")
  dt <- as.data.table(dt)
  dt <- dt[!is.na(dt[[1]]) & dt[[1]] != ""]
  dt
}

## Rename columns by ASCII regex (robust to non-UTF-8 locales: rho, superscript 2)
setnames_re <- function(dt, patterns, new) {
  stopifnot(length(patterns) == length(new))
  for (i in seq_along(patterns)) {
    j <- grep(patterns[i], names(dt), perl = TRUE)
    if (length(j) != 1) stop("Column matching '", patterns[i], "' found ", length(j), " times in: ",
                             paste(names(dt), collapse = " | "))
    setnames(dt, j, new[i])
  }
  invisible(dt)
}

## Frozen panels (S6) and gene partition (S5)
load_panels <- function() {
  s6 <- read_supp("S6_frozen_panels")
  setnames(s6, c("Program code", "Program", "Compartment", "Gene", "Matrix row ID",
                 "Coefficient", "Scaling mean", "Scaling SD", "Intercept", "Evaluation genes"),
           c("program", "program_name", "compartment", "symbol", "row_id",
             "coef", "mu", "sd", "intercept", "n_eval"))
  num <- c("coef", "mu", "sd", "intercept", "n_eval")
  s6[, (num) := lapply(.SD, as.numeric), .SDcols = num]
  split(s6, by = "program", keep.by = TRUE)
}
load_partition <- function() {
  s5 <- read_supp("S5_gene_partition")
  setnames(s5, c("Program code", "Program", "Compartment", "Ensembl ID", "Gene symbol",
                 "Matrix row ID", "Expressed", "Role"),
           c("program", "program_name", "compartment", "ensembl", "symbol", "row_id",
             "expressed", "role"))
  s5[]
}
program_compartment <- function(panels) {
  sapply(panels, function(p) {
    cmp <- unique(p$compartment)
    if (cmp %in% c("Lymphoid", "CD8", "CD8 T cells")) "Lymphoid" else cmp
  })
}

## ---------------------------------------------------------------------------
## Pseudobulk loading
## ---------------------------------------------------------------------------
load_pb <- function(cohort, compartment, min_cells = MIN_CELLS, quiet = FALSE) {
  f <- PB_FILES[[cohort]][[compartment]]
  if (is.null(f) || is.na(f)) return(NULL)
  if (!file.exists(f)) stop("Pseudobulk file not found: ", f, "\n  Edit PB_FILES in 50_revision_config.R")
  obj <- readRDS(f)
  n_cells <- NULL
  if (is.list(obj) && !is.data.frame(obj)) {
    cnt <- obj$counts %||% obj$pb %||% obj$mat %||% obj[[which(sapply(obj, function(x) is.matrix(x) || inherits(x, "Matrix")))[1]]]
    n_cells <- obj$n_cells %||% obj$ncells %||% obj$cells
    if (is.data.frame(n_cells)) n_cells <- setNames(n_cells[[2]], n_cells[[1]])
  } else cnt <- obj
  if (is.data.frame(cnt)) cnt <- as.matrix(cnt)
  cnt <- as.matrix(cnt)
  storage.mode(cnt) <- "double"
  if (any(abs(cnt - round(cnt)) > 1e-8, na.rm = TRUE)) stop(f, ": counts are not integers (normalised matrix?)")
  if (!is.null(n_cells)) {
    keep <- colnames(cnt)[n_cells[colnames(cnt)] >= min_cells]
    if (!quiet) msg(cohort, compartment, ": kept", length(keep), "of", ncol(cnt), "donors with >=", min_cells, "cells")
    cnt <- cnt[, keep, drop = FALSE]
  } else if (!quiet) msg(cohort, compartment, ": no cell counts in file; assuming all", ncol(cnt), "donors are eligible")
  cnt
}

## log2-CPM, library size only, prior count 2 (protocol v0.7)
logcpm <- function(counts, prior = PRIOR_COUNT) {
  if (requireNamespace("edgeR", quietly = TRUE)) {
    return(edgeR::cpm(counts, log = TRUE, prior.count = prior))
  }
  lib <- colSums(counts)
  adj <- prior * lib / mean(lib)
  lib2 <- lib + 2 * adj
  t(log2(t(counts) + adj) - log2(lib2) + log2(1e6))
}

## Map keys: discovery uses row IDs; validation uses symbols
key_col <- function(cohort) if (cohort == DISCOVERY) "row_id" else "symbol"

## ---------------------------------------------------------------------------
## Targets and panel predictions on the FROZEN discovery scale
## ---------------------------------------------------------------------------
# Discovery statistics of the evaluation genes (mean/SD over eligible discovery donors)
discovery_eval_stats <- function(lc_disc, eval_rows) {
  g <- intersect(eval_rows, rownames(lc_disc))
  x <- lc_disc[g, , drop = FALSE]
  data.table(row_id = g, mu = rowMeans(x), sd = apply(x, 1, sd))
}

# Target = mean z-score of evaluation genes present (z with discovery mean/SD)
target_score <- function(lc, eval_tab, key = "row_id") {
  et <- eval_tab[!is.na(get(key)) & sd > 0]
  present <- et[get(key) %in% rownames(lc)]
  cover <- nrow(present) / nrow(eval_tab)
  if (nrow(present) == 0) return(list(y = rep(NA_real_, ncol(lc)), coverage = 0))
  z <- (lc[present[[key]], , drop = FALSE] - present$mu) / present$sd
  list(y = colMeans(z), coverage = cover)
}

# Frozen prediction = intercept + sum(coef * (x - mu) / sd)
panel_predict <- function(lc, panel, key = "row_id") {
  k <- panel[[key]]
  if (!all(k %in% rownames(lc))) return(NULL)   # missing panel gene -> not evaluable
  z <- (lc[k, , drop = FALSE] - panel$mu) / panel$sd
  as.numeric(panel$intercept[1] + colSums(z * panel$coef))
}

mae_ratio <- function(y, yhat, base = 0) mean(abs(yhat - y)) / mean(abs(y - base))

boot_ratio <- function(y, yhat, B = N_BOOT, base = 0) {
  n <- length(y)
  r <- replicate(B, { i <- sample.int(n, n, replace = TRUE); mae_ratio(y[i], yhat[i], base) })
  r <- r[is.finite(r)]
  list(ci = unname(quantile(r, c(0.025, 0.975))), draws = r,
       p_one_sided = (1 + sum(r >= 1)) / (1 + length(r)))   # H0: ratio >= 1
}

## Assemble everything needed to evaluate one program in one cohort
eval_table <- function(prog, partition, disc_stats) {
  merge(disc_stats[[prog]], partition[program == prog & role == "evaluation", .(row_id, symbol)],
        by = "row_id", all.x = TRUE)
}
eval_program_cohort <- function(prog, cohort, panels, partition, disc_stats, lc_cache, panel = NULL) {
  panel <- panel %||% panels[[prog]]
  cmp <- program_compartment(panels)[[prog]]
  lc <- lc_cache[[cohort]][[cmp]]
  if (is.null(lc)) return(NULL)
  key <- key_col(cohort)
  et <- eval_table(prog, partition, disc_stats)
  tg <- target_score(lc, et, key)
  yhat <- panel_predict(lc, panel, key)
  list(program = prog, cohort = cohort, y = setNames(tg$y, colnames(lc)), yhat = yhat,
       coverage = tg$coverage, n = ncol(lc),
       evaluable = !is.null(yhat) && tg$coverage >= MIN_EVAL_COVER && ncol(lc) >= MIN_DONORS)
}

## Build log-CPM cache for every cohort x compartment needed
build_lc_cache <- function(cohorts, compartments) {
  cache <- list()
  for (co in cohorts) {
    cache[[co]] <- list()
    for (cmp in compartments) {
      cnt <- tryCatch(load_pb(co, cmp, quiet = TRUE), error = function(e) { msg("WARN", conditionMessage(e)); NULL })
      if (!is.null(cnt) && ncol(cnt) > 0) cache[[co]][[cmp]] <- logcpm(cnt)
    }
  }
  cache
}

all_disc_stats <- function(panels, partition, lc_cache) {
  cmp <- program_compartment(panels)
  setNames(lapply(names(panels), function(p) {
    discovery_eval_stats(lc_cache[[DISCOVERY]][[cmp[[p]]]],
                         partition[program == p & role == "evaluation", row_id])
  }), names(panels))
}

## ---------------------------------------------------------------------------
## Relaxed elastic net (protocol v0.8): select 1-4 genes, OLS refit
## ---------------------------------------------------------------------------
select_relaxed_en <- function(X, y, max_genes = 4, alpha = 0.5, nfolds = 5) {
  # X: donors x genes, already z-scored with training statistics
  if (ncol(X) < 2) return(colnames(X))
  cv <- glmnet::cv.glmnet(X, y, alpha = alpha, standardize = FALSE, nfolds = min(nfolds, nrow(X)))
  df <- cv$glmnet.fit$df
  ok <- which(df >= 1 & df <= max_genes)
  if (length(ok) == 0) ok <- which(df >= 1)[1]
  best <- ok[which.min(cv$cvm[match(cv$glmnet.fit$lambda[ok], cv$lambda)])]
  b <- as.matrix(glmnet::coef.glmnet(cv$glmnet.fit, s = cv$glmnet.fit$lambda[best]))[-1, 1]
  sel <- names(b)[b != 0]
  if (length(sel) > max_genes) sel <- names(sort(abs(b[sel]), decreasing = TRUE))[1:max_genes]
  sel
}
ols_fit <- function(X, y) {
  fit <- lm.fit(cbind(1, X), y)
  b <- fit$coefficients; b[is.na(b)] <- 0
  list(intercept = unname(b[1]), coef = setNames(b[-1], colnames(X)))
}
# z-score helpers
zstats <- function(M) list(mu = rowMeans(M), sd = pmax(apply(M, 1, sd), 1e-8))
zapply <- function(M, s) (M - s$mu) / s$sd

## Fit a FIXED gene set on discovery donors and return a frozen-style panel
fit_fixed_panel <- function(lc_disc, genes, y_disc, sym_map) {
  s <- zstats(lc_disc[genes, , drop = FALSE])
  X <- t(zapply(lc_disc[genes, , drop = FALSE], s))
  f <- ols_fit(X, y_disc)
  data.table(row_id = genes, symbol = unname(sym_map[genes]), coef = unname(f$coef[genes]),
             mu = unname(s$mu), sd = unname(s$sd), intercept = f$intercept)
}

## Discovery feature map row_id -> symbol (S5 genes always; all genes if file given)
load_symbol_map <- function(partition, lc_disc = NULL) {
  m <- setNames(partition$symbol, partition$row_id)
  f <- DISCOVERY_FEATURES
  if (!is.na(f) && file.exists(f)) {
    ft <- fread(f); setnames(ft, 1:2, c("row_id", "symbol"))
    m <- c(m, setNames(ft$symbol, ft$row_id)); m <- m[!duplicated(names(m))]
  } else if (!is.null(lc_disc) && !any(grepl("^ENSG", rownames(lc_disc)))) {
    m <- c(m, setNames(rownames(lc_disc), rownames(lc_disc))); m <- m[!duplicated(names(m))]
  }
  m
}

## Donor-grouped CV (one profile per donor), optional selection inside folds.
## Returns out-of-fold predictions averaged over repeats, and CV MAE ratio.
cv_panel <- function(lc_disc, eval_rows, genes_pool, select = TRUE, fixed_genes = NULL,
                     folds = CV_FOLDS, repeats = CV_REPEATS, strata = NULL, max_genes = 4) {
  donors <- colnames(lc_disc); n <- length(donors)
  oof <- matrix(NA_real_, n, repeats, dimnames = list(donors, NULL))
  base <- matrix(NA_real_, n, repeats, dimnames = list(donors, NULL))
  ytrue <- matrix(NA_real_, n, repeats, dimnames = list(donors, NULL))
  for (r in seq_len(repeats)) {
    fid <- integer(n)
    if (is.null(strata)) fid <- sample(rep_len(seq_len(folds), n)) else {
      for (s in unique(strata)) { i <- which(strata == s); fid[i] <- sample(rep_len(seq_len(folds), length(i))) }
    }
    for (k in seq_len(folds)) {
      tr <- donors[fid != k]; te <- donors[fid == k]
      es <- zstats(lc_disc[eval_rows, tr, drop = FALSE]); ok <- es$sd > 1e-8
      ytr <- colMeans(zapply(lc_disc[eval_rows[ok], tr, drop = FALSE], lapply(es, `[`, ok)))
      yte <- colMeans(zapply(lc_disc[eval_rows[ok], te, drop = FALSE], lapply(es, `[`, ok)))
      if (select) {
        ps <- zstats(lc_disc[genes_pool, tr, drop = FALSE]); okp <- ps$sd > 1e-8
        gp <- genes_pool[okp]
        Xtr <- t(zapply(lc_disc[gp, tr, drop = FALSE], lapply(ps, `[`, okp)))
        sel <- select_relaxed_en(Xtr, ytr, max_genes = max_genes)
      } else sel <- fixed_genes
      gs <- zstats(lc_disc[sel, tr, drop = FALSE])
      f <- ols_fit(t(zapply(lc_disc[sel, tr, drop = FALSE], gs)), ytr)
      pred <- f$intercept + colSums(zapply(lc_disc[sel, te, drop = FALSE], gs) * f$coef)
      oof[te, r] <- pred; base[te, r] <- mean(ytr); ytrue[te, r] <- yte
    }
  }
  ratio_rep <- sapply(seq_len(repeats), function(r) mean(abs(oof[, r] - ytrue[, r])) / mean(abs(ytrue[, r] - base[, r])))
  list(ratio = mean(ratio_rep), ratio_rep = ratio_rep, oof = oof, base = base, ytrue = ytrue)
}

## ---------------------------------------------------------------------------
## Expression-matched random gene draws (script 53, 59)
## ---------------------------------------------------------------------------
expr_bins <- function(lc, genes, n_mean_bins = 10, n_sd_bins = 3) {
  m <- rowMeans(lc[genes, , drop = FALSE]); s <- apply(lc[genes, , drop = FALSE], 1, sd)
  mb <- cut(m, unique(quantile(m, seq(0, 1, length.out = n_mean_bins + 1))), include.lowest = TRUE, labels = FALSE)
  sb <- integer(length(genes))
  for (b in unique(mb)) { i <- which(mb == b)
    sb[i] <- cut(s[i], unique(quantile(s[i], seq(0, 1, length.out = n_sd_bins + 1))), include.lowest = TRUE, labels = FALSE) }
  setNames(paste(mb, sb, sep = "_"), genes)
}
draw_matched <- function(target_genes, bins, universe) {
  ub <- bins[universe]
  vapply(target_genes, function(g) {
    cand <- universe[ub == bins[[g]]]
    if (length(cand) == 0) cand <- universe
    cand[sample.int(length(cand), 1)]
  }, character(1))
}

## ---------------------------------------------------------------------------
## Gate status shared between scripts
## ---------------------------------------------------------------------------
gate_file <- function() file.path(OUT_DIR, "52_gate_status.txt")
require_gate <- function() {
  if (!GATE_STRICT) { msg("GATE_STRICT = FALSE: running without the reproduction gate (testing only)"); return(invisible()) }
  if (!file.exists(gate_file()) || !identical(readLines(gate_file())[1], "PASS"))
    stop("Reproduction gate has not passed. Run 52_reproduce_published.R first and resolve differences.")
}

## DerSimonian-Laird pooling (identical to scripts 16/23)
dl_pool <- function(yi, sei) {
  w <- 1 / sei^2; k <- length(yi)
  fe <- sum(w * yi) / sum(w); Q <- sum(w * (yi - fe)^2)
  tau2 <- if (k > 1) max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))) else 0
  ws <- 1 / (sei^2 + tau2); est <- sum(ws * yi) / sum(ws); se <- sqrt(1 / sum(ws))
  I2 <- if (k > 1 && Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0
  list(est = est, se = se, lo = est - 1.959964 * se, hi = est + 1.959964 * se, tau2 = tau2, I2 = I2, k = k)
}
se_from_ci <- function(lo, hi) (log(hi) - log(lo)) / (2 * 1.959964)

write_out <- function(dt, name) {
  f <- file.path(OUT_DIR, paste0(name, ".csv")); fwrite(dt, f); msg("wrote", f); invisible(f)
}
