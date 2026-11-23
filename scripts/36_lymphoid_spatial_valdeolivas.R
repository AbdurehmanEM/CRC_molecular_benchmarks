
# 36_lymphoid_spatial_valdeolivas.R  (protocol v1.9, section 14.11) - criteria A and B for LYM_IL2_STAT5.
# Mirrors 19_spatial_support.R. Composition covariate is T_NK (declared proxy, s14.11.1).
# Criterion B uses T/NK proportion quartiles (s14.11.2). Writes to a separate directory.
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(digest) })
set.seed(20261013)
NB <- 1000; BLK <- 8
frz_file <- "04_processed_data/frozen_panel/frozen_panel_LYM_v1.rds"
if (digest(file = frz_file, algo = "sha256") != "7dc3256d97477a796fcacd1df71021740888f7978a01ea36caf70787da7f513f")
  stop("Lymphoid frozen panel hash mismatch - STOP")
frz <- readRDS(frz_file)
dst <- "04_processed_data/Valdeolivas_Visium/unzipped"; rdir <- "04_processed_data/Valdeolivas_Visium/rctd"
res_dir <- "06_results/spatial_support_lymphoid"; dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
comp2type <- c(Lymphoid = "T_NK")
need <- unique(unlist(lapply(frz$programs, function(P) c(P$panel$gene_symbol, P$evaluation_genes$gene_symbol))))
zs <- function(M) { Z <- scale(M); Z[!is.finite(Z)] <- 0; Z }
fitb <- function(y, p, Wm) { X <- cbind(1, p, Wm); b <- qr.coef(qr(X), y); unname(b[2]) }
EXCLUDE_PATIENTS <- c("A798015")
spots <- list(); A_sec <- list(); boots <- list(); cover <- list()
for (rf in list.files(rdir, "_rctd.rds$", full.names = TRUE)) {
  x <- readRDS(rf); stopifnot(isTRUE(x$labels_verified)); if (x$patient %in% EXCLUDE_PATIENTS) next
  s <- x$section; W <- x$weights; bc <- rownames(W)
  cnt <- Seurat::Read10X_h5(file.path(dst, s, "filtered_feature_bc_matrix.h5"))[, bc]
  lib <- Matrix::colSums(cnt); g <- intersect(need, rownames(cnt))
  L <- edgeR::cpm(as.matrix(cnt[g, , drop = FALSE]), lib.size = lib, log = TRUE, prior.count = 2)
  pos <- x$pos[match(bc, barcode)]; blk <- paste(pos$array_row %/% BLK, pos$array_col %/% BLK)
  Wm <- W[, setdiff(colnames(W), "Mast"), drop = FALSE]
  for (pg in names(frz$programs)) {
    P <- frz$programs[[pg]]; ev <- intersect(P$evaluation_genes$gene_symbol, rownames(L))
    evc <- length(ev) / nrow(P$evaluation_genes); pan_ok <- all(P$panel$gene_symbol %in% rownames(L))
    cover[[paste(s, pg)]] <- data.table(section = s, program = pg,
      panel_present = sum(P$panel$gene_symbol %in% rownames(L)), panel_genes = nrow(P$panel),
      eval_coverage = round(evc, 3))
    if (!pan_ok || evc < 0.80) next
    y <- rowMeans(zs(t(L[ev, , drop = FALSE])))
    p <- as.numeric(P$intercept + zs(t(L[P$panel$gene_symbol, , drop = FALSE])) %*% P$panel$coefficient)
    ys <- as.numeric(scale(y)); ps <- if (sd(p) > 0) as.numeric(scale(p)) else rep(0, length(p))
    comp_prop <- W[, comp2type[[P$compartment]]]
    spots[[paste(s, pg)]] <- data.table(section = s, patient = x$patient, program = pg, barcode = bc,
      array_row = pos$array_row, array_col = pos$array_col, pxl_row = pos$pxl_row, pxl_col = pos$pxl_col,
      target = y, panel = p, comp_prop = comp_prop, tnk_prop = W[, "T_NK"])
    est <- fitb(ys, ps, Wm); ub <- unique(blk); bi <- split(seq_along(blk), blk)
    bs <- replicate(NB, { i <- unlist(bi[sample(ub, length(ub), TRUE)], use.names = FALSE)
      fitb(ys[i], ps[i], Wm[i, , drop = FALSE]) })
    q <- comp_prop >= quantile(comp_prop, 0.75)
    est_q <- if (sum(q) >= 30) fitb(ys[q], ps[q], Wm[q, , drop = FALSE]) else NA_real_
    A_sec[[paste(s, pg)]] <- data.table(section = s, patient = x$patient, program = pg, n_spots = length(ys),
      effect = est, lo = quantile(bs, 0.025, na.rm = TRUE), hi = quantile(bs, 0.975, na.rm = TRUE),
      effect_topQ_comp = est_q, raw_spearman = cor(y, p, method = "spearman"))
    boots[[paste(s, pg)]] <- bs }
  message("done ", s) }
spots <- rbindlist(spots); A_sec <- rbindlist(A_sec); cover <- rbindlist(cover)
saveRDS(spots, "04_processed_data/Valdeolivas_Visium/spot_program_scores_lymphoid.rds")
A_pat <- A_sec[, { k <- paste(section, program); B <- do.call(cbind, boots[k]); bm <- rowMeans(B, na.rm = TRUE)
  .(n_sections = .N, n_spots = sum(n_spots), effect = mean(effect), se = sd(bm, na.rm = TRUE),
    effect_topQ_comp = mean(effect_topQ_comp, na.rm = TRUE)) }, by = .(patient, program)]
dl <- function(est, se) { w <- 1 / se^2; fe <- sum(w * est) / sum(w); k <- length(est); Q <- sum(w * (est - fe)^2)
  t2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))); wr <- 1 / (se^2 + t2)
  m <- sum(wr * est) / sum(wr); s <- sqrt(1 / sum(wr))
  c(est = m, lo = m - 1.96 * s, hi = m + 1.96 * s, I2 = if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0) }
Bp <- spots[, { rule <- "tnk_quartile"
  hi <- ave(tnk_prop, section, FUN = function(v) v >= quantile(v, 0.75)) == 1
  lo <- ave(tnk_prop, section, FUN = function(v) v <= quantile(v, 0.25)) == 1
  .(rule = rule, n_hi = sum(hi), n_lo = sum(lo),
    diff = if (sum(hi) >= 10 && sum(lo) >= 10) mean(target[hi]) - mean(target[lo]) else NA_real_)
}, by = .(patient, program)]
S <- rbindlist(lapply(names(frz$programs), function(pg) {
  a <- A_pat[program == pg & is.finite(se) & se > 0]; b <- Bp[program == pg & is.finite(diff)]
  np <- nrow(a); need_k <- ceiling(2 / 3 * max(np, 1))
  if (np < 3) return(data.table(program = pg, n_patients = np, spatial_label = "Not evaluable"))
  m <- dl(a$effect, a$se); A_met <- m[["lo"]] > 0 && sum(a$effect > 0) >= need_k
  B_met <- nrow(b) >= 3 && sum(b$diff > 0) >= ceiling(2 / 3 * nrow(b))
  data.table(program = pg, n_patients = np, pooled_effect = m[["est"]], lo = m[["lo"]], hi = m[["hi"]], I2 = m[["I2"]],
    patients_positive_A = sum(a$effect > 0), A_met = A_met, B_patients = nrow(b),
    patients_positive_B = sum(b$diff > 0), B_met = B_met,
    pooled_topQ_effect = mean(a$effect_topQ_comp, na.rm = TRUE),
    spatial_label = if (A_met && B_met) "Spatial support" else if (A_met || B_met) "Partial spatial support" else "No spatial support")
}), fill = TRUE)
rnd <- function(X) { nn <- names(X)[vapply(X, is.numeric, logical(1))]; X[, (nn) := lapply(.SD, round, 3), .SDcols = nn][] }
fwrite(rnd(cover), file.path(res_dir, "spatial_gene_coverage.csv"))
fwrite(rnd(A_sec), file.path(res_dir, "criterionA_by_section.csv"))
fwrite(rnd(A_pat), file.path(res_dir, "criterionA_by_patient.csv"))
fwrite(rnd(Bp), file.path(res_dir, "criterionB_by_patient.csv"))
fwrite(rnd(S), file.path(res_dir, "spatial_support_summary.csv"))
cat("\n=== Coverage failures (if any) ===\n"); print(cover[panel_present < panel_genes | eval_coverage < 0.8])
cat("\n=== Criterion A by patient ===\n"); print(dcast(rnd(A_pat), program ~ patient, value.var = "effect"))
cat("\n=== Criterion B by patient (T/NK quartile contrast) ===\n"); print(dcast(rnd(Bp), program ~ patient, value.var = "diff"))
cat("\n=== LYMPHOID SPATIAL SUMMARY (Valdeolivas, descriptive) ===\n"); print(rnd(S))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_36.txt"))
 
