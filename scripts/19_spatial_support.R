
# 19_spatial_support.R  (protocol v1.2, section 10d) - criteria A (composition-adjusted concordance) and B (anatomy)
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(digest) })
set.seed(20261013)
NB <- 1000; BLK <- 8
frz_file <- "04_processed_data/frozen_panel/frozen_panel_v2.rds"
if (digest(file = frz_file, algo = "sha256") != "a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb") stop("Frozen panel hash mismatch - STOP")
frz <- readRDS(frz_file)
dst <- "04_processed_data/Valdeolivas_Visium/unzipped"; rdir <- "04_processed_data/Valdeolivas_Visium/rctd"
res_dir <- "06_results/spatial_support"; dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
comp2type <- c(Epithelial = "Malignant_epithelial", Stromal = "Fibroblast", Myeloid = "Myeloid", Endothelial = "Endothelial")
need <- unique(unlist(lapply(frz$programs, function(P) c(P$panel$gene_symbol, P$evaluation_genes$gene_symbol))))
pa_files <- list.files(dst, pattern = "^Pathologist_Annotations_.*\\.csv$", recursive = TRUE, full.names = TRUE)
normlab <- function(x) { x <- tolower(trimws(x)); x <- gsub("_", " ", x); x <- gsub("\\s+", " ", x)
  x <- gsub("epitehlium", "epithelium", x); x <- gsub("aggr[a-z]*", "aggregate", x); x }
zs <- function(M) { Z <- scale(M); Z[!is.finite(Z)] <- 0; Z }
fitb <- function(y, p, Wm) { X <- cbind(1, p, Wm); b <- qr.coef(qr(X), y); unname(b[2]) }
EXCLUDE_PATIENTS <- c("A798015")   # non-neoplastic tissue (Valdeolivas clinical table) - fails eligibility, protocol 10d
spots <- list(); A_sec <- list(); boots <- list(); cover <- list()
for (rf in list.files(rdir, "_rctd.rds$", full.names = TRUE)) {
  x <- readRDS(rf); stopifnot(isTRUE(x$labels_verified)); if (x$patient %in% EXCLUDE_PATIENTS) next; s <- x$section; W <- x$weights; bc <- rownames(W)
  cnt <- Seurat::Read10X_h5(file.path(dst, s, "filtered_feature_bc_matrix.h5"))[, bc]
  lib <- Matrix::colSums(cnt); g <- intersect(need, rownames(cnt))
  L <- cpm(as.matrix(cnt[g, , drop = FALSE]), lib.size = lib, log = TRUE, prior.count = 2)
  pos <- x$pos[match(bc, barcode)]; blk <- paste(pos$array_row %/% BLK, pos$array_col %/% BLK)
  pf <- grep(paste0("_", s, "\\.csv$"), pa_files, value = TRUE)
  lab <- if (length(pf)) { a <- fread(pf[1]); normlab(a[[2]][match(bc, a[[1]])]) } else rep(NA_character_, length(bc))
  Wm <- W[, setdiff(colnames(W), "Mast"), drop = FALSE]
  for (pg in names(frz$programs)) {
    P <- frz$programs[[pg]]; ev <- intersect(P$evaluation_genes$gene_symbol, rownames(L))
    evc <- length(ev) / nrow(P$evaluation_genes); pan_ok <- all(P$panel$gene_symbol %in% rownames(L))
    cover[[paste(s, pg)]] <- data.table(section = s, program = pg, panel_present = sum(P$panel$gene_symbol %in% rownames(L)),
      panel_genes = nrow(P$panel), eval_coverage = round(evc, 3))
    if (!pan_ok || evc < 0.80) next
    y <- rowMeans(zs(t(L[ev, , drop = FALSE])))
    p <- as.numeric(P$intercept + zs(t(L[P$panel$gene_symbol, , drop = FALSE])) %*% P$panel$coefficient)
    ys <- as.numeric(scale(y)); ps <- if (sd(p) > 0) as.numeric(scale(p)) else rep(0, length(p))
    comp_prop <- W[, comp2type[[P$compartment]]]
    spots[[paste(s, pg)]] <- data.table(section = s, patient = x$patient, program = pg, barcode = bc,
      array_row = pos$array_row, array_col = pos$array_col, pxl_row = pos$pxl_row, pxl_col = pos$pxl_col,
      label = lab, target = y, panel = p, comp_prop = comp_prop, myeloid_prop = W[, "Myeloid"], endo_prop = W[, "Endothelial"])
    est <- fitb(ys, ps, Wm); ub <- unique(blk); bi <- split(seq_along(blk), blk)
    bs <- replicate(NB, { i <- unlist(bi[sample(ub, length(ub), TRUE)], use.names = FALSE); fitb(ys[i], ps[i], Wm[i, , drop = FALSE]) })
    q <- comp_prop >= quantile(comp_prop, 0.75)
    est_q <- if (sum(q) >= 30) fitb(ys[q], ps[q], Wm[q, , drop = FALSE]) else NA_real_
    A_sec[[paste(s, pg)]] <- data.table(section = s, patient = x$patient, program = pg, n_spots = length(ys),
      effect = est, lo = quantile(bs, 0.025, na.rm = TRUE), hi = quantile(bs, 0.975, na.rm = TRUE), effect_topQ_comp = est_q,
      raw_spearman = cor(y, p, method = "spearman"))
    boots[[paste(s, pg)]] <- bs
  }
  message("done ", s)
}
spots <- rbindlist(spots); A_sec <- rbindlist(A_sec); cover <- rbindlist(cover)
saveRDS(spots, file.path("04_processed_data/Valdeolivas_Visium", "spot_program_scores.rds"))
A_pat <- A_sec[, { k <- paste(section, program); B <- do.call(cbind, boots[k])
  bm <- rowMeans(B, na.rm = TRUE)
  .(n_sections = .N, n_spots = sum(n_spots), effect = mean(effect), se = sd(bm, na.rm = TRUE),
    effect_topQ_comp = mean(effect_topQ_comp, na.rm = TRUE)) }, by = .(patient, program)]
dl <- function(est, se) { w <- 1 / se^2; fe <- sum(w * est) / sum(w); k <- length(est); Q <- sum(w * (est - fe)^2)
  t2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))); wr <- 1 / (se^2 + t2); m <- sum(wr * est) / sum(wr); s <- sqrt(1 / sum(wr))
  c(est = m, lo = m - 1.96 * s, hi = m + 1.96 * s, I2 = if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0) }
Bp <- spots[, {
  comp <- frz$programs[[program[1]]]$compartment; hi <- lo <- rep(FALSE, .N); rule <- "label"
  if (comp == "Epithelial") { ok <- !is.na(label) & label != "" & label != "exclude"; hi <- ok & grepl("^tumor", label); lo <- ok & !grepl("^tumor", label) }
  if (comp == "Stromal") { hi <- !is.na(label) & grepl("^stroma", label); lo <- !is.na(label) & label == "tumor" }
  if (comp == "Myeloid") { agg <- grepl("aggregate", label)
    hi <- !is.na(label) & !agg & grepl("ic high|ic med to high", label); lo <- !is.na(label) & !agg & grepl("ic low", label)
    if (sum(hi) < 10 || sum(lo) < 10) { rule <- "myeloid_quartile"
      qh <- ave(myeloid_prop, section, FUN = function(v) v >= quantile(v, 0.75)) == 1
      ql <- ave(myeloid_prop, section, FUN = function(v) v <= quantile(v, 0.25)) == 1; hi <- qh; lo <- ql } }
  if (comp == "Endothelial") { rule <- "endothelial_quartile"
    hi <- ave(endo_prop, section, FUN = function(v) v >= quantile(v, 0.75)) == 1
    lo <- ave(endo_prop, section, FUN = function(v) v <= quantile(v, 0.25)) == 1 }
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
    patients_positive_A = sum(a$effect > 0), A_met = A_met, B_patients = nrow(b), patients_positive_B = sum(b$diff > 0), B_met = B_met,
    pooled_topQ_effect = mean(a$effect_topQ_comp, na.rm = TRUE),
    spatial_label = if (A_met && B_met) "Spatial support" else if (A_met || B_met) "Partial spatial support" else "No spatial support")
}), fill = TRUE)
rnd <- function(X) { nn <- names(X)[vapply(X, is.numeric, logical(1))]; X[, (nn) := lapply(.SD, round, 3), .SDcols = nn][] }
fwrite(rnd(cover), file.path(res_dir, "spatial_gene_coverage.csv")); fwrite(rnd(A_sec), file.path(res_dir, "criterionA_by_section.csv"))
fwrite(rnd(A_pat), file.path(res_dir, "criterionA_by_patient.csv")); fwrite(rnd(Bp), file.path(res_dir, "criterionB_by_patient.csv"))
fwrite(rnd(S), file.path(res_dir, "spatial_support_summary.csv"))
cat("\n=== Coverage failures (if any) ===\n"); print(cover[panel_present < panel_genes | eval_coverage < 0.8])
cat("\n=== Criterion A by patient (standardised panel coefficient, composition-adjusted) ===\n")
print(dcast(rnd(A_pat), program ~ patient, value.var = "effect"))
cat("\n=== Criterion B by patient (mean target difference, expected-high minus comparison) ===\n")
print(dcast(rnd(Bp), program ~ patient, value.var = "diff")); print(unique(Bp[, .(program, rule)]))
cat("\n=== SPATIAL SUPPORT SUMMARY ===\n"); print(rnd(S))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_19.txt"))
 
