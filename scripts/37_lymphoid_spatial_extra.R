
# 37_lymphoid_spatial_extra.R  (protocol v1.9, section 14.11) - LYM_IL2_STAT5 in GSE280318 + GSE226997.
# Scoring only: reuses the RCTD weights already saved by script 22. No deconvolution is re-run.
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(digest) })
set.seed(20261013); NB <- 1000; BLK <- 8
frz_file <- "04_processed_data/frozen_panel/frozen_panel_LYM_v1.rds"
if (digest(file = frz_file, algo = "sha256") != "7dc3256d97477a796fcacd1df71021740888f7978a01ea36caf70787da7f513f")
  stop("Lymphoid frozen panel hash mismatch - STOP")
frz <- readRDS(frz_file)
odir <- "04_processed_data/spatial_extra"
res_dir <- "06_results/spatial_support_lymphoid"; dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
rf <- list.files(odir, "_rctd.rds$", full.names = TRUE)
if (!length(rf)) stop("no saved RCTD objects in ", odir, " - run script 22 first")
cat("sections with saved RCTD weights:", length(rf), "\n")
comp2type <- c(Lymphoid = "T_NK")
need <- unique(unlist(lapply(frz$programs, function(P) c(P$panel$gene_symbol, P$evaluation_genes$gene_symbol))))
zs <- function(M) { Z <- scale(M); Z[!is.finite(Z)] <- 0; Z }
fitb <- function(y, p, Wm) { X <- cbind(1, p, Wm); unname(qr.coef(qr(X), y)[2]) }
A <- list(); Bq <- list(); cover <- list()
for (p in rf) {
  x <- readRDS(p); W <- x$weights; bc <- rownames(W); cnt <- x$cnt[, bc]
  L <- edgeR::cpm(as.matrix(cnt[intersect(need, rownames(cnt)), , drop = FALSE]),
                  lib.size = Matrix::colSums(cnt), log = TRUE, prior.count = 2)
  pos <- x$pos[match(bc, barcode)]; blk <- paste(pos$array_row %/% BLK, pos$array_col %/% BLK)
  Wm <- W[, setdiff(colnames(W), "Mast"), drop = FALSE]
  for (pg in names(frz$programs)) {
    P <- frz$programs[[pg]]; ev <- intersect(P$evaluation_genes$gene_symbol, rownames(L))
    evc <- length(ev) / nrow(P$evaluation_genes); pan_ok <- all(P$panel$gene_symbol %in% rownames(L))
    cover[[paste(x$section, pg)]] <- data.table(section = x$section, cohort = x$cohort, program = pg,
      panel_present = sum(P$panel$gene_symbol %in% rownames(L)), panel_genes = nrow(P$panel),
      missing_panel = paste(setdiff(P$panel$gene_symbol, rownames(L)), collapse = ";"), eval_coverage = round(evc, 3))
    if (!pan_ok || evc < 0.80) next
    y <- rowMeans(zs(t(L[ev, , drop = FALSE])))
    pz <- as.numeric(P$intercept + zs(t(L[P$panel$gene_symbol, , drop = FALSE])) %*% P$panel$coefficient)
    ys <- as.numeric(scale(y)); ps <- if (sd(pz) > 0) as.numeric(scale(pz)) else rep(0, length(pz))
    est <- fitb(ys, ps, Wm); ub <- unique(blk); bi <- split(seq_along(blk), blk)
    bs <- replicate(NB, { i <- unlist(bi[sample(ub, length(ub), TRUE)], use.names = FALSE)
      fitb(ys[i], ps[i], Wm[i, , drop = FALSE]) })
    cp <- W[, comp2type[[P$compartment]]]; hi <- cp >= quantile(cp, 0.75); lo <- cp <= quantile(cp, 0.25)
    A[[paste(x$section, pg)]] <- data.table(cohort = x$cohort, section = x$section, patient = x$patient,
      platform = x$platform, role = x$role, program = pg, n_spots = length(ys), effect = est,
      se = sd(bs, na.rm = TRUE), lo = quantile(bs, 0.025, na.rm = TRUE), hi = quantile(bs, 0.975, na.rm = TRUE))
    Bq[[paste(x$section, pg)]] <- data.table(cohort = x$cohort, section = x$section, patient = x$patient,
      role = x$role, program = pg, rule = "composition_quartile", n_hi = sum(hi), n_lo = sum(lo),
      diff = mean(y[hi]) - mean(y[lo])) }
  message("scored ", x$section) }
A <- rbindlist(A); Bq <- rbindlist(Bq); cover <- rbindlist(cover)
rnd <- function(X) { nn <- names(X)[vapply(X, is.numeric, logical(1))]; X[, (nn) := lapply(.SD, round, 3), .SDcols = nn][] }
fwrite(rnd(A), file.path(res_dir, "criterionA_by_section_extra.csv"))
fwrite(rnd(Bq), file.path(res_dir, "criterionB_by_section_extra.csv"))
fwrite(cover, file.path(res_dir, "spatial_gene_coverage_extra.csv"))
cat("\n=== Coverage failures ===\n"); print(cover[panel_present < panel_genes | eval_coverage < 0.8])
cat("\n=== Criterion A by section ===\n")
print(dcast(A[, .(program, section, v = sprintf("%.3f [%.3f, %.3f]", effect, lo, hi))], program ~ section, value.var = "v"))
cat("\n=== Criterion B by section (T/NK composition quartiles) ===\n")
print(dcast(rnd(Bq), program ~ section, value.var = "diff"))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_37.txt"))
 
