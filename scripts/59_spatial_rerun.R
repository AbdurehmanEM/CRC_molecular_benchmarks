###############################################################################
# 59_spatial_rerun.R - Reviewer points 9 and 10 (pseudoreplication, circularity)
#
# Re-implements criterion A at spot level so that it can be varied:
#   (a) block size of the spatial block bootstrap: 4, 8 (published), 12, 16 spots;
#   (b) negative-control targets: random non-program genes matched to the
#       evaluation genes on within-section expression -> null distribution of
#       the composition-adjusted panel coefficient;
#   (c) residualisation: criterion A with log total UMI and a matched
#       negative-control gene score added as covariates;
#   (d) alternative RCTD weights (independent reference from script 60, or the
#       random-selection reference from script 42) via WEIGHT_SOURCES.
# GATE: with block size 8 and the published weights, patient coefficients must
# match Supplementary Table S10a (tolerance 0.01). If they do not, set the
# SPATIAL_* options in the config to match script 19/22 before reading results.
#
# SECTION 0 (run once): export one .rds per section from your script 19/22
# objects into SPATIAL_SECTION_DIR with this structure:
#   list(section = "SN048_A416371_Rep1", patient = "A416371", cohort = "Valdeolivas",
#        grid = "hex" | "square",          # Visium = hex, Visium HD pseudo-spots = square
#        counts = <dgCMatrix genes x spots, rownames = gene symbols, QC-passed spots>,
#        coords = data.frame(spot, array_row, array_col),
#        weights = list(published = <spots x 9 matrix, RCTD proportions, colnames = cell types>,
#                       randomsel = <optional>, independent = <added by script 60>))
# A template is given in export_section_template() below; adapt the three
# marked lines to the object names used in your scripts.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision()

export_section_template <- function() {
  # ---- ADAPT: loop over the objects your script 19 / 22 creates --------------
  # for (sec in sections) {
  #   obj <- readRDS(<your section object>)                      # ADAPT
  #   w   <- readRDS(<your RCTD weights for this section>)       # ADAPT (normalised proportions)
  #   saveRDS(list(section = sec, patient = <id>, cohort = <cohort>, grid = "hex",
  #                counts = obj$counts, coords = obj$coords,     # ADAPT
  #                weights = list(published = w)),
  #           file.path(SPATIAL_SECTION_DIR, paste0(sec, ".rds")))
  # }
  invisible(NULL)
}

WEIGHT_SOURCES <- c("published", "randomsel", "independent")   # those absent in a file are skipped

panels <- load_panels(); partition <- load_partition()
files <- list.files(SPATIAL_SECTION_DIR, pattern = "\\.rds$", full.names = TRUE)
if (!length(files)) stop("No section files in ", SPATIAL_SECTION_DIR, ". Run SECTION 0 (export) first.")

## ---- spot-level helpers ---------------------------------------------------------
spot_logcpm <- function(counts, genes, lib) {
  x <- as.matrix(counts[genes, , drop = FALSE])
  if (SPATIAL_LOGNORM == "log1p_cp10k") return(log1p(t(t(x) / lib) * 1e4))
  adj <- PRIOR_COUNT * lib / mean(lib)
  t(log2(t(x) + adj) - log2(lib + 2 * adj) + log2(1e6))
}
zrows <- function(M) { s <- apply(M, 1, sd); ok <- s > 0; (M[ok, , drop = FALSE] - rowMeans(M[ok, , drop = FALSE])) / s[ok] }
stdz <- function(v) (v - mean(v)) / sd(v)
section_scores <- function(sec, eval_sym, panel) {
  lib <- colSums(sec$counts)
  ev <- intersect(eval_sym, rownames(sec$counts)); cover <- length(ev) / length(eval_sym)
  if (!all(panel$symbol %in% rownames(sec$counts)) || cover < MIN_EVAL_COVER) return(NULL)
  tz <- zrows(spot_logcpm(sec$counts, ev, lib))
  pz <- spot_logcpm(sec$counts, panel$symbol, lib); pz <- (pz - rowMeans(pz)) / pmax(apply(pz, 1, sd), 1e-8)
  list(target = stdz(colMeans(tz)), panel = stdz(colSums(pz * panel$coef)), lib = lib, cover = cover)
}
design <- function(sc, W, extra = NULL) {
  W <- W[, setdiff(colnames(W), SPATIAL_DROP_TYPE), drop = FALSE]
  cbind(1, panel = sc$panel, W, extra)
}
coefA <- function(X, y) .lm.fit(X, y)$coefficients[2]
blocks_of <- function(coords, grid, bs) {
  cs <- if (grid == "hex") 2 * bs else bs           # hex array_col steps by 2
  paste(floor(coords$array_row / bs), floor(coords$array_col / cs))
}
block_boot <- function(X, y, blocks, B) {
  ub <- unique(blocks); idx <- split(seq_along(blocks), blocks)
  replicate(B, { i <- unlist(idx[sample(ub, length(ub), replace = TRUE)], use.names = FALSE); coefA(X[i, , drop = FALSE], y[i]) })
}

## ---- load sections -----------------------------------------------------------
secs <- lapply(files, readRDS)
secs <- secs[!vapply(secs, function(s) s$patient %in% SPATIAL_EXCLUDE_PATIENTS, logical(1))]
msg(length(secs), "sections from", length(unique(sapply(secs, `[[`, "patient"))), "patients")

set.seed(SEED)
res <- list(); nulls <- list()
for (p in names(panels)) {
  panel <- panels[[p]]; eval_sym <- partition[program == p & role == "evaluation", symbol]
  all_prog <- unique(partition$symbol)
  for (s in secs) {
    sc <- section_scores(s, eval_sym, panel); if (is.null(sc)) next
    spots <- colnames(s$counts)
    # matched negative-control gene sets (within-section mean expression deciles)
    nonprog <- setdiff(rownames(s$counts), all_prog)
    m_all <- Matrix::rowMeans(s$counts[c(intersect(eval_sym, rownames(s$counts)), nonprog), , drop = FALSE] > 0)
    dec <- cut(m_all, unique(quantile(m_all, 0:10 / 10)), include.lowest = TRUE, labels = FALSE); names(dec) <- names(m_all)
    draw_nc <- function() { ev <- intersect(eval_sym, rownames(s$counts))
      vapply(ev, function(g) { c <- nonprog[dec[nonprog] == dec[[g]]]; if (!length(c)) c <- nonprog; c[sample.int(length(c), 1)] }, "") }
    nc_score <- function(g) stdz(colMeans(zrows(spot_logcpm(s$counts, unique(g), sc$lib))))
    for (ws in intersect(WEIGHT_SOURCES, names(s$weights))) {
      W <- as.matrix(s$weights[[ws]]); keep <- intersect(spots, rownames(W))
      if (length(keep) < 50) next
      scw <- list(target = sc$target[keep], panel = sc$panel[keep], lib = sc$lib[keep])
      W <- W[keep, , drop = FALSE]
      X <- design(scw, W); y <- scw$target
      crd <- s$coords[match(keep, s$coords$spot), ]
      est <- coefA(X, y)
      for (bs in SPATIAL_BLOCK_SIZES) {
        if (ws != "published" && bs != 8) next
        bt <- block_boot(X, y, blocks_of(crd, s$grid %||% "hex", bs), SPATIAL_N_BOOT)
        res[[length(res) + 1]] <- data.table(program = p, cohort = s$cohort, patient = s$patient, section = s$section,
                                             weights = ws, model = "criterion A", block = bs, est = est, boot = list(bt))
      }
      if (ws == "published") {
        ## residualised model: + log total UMI + one matched negative-control score
        ncs <- nc_score(draw_nc())[keep]
        Xr <- design(scw, W, cbind(loglib = log(scw$lib), negctrl = ncs))
        bt <- block_boot(Xr, y, blocks_of(crd, s$grid %||% "hex", 8), SPATIAL_N_BOOT)
        res[[length(res) + 1]] <- data.table(program = p, cohort = s$cohort, patient = s$patient, section = s$section,
                                             weights = ws, model = "criterion A + log UMI + negative-control score",
                                             block = 8, est = coefA(Xr, y), boot = list(bt))
        ## negative-control targets: same panel, random matched target genes
        nl <- replicate(SPATIAL_N_NEGCTRL, { yn <- nc_score(draw_nc())[keep]; coefA(X, yn) })
        nulls[[length(nulls) + 1]] <- data.table(program = p, patient = s$patient, section = s$section, draw = seq_along(nl), null_est = nl)
      }
    }
  }
  msg("done", p)
}
res <- rbindlist(res); nulls <- rbindlist(nulls)

## ---- patient level: average sections; SE = SD of averaged bootstrap replicates
pat <- res[, { bm <- Reduce(`+`, boot) / length(boot); .(est = mean(est), se = sd(bm), n_sections = .N) },
           by = .(program, cohort, patient, weights, model, block)]
## ---- gate against S10a (published weights, block 8, criterion A)
A <- read_supp("S10a_spatial_criterionA")
setnames(A, c("Program code", "Patient", "Criterion A coefficient", "Bootstrap SE"), c("program", "patient", "coef_pub", "se_pub"))
g <- merge(pat[weights == "published" & model == "criterion A" & block == 8], A[, .(program, patient, coef_pub = as.numeric(coef_pub), se_pub = as.numeric(se_pub))],
           by = c("program", "patient"))
g[, `:=`(d_coef = round(est - coef_pub, 4), d_se = round(se - se_pub, 4))]
write_out(g, "59_spatial_gate")
gate_ok <- nrow(g) > 0 && all(abs(g$d_coef) <= 0.01)
msg("Spatial gate:", if (gate_ok) "PASS" else "FAIL (adjust SPATIAL_* options to match script 19/22)")

## ---- pooled summaries per scenario
pooled <- pat[, { d <- dl_pool(est, pmax(se, 1e-6)); thr <- ceiling(2/3 * .N - 1e-9)
  .(n_patients = .N, pos = sum(est > 0), pooled = d$est, lo = d$lo, hi = d$hi, I2 = d$I2, metA = sum(est > 0) >= thr) },
  by = .(program, weights, model, block)]
## negative-control comparison: observed patient-mean coefficient vs null patient-mean distribution
obs <- pat[weights == "published" & model == "criterion A" & block == 8, .(obs_mean = mean(est)), by = program]
nm <- nulls[, .(null_est = mean(null_est)), by = .(program, patient, draw)][, .(null_mean = mean(null_est)), by = .(program, draw)]
negc <- nm[obs, on = "program"][, .(obs_mean = obs_mean[1], null_median = median(null_mean),
                                    null_q975 = quantile(null_mean, 0.975), p_empirical = (1 + sum(null_mean >= obs_mean[1])) / (1 + .N)), by = program]
num <- c("pooled", "lo", "hi", "I2"); pooled[, (num) := lapply(.SD, round, 3), .SDcols = num]
write_out(pooled, "59_spatial_rerun_pooled")
write_out(pat[, .(program, cohort, patient, weights, model, block, est = round(est, 4), se = round(se, 4), n_sections)], "59_spatial_rerun_patients")
write_out(negc, "59_spatial_negative_control")
print(pooled); print(negc)
