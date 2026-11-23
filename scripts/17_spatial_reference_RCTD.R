
# 17_spatial_reference_RCTD.R  (protocol v1.2, section 10d) - RCTD reference from GSE178341 (composition only)
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(hdf5r) })
set.seed(20261012)
dst <- "04_processed_data/Valdeolivas_Visium/unzipped"
secs <- grep("^SN", list.dirs(dst, recursive = FALSE, full.names = FALSE), value = TRUE)
have_h5 <- file.exists(file.path(dst, secs, "filtered_feature_bc_matrix.h5"))
cat("Sections found:", length(secs), "| with h5:", sum(have_h5), "\n")
if (length(secs) < 14 || !all(have_h5)) stop("Some sections missing - rerun the unzip loop")
# 19 Oct 2026: the GitHub fallback install (dmcable/spacexr) was removed. That version
# provides create.RCTD()/run.RCTD(), whereas scripts 18 and 22 call the Bioconductor
# API createRctd()/runRctd(); on a clean machine stage 18 would have failed.
if (!requireNamespace("spacexr", quietly = TRUE) || !exists("createRctd", envir = asNamespace("spacexr")))
  stop("Bioconductor spacexr with createRctd() is required - run 00_setup_and_versions.R")
cat("spacexr version:", as.character(packageVersion("spacexr")), "\n")
# FIX_SELECTION_ORDER: the original code reorders the 1,200 sampled cells per type by
# h5 position and then keeps the FIRST 800 that pass QC, so donors stored early in the
# h5 file are over-represented. FALSE reproduces the published reference exactly;
# TRUE keeps 800 passing cells in random order and writes a SEPARATE file
# (..._randomsel.rds) for a sensitivity analysis - the published file is not touched.
if (!exists("FIX_SELECTION_ORDER")) FIX_SELECTION_ORDER <- FALSE
in_dir <- "03_original_data/GSE178341"; out_dir <- "04_processed_data/spatial_reference"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
meta <- fread(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"))
clus <- fread(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"))
d <- merge(meta, clus, by.x = "cellID", by.y = "sampleID")
cat("\nclMidwayPr labels:\n"); print(d[, .N, by = clMidwayPr][order(-N)])
Tm <- d$SPECIMEN_TYPE == "T"
d[, ref_type := fcase(
  clMidwayPr == "Epi" & Tm & cl295v11SubShort %in% c("cE01", "cE02", "cE03"), "Malignant_epithelial",
  clMidwayPr == "Epi" & SPECIMEN_TYPE == "N", "Normal_epithelial",
  clMidwayPr == "Fibro" & Tm, "Fibroblast",
  clMidwayPr %in% c("Peri", "SmoothMuscle") & Tm, "Pericyte_SMC",
  clMidwayPr == "Endo" & Tm, "Endothelial",
  clMidwayPr %in% c("Mono", "Macro", "DC", "Granulo") & Tm, "Myeloid",
  grepl("Mast", clMidwayPr) & Tm, "Mast",
  clMidwayPr %in% c("B", "Plasma") & Tm, "B_Plasma",
  grepl("^T|NK|ILC", clMidwayPr) & Tm, "T_NK", default = NA_character_)]
cand <- d[!is.na(ref_type), .SD[sample(.N, min(.N, 1200))], by = ref_type]
h5 <- H5File$new(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit.h5"), mode = "r")
bc <- h5[["matrix/barcodes"]]$read(); gnm <- h5[["matrix/features/name"]]$read()
ip <- h5[["matrix/indptr"]]$read(); dsx <- h5[["matrix/data"]]; dsi <- h5[["matrix/indices"]]
idx <- match(cand$cellID, bc); stopifnot(!anyNA(idx))
ord <- order(idx); idx <- idx[ord]; cand <- cand[ord]
h5$close_all(); h5path <- file.path(in_dir, 'GSE178341_crc10x_full_c295v4_submit.h5')
WIN <- 40000; ii <- vector('list', length(idx)); xx <- ii
for (w0 in seq(1, length(ip) - 1, by = WIN)) { w1 <- min(w0 + WIN - 1, length(ip) - 1)
  k <- which(idx >= w0 & idx <= w1); if (!length(k)) next
  j1 <- min(idx[k]); j2 <- max(idx[k]); st <- ip[j1] + 1; n <- ip[j2 + 1] - ip[j1]
  I <- as.integer(rhdf5::h5read(h5path, 'matrix/indices', start = st, count = n))
  X <- as.numeric(rhdf5::h5read(h5path, 'matrix/data', start = st, count = n)); off <- ip[j1]
  for (kk in k) { j <- idx[kk]; s <- ip[j] - off + 1; e <- ip[j + 1] - off
    if (e >= s) { ii[[kk]] <- I[s:e]; xx[[kk]] <- X[s:e] } else { ii[[kk]] <- integer(0); xx[[kk]] <- numeric(0) } }
  rm(I, X); invisible(gc()); message('cells ', w0, '-', w1, ' done') }
cnt <- new("dgCMatrix", i = as.integer(unlist(ii)), p = as.integer(c(0, cumsum(lengths(ii)))),
  x = as.numeric(unlist(xx)), Dim = as.integer(c(length(gnm), length(idx))))
dimnames(cnt) <- list(make.unique(gnm), cand$cellID)
mt <- grepl("^MT-", rownames(cnt))
cand[, `:=`(nCount = Matrix::colSums(cnt), nFeature = as.numeric(diff(cnt@p)),
  pct_mt = 100 * Matrix::colSums(cnt[mt, , drop = FALSE]) / pmax(Matrix::colSums(cnt), 1))]
cand[, pass := nCount >= 500 & nFeature >= 200 & pct_mt < 25]
if (FIX_SELECTION_ORDER) {
  cand[, rk := { o <- sample(.N); r <- integer(.N); r[o] <- cumsum(pass[o]); r }, by = ref_type]
  out_file <- "GSE178341_reference_counts_randomsel.rds"
} else {
  cand[, rk := cumsum(pass), by = ref_type]            # published behaviour (h5 order)
  out_file <- "GSE178341_reference_counts.rds"
}
sel <- cand[pass & rk <= 800]
cnt <- cnt[, sel$cellID]
saveRDS(list(counts = cnt, cell_types = setNames(sel$ref_type, sel$cellID), meta = sel[, .(cellID, PID, ref_type, clMidwayPr, cl295v11SubShort)]),
  file.path(out_dir, out_file))
cat("\nSaved", out_file, "\n")
cat("\n=== Reference cells per type (donors) ===\n")
print(sel[, .(cells = .N, donors = uniqueN(PID), median_UMI = median(nCount)), by = ref_type][order(ref_type)])
 
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_17.txt"))
