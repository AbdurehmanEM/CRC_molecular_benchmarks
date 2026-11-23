# 08_GSE178341_pseudobulk.R
# Discovery cohort: donor x compartment pseudobulk raw-count matrices (protocol v0.5).
# Reads the h5 in blocks of cells, so the full ~9 GB matrix is never in memory.
#   Pass 1: per-cell QC metrics -> per-channel 3-MAD flags
#   Pass 2: sum raw counts of retained cells per donor for each profile

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(hdf5r); library(digest) })
set.seed(20261008)
in_dir  <- "03_original_data/GSE178341"
out_dir <- "04_processed_data/GSE178341"
res_dir <- "06_results/GSE178341"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
BLOCK <- 10000                        # lower to 5000 if memory is short

meta <- fread(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"))
clus <- fread(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"))
d <- merge(meta, clus, by.x = "cellID", by.y = "sampleID")
d[, compartment := fcase(
  clMidwayPr == "Endo" & grepl("lymphatic", cl295v11SubFull, ignore.case = TRUE), "Lymphatic_EC",
  clMidwayPr == "Endo",                              "Endothelial",
  clMidwayPr == "Fibro",                             "Stromal",
  clMidwayPr %in% c("Peri", "SmoothMuscle"),         "Mural_other",
  clMidwayPr == "Schwann",                           "Glia",
  clMidwayPr %in% c("Mono", "Macro"),                "Myeloid",
  clMidwayPr %in% c("DC", "Granulo"),                "Myeloid_other",
  clMidwayPr == "Epi",                               "Epithelial",
  default = "Other_context")]
d[, stemTA := cl295v11SubShort %in% c("cE01", "cE02", "cE03")]
d[, primary_fraction := PROCESSING_TYPE %in% c("unsorted", "LiveMACS")]
d[, tumour := SPECIMEN_TYPE == "T"]

h5   <- H5File$new(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit.h5"), mode = "r")
bc   <- h5[["matrix/barcodes"]]$read()
gid  <- h5[["matrix/features/id"]]$read()
gnm  <- h5[["matrix/features/name"]]$read()
ip   <- h5[["matrix/indptr"]]$read()
dsx  <- h5[["matrix/data"]]; dsi <- h5[["matrix/indices"]]
G <- length(gid); N <- length(bc)
d <- d[match(bc, cellID)]; stopifnot(identical(d$cellID, bc))
mt <- grepl("^MT-", gnm)
fwrite(data.table(ensembl_id = gid, symbol = gnm, mito = mt), file.path(out_dir, "genes.csv"))

read_block <- function(j1, j2) {
  lo <- ip[j1] + 1; hi <- ip[j2 + 1]
  p  <- ip[j1:(j2 + 1)] - ip[j1]
  if (hi < lo) return(sparseMatrix(i = integer(0), j = integer(0), dims = c(G, j2 - j1 + 1)))
  new("dgCMatrix", i = as.integer(dsi[lo:hi]), p = as.integer(p),
      x = as.numeric(dsx[lo:hi]), Dim = as.integer(c(G, j2 - j1 + 1)))
}
starts <- seq(1, N, by = BLOCK)
need_any <- d$tumour

# ---- pass 1: QC metrics
qc <- data.table(cellID = bc, nCount = NA_real_, nFeature = NA_real_, pct_mt = NA_real_)
for (s in starts) {
  e <- min(s + BLOCK - 1, N)
  if (!any(need_any[s:e])) next
  m <- read_block(s, e)
  set(qc, s:e, "nCount",   Matrix::colSums(m))
  set(qc, s:e, "nFeature", as.numeric(diff(m@p)))
  set(qc, s:e, "pct_mt",   100 * Matrix::colSums(m[mt, , drop = FALSE]) / pmax(Matrix::colSums(m), 1))
  message("QC pass: cells ", e, " / ", N)
}
d <- cbind(d, qc[, .(nCount, nFeature, pct_mt)])
mad_flag <- function(x, n = 3, side = c("both", "high")) {
  side <- match.arg(side); m <- median(x); s <- mad(x)
  if (!is.finite(s) || s == 0) return(rep(FALSE, length(x)))
  if (side == "high") x > m + n * s else x < m - n * s | x > m + n * s
}
d[tumour == TRUE, `:=`(
  flag_count = mad_flag(log10(nCount + 1)),
  flag_genes = mad_flag(log10(nFeature + 1)),
  flag_mt    = mad_flag(pct_mt, side = "high")), by = batchID]
d[tumour == TRUE, qc_pass := !(flag_count | flag_genes | flag_mt)]
qc_ch <- d[tumour == TRUE, .(cells = .N, removed = sum(!qc_pass),
                             pct_removed = round(100 * mean(!qc_pass), 1),
                             median_counts = median(nCount), median_genes = median(nFeature),
                             median_pct_mt = round(median(pct_mt), 2)),
           by = .(PID, batchID, PROCESSING_TYPE, SINGLECELL_TYPE)]
fwrite(qc_ch, file.path(res_dir, "qc_per_channel_tumour.csv"))
cat("\nCells removed by QC flags (tumour):", sum(!d$qc_pass, na.rm = TRUE), "of", sum(d$tumour), "\n")

# ---- profiles (primary + sensitivity)
ok <- d$tumour & d$qc_pass %in% TRUE
profiles <- list(
  Epithelial      = ok & d$primary_fraction & d$compartment == "Epithelial" & d$stemTA,
  Stromal         = ok & d$primary_fraction & d$compartment == "Stromal",
  Myeloid         = ok & d$primary_fraction & d$compartment == "Myeloid",
  Endothelial     = ok & d$primary_fraction & d$compartment == "Endothelial",
  Epithelial_all  = ok & d$primary_fraction & d$compartment == "Epithelial",
  Epithelial_allfrac  = ok & d$compartment == "Epithelial" & d$stemTA,
  Stromal_allfrac     = ok & d$compartment == "Stromal",
  Myeloid_allfrac     = ok & d$compartment == "Myeloid",
  Endothelial_allfrac = ok & d$compartment == "Endothelial")
donors <- sort(unique(d$PID))
acc   <- lapply(profiles, function(z) matrix(0, G, length(donors), dimnames = list(gid, donors)))
ncell <- lapply(profiles, function(z) setNames(integer(length(donors)), donors))

# ---- pass 2: aggregate
for (s in starts) {
  e <- min(s + BLOCK - 1, N)
  idx <- s:e
  if (!any(vapply(profiles, function(z) any(z[idx]), logical(1)))) next
  m <- read_block(s, e)
  for (k in names(profiles)) {
    use <- which(profiles[[k]][idx])
    if (!length(use)) next
    f <- factor(d$PID[idx][use], levels = donors)
    ind <- sparseMatrix(i = seq_along(use), j = as.integer(f), x = 1,
                        dims = c(length(use), length(donors)))
    acc[[k]] <- acc[[k]] + as.matrix(m[, use, drop = FALSE] %*% ind)
    ncell[[k]] <- ncell[[k]] + tabulate(as.integer(f), nbins = length(donors))
  }
  message("Aggregate pass: cells ", e, " / ", N)
}
h5$close_all()

# ---- donor metadata + save
dm <- d[tumour == TRUE, .(hospital = SOURCE_HOSPITAL[1], MMRStatus = MMRStatus[1],
                          site = TissueSiteSimple[1], metastasis = MetastasisStatus[1],
                          sex = Sex[1], age = Age[1],
                          n_tumour_specimens = uniqueN(PatientTypeID),
                          chemistry = names(which.max(table(SINGLECELL_TYPE))),
                          has_primary_fraction = any(primary_fraction)), by = PID]
for (k in names(profiles)) {
  dm[[paste0("n_", k)]] <- ncell[[k]][dm$PID]
  cnt <- Matrix(acc[[k]], sparse = TRUE)
  stopifnot(all(abs(cnt@x - round(cnt@x)) < 1e-6))
  saveRDS(list(counts = cnt, n_cells = ncell[[k]], profile = k,
               created = as.character(Sys.time()), protocol = "v0.5"),
          file.path(out_dir, paste0("pseudobulk_", k, ".rds")))
}
fwrite(dm, file.path(out_dir, "donor_metadata_tumour.csv"))

elig <- rbindlist(lapply(c(10, 20, 30, 50), function(k)
  data.table(min_cells = k, t(vapply(names(profiles), function(p) sum(ncell[[p]] >= k), numeric(1))))))
fwrite(elig, file.path(res_dir, "pseudobulk_eligibility.csv"))
cat("\n=== Eligible donors per profile (after QC) ===\n"); print(elig)

h <- digest(lapply(acc, function(a) round(colSums(a))))
cat("\nPseudobulk hash:", h, "\n")
cat(format(Sys.time()), h, "\n", file = file.path(res_dir, "pseudobulk_hashes.txt"), append = TRUE)
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_08.txt"))