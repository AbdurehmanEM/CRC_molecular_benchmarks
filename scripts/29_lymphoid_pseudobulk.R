
# 29_lymphoid_pseudobulk.R  (protocol v1.8, section 14) - CD8 compartment pseudobulk
# Mirrors 08_GSE178341_pseudobulk.R exactly: same QC, same primary_fraction rule,
# same aggregation. Existing four-compartment outputs are NOT touched.
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(hdf5r); library(digest) })
set.seed(20261008)
in_dir  <- "03_original_data/GSE178341"
out_dir <- "04_processed_data/GSE178341"
res_dir <- "06_results/GSE178341"
BLOCK <- 10000

meta <- fread(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"))
clus <- fread(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"))
d <- merge(meta, clus, by.x = "cellID", by.y = "sampleID")
d[, primary_fraction := PROCESSING_TYPE %in% c("unsorted", "LiveMACS")]
d[, tumour := SPECIMEN_TYPE == "T"]
d[, is_cd8 := clMidwayPr == "TCD8"]

h5  <- H5File$new(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit.h5"), mode = "r")
bc  <- h5[["matrix/barcodes"]]$read(); gid <- h5[["matrix/features/id"]]$read()
gnm <- h5[["matrix/features/name"]]$read(); ip <- h5[["matrix/indptr"]]$read()
dsx <- h5[["matrix/data"]]; dsi <- h5[["matrix/indices"]]
G <- length(gid); N <- length(bc)
d <- d[match(bc, cellID)]; stopifnot(identical(d$cellID, bc))
mt <- grepl("^MT-", gnm)

read_block <- function(j1, j2) {
  lo <- ip[j1] + 1; hi <- ip[j2 + 1]; p <- ip[j1:(j2 + 1)] - ip[j1]
  if (hi < lo) return(sparseMatrix(i = integer(0), j = integer(0), dims = c(G, j2 - j1 + 1)))
  new("dgCMatrix", i = as.integer(dsi[lo:hi]), p = as.integer(p),
      x = as.numeric(dsx[lo:hi]), Dim = as.integer(c(G, j2 - j1 + 1))) }
starts <- seq(1, N, by = BLOCK)

# ---- pass 1: QC (all tumour cells, so MAD matches script 08 exactly) ----
qcf <- file.path(out_dir, "per_cell_qc_tumour.csv")
if (file.exists(qcf)) { qc <- fread(qcf); stopifnot(identical(qc$cellID, bc)); message("re-using cached QC")
} else {
  qc <- data.table(cellID = bc, nCount = NA_real_, nFeature = NA_real_, pct_mt = NA_real_)
  for (s in starts) { e <- min(s + BLOCK - 1, N)
    if (!any(d$tumour[s:e])) next
    m <- read_block(s, e)
    set(qc, s:e, "nCount", Matrix::colSums(m)); set(qc, s:e, "nFeature", as.numeric(diff(m@p)))
    set(qc, s:e, "pct_mt", 100 * Matrix::colSums(m[mt, , drop = FALSE]) / pmax(Matrix::colSums(m), 1))
    message("QC pass: ", e, " / ", N) }
  fwrite(qc, qcf) }
d <- cbind(d, qc[, .(nCount, nFeature, pct_mt)])
mad_flag <- function(x, n = 3, side = c("both", "high")) { side <- match.arg(side)
  m <- median(x); s <- mad(x); if (!is.finite(s) || s == 0) return(rep(FALSE, length(x)))
  if (side == "high") x > m + n * s else x < m - n * s | x > m + n * s }
d[tumour == TRUE, `:=`(flag_count = mad_flag(log10(nCount + 1)),
  flag_genes = mad_flag(log10(nFeature + 1)),
  flag_mt = mad_flag(pct_mt, side = "high")), by = batchID]
d[tumour == TRUE, qc_pass := !(flag_count | flag_genes | flag_mt)]

# ---- integrity check against script 08 ----
new_ch <- d[tumour == TRUE, .(removed = sum(!qc_pass)), by = .(PID, batchID)]
old <- fread(file.path(res_dir, "qc_per_channel_tumour.csv"))[, .(PID, batchID, removed)]
chk <- merge(old, new_ch, by = c("PID", "batchID"), suffixes = c("_08", "_29"))
if (nrow(chk) != nrow(old) || any(chk$removed_08 != chk$removed_29))
  stop("QC does not reproduce script 08 - STOP, do not use this pseudobulk")
cat("QC reproduces script 08 exactly across", nrow(chk), "channels\n")

# ---- pass 2: aggregate CD8 only ----
ok <- d$tumour & d$qc_pass %in% TRUE
profiles <- list(Lymphoid_CD8         = ok & d$primary_fraction & d$is_cd8,
                 Lymphoid_CD8_allfrac = ok & d$is_cd8)
donors <- sort(unique(d$PID))
acc <- lapply(profiles, function(z) matrix(0, G, length(donors), dimnames = list(gid, donors)))
ncell <- lapply(profiles, function(z) setNames(integer(length(donors)), donors))
for (s in starts) { e <- min(s + BLOCK - 1, N); idx <- s:e
  if (!any(vapply(profiles, function(z) any(z[idx]), logical(1)))) next
  m <- read_block(s, e)
  for (k in names(profiles)) { use <- which(profiles[[k]][idx]); if (!length(use)) next
    fct <- factor(d$PID[idx][use], levels = donors)
    ind <- sparseMatrix(i = seq_along(use), j = as.integer(fct), x = 1, dims = c(length(use), length(donors)))
    acc[[k]] <- acc[[k]] + as.matrix(m[, use, drop = FALSE] %*% ind)
    ncell[[k]] <- ncell[[k]] + tabulate(as.integer(fct), nbins = length(donors)) }
  message("Aggregate pass: ", e, " / ", N) }
h5$close_all()

for (k in names(profiles)) { cnt <- Matrix(acc[[k]], sparse = TRUE)
  stopifnot(all(abs(cnt@x - round(cnt@x)) < 1e-6))
  saveRDS(list(counts = cnt, n_cells = ncell[[k]], profile = k,
               created = as.character(Sys.time()), protocol = "v1.8"),
          file.path(out_dir, paste0("pseudobulk_", k, ".rds"))) }
fwrite(data.table(PID = donors, n_CD8 = ncell[["Lymphoid_CD8"]],
                  n_CD8_allfrac = ncell[["Lymphoid_CD8_allfrac"]]),
       file.path(out_dir, "donor_CD8_cellcounts.csv"))
elig <- rbindlist(lapply(c(10, 20, 30, 50), function(k)
  data.table(min_cells = k, t(vapply(names(profiles), function(p) sum(ncell[[p]] >= k), numeric(1))))))
fwrite(elig, file.path(res_dir, "pseudobulk_eligibility_lymphoid.csv"))
cat("\n=== Eligible donors, lymphoid (after QC) ===\n"); print(elig)
cat("\nLymphoid pseudobulk hash:", digest(lapply(acc, function(a) round(colSums(a)))), "\n")
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_29.txt"))
 
