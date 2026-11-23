# 01_download_GSE132465.R
# Download the processed GSE132465 supplementary files into 03_original_data/
# and record checksums. Files in 03_original_data are NEVER modified afterwards.
#
# Expected files (confirm on the GEO page):
#   GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz   <- raw counts (use this)
#   GSE132465_GEO_processed_CRC_10X_natural_log_TPM_matrix.txt.gz <- normalized; NOT counts
#   GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz

options(timeout = 60 * 60)   # large files
out_dir <- "03_original_data/GSE132465"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/"
files <- c(
  "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz",
  "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
)

for (f in files) {
  dest <- file.path(out_dir, f)
  if (!file.exists(dest)) {
    message("Downloading ", f)
    download.file(paste0(base, f), dest, mode = "wb")
  } else message("Already present: ", f)
}

# Series-level sample metadata (characteristics, sample titles)
if (requireNamespace("GEOquery", quietly = TRUE)) {
  gse <- GEOquery::getGEO("GSE132465", GSEMatrix = TRUE, destdir = out_dir)
  pd <- Biobase::pData(gse[[1]])
  data.table::fwrite(pd, "02_metadata/GSE132465_geo_sample_metadata.csv")
}

# Checksums: re-run later to prove the originals are unchanged
md5 <- tools::md5sum(list.files(out_dir, full.names = TRUE))
write.csv(data.frame(file = basename(names(md5)), md5 = unname(md5),
                     downloaded = as.character(Sys.Date())),
          "02_metadata/GSE132465_checksums.csv", row.names = FALSE)
print(md5)

dir.create("06_results/sessionInfo", recursive = TRUE, showWarnings = FALSE); writeLines(capture.output(sessionInfo()), "06_results/sessionInfo/sessionInfo_01.txt")
