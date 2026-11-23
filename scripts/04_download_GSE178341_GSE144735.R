# 04_download_GSE178341_GSE144735.R
# Download the processed supplementary files for the discovery cohort (GSE178341,
# Pelka 2021) and the KUL3 validation cohort (GSE144735, Lee 2020) into
# 03_original_data/, and record checksums. Files in 03_original_data are NEVER
# modified afterwards.
#
# RECONSTRUCTED 19 Oct 2026. The project copy of this script had been replaced by
# an editor helper (setwd + file.create + file.edit). file.create() TRUNCATES an
# existing file, so sourcing that helper erased 06_audit_GSE178341.R. The file
# names below are the ones read by scripts 06, 07, 08, 12/12b, 17 and 29.
#
# If checksums from the original download exist in 02_metadata/, the files on disk
# are compared against them (originals-unchanged check) instead of overwriting.

options(timeout = 6 * 60 * 60)
dir.create("02_metadata", showWarnings = FALSE)

FILES <- list(
  GSE178341 = c("GSE178341_crc10x_full_c295v4_submit.h5",
                "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz",
                "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"),
  GSE144735 = c("GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz",
                "GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz"))

geo_base <- function(acc) sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/suppl/",
                                  substr(acc, 1, nchar(acc) - 3), acc)

for (acc in names(FILES)) {
  out_dir <- file.path("03_original_data", acc)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (f in FILES[[acc]]) {
    dest <- file.path(out_dir, f)
    if (file.exists(dest)) { message("Already present: ", f); next }
    message("Downloading ", f)
    ok <- tryCatch({ download.file(paste0(geo_base(acc), f), dest, mode = "wb"); TRUE },
                   error = function(e) { message("  FAILED: ", conditionMessage(e)); FALSE })
    if (!ok) { if (file.exists(dest)) file.remove(dest); stop("Download failed: ", acc, "/", f) }
  }

  fl   <- file.path(out_dir, FILES[[acc]])
  info <- data.frame(file = FILES[[acc]], size_MB = round(file.size(fl) / 1e6, 1),
                     md5 = unname(tools::md5sum(fl)), checked = as.character(Sys.Date()))
  ck <- file.path("02_metadata", paste0(acc, "_checksums.csv"))
  if (file.exists(ck)) {
    old <- read.csv(ck, stringsAsFactors = FALSE)
    cmp <- merge(old[, c("file", "md5")], info[, c("file", "md5")], by = "file",
                 suffixes = c("_recorded", "_now"), all.y = TRUE)
    cmp$unchanged <- cmp$md5_recorded == cmp$md5_now
    print(cmp)
    if (any(cmp$unchanged %in% FALSE)) stop(acc, ": a file differs from its recorded checksum - STOP")
    write.csv(info, file.path("02_metadata", paste0(acc, "_checksums_recheck.csv")), row.names = FALSE)
    message(acc, ": all files match the recorded checksums")
  } else {
    write.csv(info, ck, row.names = FALSE)
    print(info)
  }
}

# Series-level sample metadata, as in script 01
if (requireNamespace("GEOquery", quietly = TRUE)) for (acc in names(FILES)) {
  gse <- tryCatch(GEOquery::getGEO(acc, GSEMatrix = TRUE, destdir = file.path("03_original_data", acc)),
                  error = function(e) NULL)
  if (!is.null(gse)) data.table::fwrite(Biobase::pData(gse[[1]]),
                                        file.path("02_metadata", paste0(acc, "_geo_sample_metadata.csv")))
}
