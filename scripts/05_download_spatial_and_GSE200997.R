# 05_download_spatial_and_GSE200997.R
# Protocol v0.4: Valdeolivas Visium (Zenodo, main spatial), GSE280318 10x CRC study
# (Visium/Visium HD/Xenium + matched single-cell; subseries GSE280311/13/15), GSE226997
# (listed only), and scRNA-seq validation cohorts GSE200997 and GSE166555.
# Run step by step. Files go to 03_original_data/ and are never modified.

suppressPackageStartupMessages({ library(GEOquery); library(jsonlite) })
options(timeout = 6 * 60 * 60)
dir.create("02_metadata", showWarnings = FALSE)

save_checksums <- function(dir, tag) {
  fl <- list.files(dir, full.names = TRUE, recursive = TRUE)
  if (!length(fl)) return(invisible())
  info <- data.frame(file = sub(paste0(dir, "/"), "", fl), size_MB = round(file.size(fl) / 1e6, 1),
                     md5 = unname(tools::md5sum(fl)), downloaded = as.character(Sys.Date()))
  write.csv(info, file.path("02_metadata", paste0(tag, "_checksums.csv")), row.names = FALSE)
  print(info[, c("file", "size_MB")])
}

safe_download <- function(url, dest) {
  if (file.exists(dest)) { message("Already present: ", basename(dest)); return(TRUE) }
  message("Downloading ", basename(dest), " ...")
  ok <- tryCatch({ download.file(url, dest, mode = "wb"); TRUE },
                 error = function(e) { message("  FAILED: ", conditionMessage(e)); FALSE })
  if (!ok && file.exists(dest)) file.remove(dest)
  ok
}

# ================================================================ 1. Valdeolivas Visium (Zenodo)
# Latest version of record 7739700. Only Space Ranger outputs and pathologist
# annotations are downloaded; the authors' derived files (deconvolution,
# intermediary objects; 4-5 GB zips) are skipped - we compute our own.
rec <- fromJSON("https://zenodo.org/api/records/7739700/versions/latest")
cat("Zenodo record:", rec$id, "| version:", rec$metadata$version, "|", rec$metadata$title, "\n")
files <- rec$files
files$size_MB <- round(files$size / 1e6, 1)
print(files[, c("key", "size_MB")])
write.csv(files[, c("key", "size_MB")], "02_metadata/Valdeolivas_zenodo_listing.csv", row.names = FALSE)

keep <- !grepl("Intermediary|Deconvolution|LiverMetastasis", files$key, ignore.case = TRUE)
out_dir <- "03_original_data/Valdeolivas_Visium"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (k in which(keep)) safe_download(files$links$self[k], file.path(out_dir, files$key[k]))
save_checksums(out_dir, "Valdeolivas_Visium")

# ================================================================ 2. GEO cohorts
geo <- c("GSE200997", "GSE166555",                       # scRNA-seq validation
         "GSE280318", "GSE280311", "GSE280313", "GSE280315",  # 10x CRC study (super + subseries)
         "GSE226997")                                     # listed only (41 GB archive)
MAX_GB <- 3.9                                             # FAT32 limit on F: drive
remote_size_GB <- function(url) {
  h <- tryCatch(curlGetHeaders(url), error = function(e) character(0))
  cl <- grep("^content-length:", h, ignore.case = TRUE, value = TRUE)
  if (!length(cl)) return(NA_real_)
  as.numeric(sub(".*:\\s*", "", tail(cl, 1))) / 1e9
}
listing <- do.call(rbind, lapply(geo, function(acc) {
  x <- tryCatch(getGEOSuppFiles(acc, makeDirectory = FALSE, fetch_files = FALSE),
                error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  data.frame(accession = acc, fname = x$fname, url = x$url)
}))
listing$size_GB <- round(vapply(listing$url, remote_size_GB, numeric(1)), 2)
print(listing[, c("accession", "fname", "size_GB")])
write.csv(listing, "02_metadata/geo_v04_supplementary_listing.csv", row.names = FALSE)

# GSE226997 is a single 41 GB RAW.tar - too big for the F: drive (FAT32, 4 GB per file).
# It is LISTED here but NOT downloaded; see the note printed below.
for (k in seq_len(nrow(listing))) {
  # RAW.tar archives are kept when small (some cohorts ship their per-sample count
  # matrices only inside RAW.tar); the size check below skips anything > 3.9 GB.
  if (listing$accession[k] == "GSE226997") {
    message("NOT downloading ", listing$fname[k], " (41 GB archive; see note).")
    next
  }
  if (!is.na(listing$size_GB[k]) && listing$size_GB[k] > MAX_GB) {
    message("NOT downloading ", listing$fname[k], " (", listing$size_GB[k],
            " GB > 4 GB FAT32 limit; needs C: or an NTFS drive).")
    next
  }
  d <- file.path("03_original_data", listing$accession[k])
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  safe_download(listing$url[k], file.path(d, listing$fname[k]))
}
for (acc in geo) {
  gse <- tryCatch(getGEO(acc, GSEMatrix = TRUE, destdir = file.path("03_original_data", acc)),
                  error = function(e) NULL)
  if (!is.null(gse)) data.table::fwrite(Biobase::pData(gse[[1]]),
                                        file.path("02_metadata", paste0(acc, "_geo_sample_metadata.csv")))
  save_checksums(file.path("03_original_data", acc), acc)
}
message("NOTE: GSE226997 per-sample files can be fetched individually from each GSM page ",
        "(GSM7089855-58) if they are provided there; check before downloading the 41 GB archive.")

dir.create("06_results/sessionInfo", recursive = TRUE, showWarnings = FALSE); writeLines(capture.output(sessionInfo()), "06_results/sessionInfo/sessionInfo_05.txt")
