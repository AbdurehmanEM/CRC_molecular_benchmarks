# 02_pilot_load_GSE132465.R
# Load the GSE132465 raw UMI matrix as a SPARSE matrix, match it to the cell
# annotation, and check the pilot go/no-go conditions.

suppressPackageStartupMessages({ library(data.table); library(Matrix) })
set.seed(20261008)

in_dir   <- "03_original_data/GSE132465"
out_dir  <- "04_processed_data/GSE132465"
res_dir  <- "06_results/pilot"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

count_file <- file.path(in_dir, "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz")
annot_file <- file.path(in_dir, "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz")

# ---------------------------------------------------------------- annotation
annot <- fread(annot_file)
cat("Annotation columns:", paste(names(annot), collapse = ", "), "\n")
print(head(annot))

# >>> EDIT THESE if the printed column names differ <<<
COL_CELL    <- "Index"
COL_PATIENT <- "Patient"
COL_CLASS   <- "Class"
COL_SAMPLE  <- "Sample"
COL_TYPE    <- "Cell_type"
COL_SUBTYPE <- "Cell_subtype"
stopifnot(all(c(COL_CELL, COL_PATIENT, COL_CLASS, COL_SAMPLE, COL_TYPE, COL_SUBTYPE) %in% names(annot)))

# ---------------------------------------------------------------- pilot donors
donor_class <- annot[, .(classes = paste(sort(unique(get(COL_CLASS))), collapse = "+"),
                         n_cells = .N), by = COL_PATIENT]
print(donor_class)
paired <- donor_class[grepl("\\+", classes)][[COL_PATIENT]]
PILOT_DONORS <- head(unique(c(paired, donor_class[[COL_PATIENT]])), 4)
cat("Pilot donors:", paste(PILOT_DONORS, collapse = ", "), "\n")
FULL_RUN <- TRUE   # set TRUE after the pilot passes to load all donors
keep_cells <- if (FULL_RUN) annot[[COL_CELL]] else annot[get(COL_PATIENT) %in% PILOT_DONORS][[COL_CELL]]

# ---------------------------------------------------------------- counts (streamed, sparse)
# Read the .gz directly, one pass, 2000 genes at a time. No decompressed copy is
# written to disk (the unzipped matrix is >4 GB, which fails on FAT32/USB drives).
con <- gzfile(count_file, open = "r")
header <- gsub("\"", "", strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]])
first  <- readLines(con, n = 1)
n_fields <- length(strsplit(first, "\t", fixed = TRUE)[[1]])
cell_ids <- if (length(header) == n_fields) header[-1] else header
stopifnot(length(cell_ids) == n_fields - 1)
sel <- which(cell_ids %in% keep_cells)
cat("Cells in matrix:", length(cell_ids), " | selected:", length(sel), "\n")

chunk <- 2000
blocks <- list(); genes <- character(0); i <- 0
pending <- first
repeat {
  lines <- c(pending, readLines(con, n = chunk - length(pending)))
  pending <- character(0)
  if (length(lines) == 0) break
  dt <- fread(text = lines, header = FALSE, sep = "\t", select = c(1, sel + 1))
  stopifnot(nrow(dt) == length(lines))
  i <- i + 1
  genes <- c(genes, dt[[1]])
  blocks[[i]] <- Matrix(as.matrix(dt[, -1]), sparse = TRUE)
  message("rows read: ", length(genes))
  if (length(lines) < chunk) break
}
close(con)
counts <- do.call(rbind, blocks)
dimnames(counts) <- list(genes, cell_ids[sel])
rm(blocks, lines, dt); invisible(gc())
cat("Genes read:", length(genes), "\n")

# ---------------------------------------------------------------- checks
checks <- list()
v <- counts@x
checks$integer_counts   <- all(abs(v - round(v)) < 1e-8)
checks$non_negative     <- all(v >= 0)
checks$max_value        <- max(v)
checks$duplicated_genes <- sum(duplicated(genes))
matched <- intersect(colnames(counts), annot[[COL_CELL]])
checks$cells_in_matrix        <- ncol(counts)
checks$cells_matched_to_annot <- length(matched)
checks$cells_without_annot    <- ncol(counts) - length(matched)
checks$donor_recoverable      <- all(!is.na(annot[match(matched, get(COL_CELL))][[COL_PATIENT]]))
checks$genes_read             <- length(genes)
checks$matrix_size_MB         <- round(as.numeric(object.size(counts)) / 1e6, 1)
print(checks)
writeLines(paste(names(checks), unlist(checks), sep = "\t"), file.path(res_dir, "pilot_load_checks.tsv"))

if (!checks$integer_counts) stop("Matrix is not integer counts - do NOT use for count-based methods.")

# Duplicated gene symbols are summed (documented in Methods)
if (checks$duplicated_genes > 0) {
  f <- factor(rownames(counts))
  agg <- sparseMatrix(i = as.integer(f), j = seq_along(f), x = 1,
                      dims = c(nlevels(f), length(f)), dimnames = list(levels(f), NULL))
  counts <- agg %*% counts
}

meta <- as.data.frame(annot[match(colnames(counts), get(COL_CELL))])
rownames(meta) <- meta[[COL_CELL]]
tag <- if (FULL_RUN) "full" else "pilot"
saveRDS(counts, file.path(out_dir, paste0("GSE132465_", tag, "_counts_sparse.rds")))
saveRDS(meta,   file.path(out_dir, paste0("GSE132465_", tag, "_cell_meta.rds")))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_02.txt"))
cat("Saved", tag, "objects to", out_dir, "\n")