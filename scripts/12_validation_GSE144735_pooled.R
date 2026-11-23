# 12_validation_GSE144735_pooled.R
# Protocol v0.9: frozen panel v2 applied to GSE144735 (KUL3, Belgium; tumour core only).
# KUL has < 10 donors -> reported descriptively; pooled SMC + KUL = pre-specified sensitivity analysis.
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(digest) })
set.seed(20261008)
MIN_CELLS <- 20; NBOOT <- 2000; MIN_DONORS <- 10; MIN_EVAL_COVER <- 0.80
frz_file <- "04_processed_data/frozen_panel/frozen_panel_v2.rds"
EXPECTED_SHA <- "a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb"
if (digest(file = frz_file, algo = "sha256") != EXPECTED_SHA) stop("Frozen panel hash mismatch - STOP")
frz <- readRDS(frz_file); cat("Frozen panel verified\n")
in_dir <- "03_original_data/GSE144735"; res_dir <- "06_results/validation_GSE144735"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("04_processed_data/GSE144735", recursive = TRUE, showWarnings = FALSE)
annot <- fread(file.path(in_dir, "GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz"))
cat("Annotation columns:", paste(names(annot), collapse = ", "), "\n"); print(annot[, .N, by = .(Class)])
print(annot[, .N, by = .(Cell_type, Cell_subtype)][order(Cell_type, -N)], nrows = 100)
con <- gzfile(file.path(in_dir, "GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz"), "r")
header <- gsub("\"", "", strsplit(readLines(con, 1), "\t", fixed = TRUE)[[1]])
first <- readLines(con, 1); nf <- length(strsplit(first, "\t", fixed = TRUE)[[1]])
cells <- if (length(header) == nf) header[-1] else header
blocks <- list(); genes <- character(0); pending <- first
repeat {
  lines <- c(pending, readLines(con, n = 2000 - length(pending))); pending <- character(0)
  if (!length(lines)) break
  dt <- fread(text = lines, header = FALSE, sep = "\t")
  genes <- c(genes, dt[[1]]); blocks[[length(blocks) + 1]] <- as(as.matrix(dt[, -1]), "CsparseMatrix"); message("KUL rows read: ", length(genes))
  if (length(lines) < 2000) break
}
close(con)
counts <- do.call(rbind, blocks); dimnames(counts) <- list(genes, cells); rm(blocks)
stopifnot(all(abs(counts@x - round(counts@x)) < 1e-8))
saveRDS(counts, "04_processed_data/GSE144735/GSE144735_counts_sparse.rds")
meta <- annot[match(cells, Index)]
cat("Cells:", ncol(counts), "| genes:", nrow(counts), "| matched to annotation:", sum(!is.na(meta$Patient)), "\n")
mt <- grepl("^MT-", genes)
meta[, nCount := Matrix::colSums(counts)]; meta[, nFeature := Matrix::colSums(counts > 0)]
meta[, pct_mt := 100 * Matrix::colSums(counts[mt, , drop = FALSE]) / nCount]
mad_flag <- function(x, n = 3, side = c("both", "high")) { side <- match.arg(side); m <- median(x); s <- mad(x)
  if (!is.finite(s) || s == 0) return(rep(FALSE, length(x)))
  if (side == "high") x > m + n * s else x < m - n * s | x > m + n * s }
meta[, qc_flag := mad_flag(log10(nCount + 1)) | mad_flag(log10(nFeature + 1)) | mad_flag(pct_mt, side = "high"), by = Sample]
meta[, val_compartment := fcase(
  Cell_subtype %in% c("Tip-like ECs", "Stalk-like ECs", "Proliferative ECs"), "Endothelial",
  Cell_subtype %in% c("Stromal 1", "Stromal 2", "Stromal 3", "Myofibroblasts"), "Stromal",
  Cell_type == "Myeloids" & !Cell_subtype %in% c("cDC", "Unknown"), "Myeloid",
  Cell_type == "Epithelial cells" & grepl("CMS", Cell_subtype), "Epithelial",
  default = NA_character_)]
use <- meta$Class == "Tumor" & !meta$qc_flag & !is.na(meta$val_compartment)
 
