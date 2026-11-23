# 06_audit_GSE178341.R
# Audit the DISCOVERY cohort (Pelka 2021) without loading the full matrix:
# metadata columns, patients, tumour/normal, MMR status, sorting/enrichment,
# h5 structure, raw-count check on a slice, cell-ID match.

if (!requireNamespace("hdf5r", quietly = TRUE)) install.packages("hdf5r", type = "binary")
suppressPackageStartupMessages({ library(data.table); library(hdf5r) })

in_dir  <- "03_original_data/GSE178341"
res_dir <- "06_results/GSE178341_audit"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(res_dir, "audit_log.txt"), split = TRUE)

meta <- fread(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"))
clus <- fread(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"))

cat("\n=== METATABLE:", nrow(meta), "rows x", ncol(meta), "columns ===\n")
print(names(meta)); print(head(meta, 3))
cat("\n=== CLUSTER TABLE:", nrow(clus), "rows x", ncol(clus), "columns ===\n")
print(names(clus)); print(head(clus, 3))

cat("\n=== Categorical columns (<= 100 distinct values) ===\n")
for (cn in names(meta)) {
  u <- uniqueN(meta[[cn]])
  if (u <= 100) { cat("\n--", cn, "(", u, "values )\n"); print(meta[, .N, by = cn][order(-N)], nrows = 100) }
}

sort_cols <- grep("sort|cd45|frac|enrich|facs|processing|platform|chem|batch|kit",
                  names(meta), ignore.case = TRUE, value = TRUE)
cat("\n=== Possible sorting/batch columns:", paste(sort_cols, collapse = ", "), "===\n")

cat("\n=== Cell-type columns in cluster table ===\n")
for (cn in names(clus)) if (uniqueN(clus[[cn]]) <= 120) {
  cat("\n--", cn, "\n"); print(clus[, .N, by = cn][order(-N)], nrows = 120)
}

h5 <- H5File$new(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit.h5"), mode = "r")
cat("\n=== h5 contents ===\n"); print(h5$ls(recursive = TRUE)[, c("name", "obj_type", "dataset.dims")])
grp   <- names(h5)[1]
shape <- h5[[paste0(grp, "/shape")]]$read()
bc    <- h5[[paste0(grp, "/barcodes")]]$read()
genes <- tryCatch(h5[[paste0(grp, "/features/name")]]$read(),
                  error = function(e) h5[[paste0(grp, "/gene_names")]]$read())
x_slice <- h5[[paste0(grp, "/data")]]$read(args = list(1:1e6))
h5$close_all()

id_col <- names(meta)[vapply(meta, function(v) is.character(v) && mean(bc[1:1000] %in% v) > 0.5, logical(1))][1]
checks <- list(
  genes_x_cells        = paste(shape, collapse = " x "),
  n_barcodes           = length(bc),
  n_gene_names         = length(genes),
  duplicated_genes     = sum(duplicated(genes)),
  integer_first_1e6    = all(abs(x_slice - round(x_slice)) < 1e-8),
  min_first_1e6        = min(x_slice),
  max_first_1e6        = max(x_slice),
  meta_cell_id_column  = id_col,
  barcodes_in_meta     = if (is.na(id_col)) NA else sum(bc %in% meta[[id_col]]),
  meta_rows_in_h5      = if (is.na(id_col)) NA else sum(meta[[id_col]] %in% bc)
)
cat("\n=== CHECKS ===\n"); print(checks)
writeLines(paste(names(checks), unlist(checks), sep = "\t"), file.path(res_dir, "h5_checks.tsv"))
sink()
cat("\nAudit saved to", res_dir, "\n")