
# 28_CD8_feasibility.R - phase 1: can a lymphoid compartment be added?
# Reads cell labels and gene-set membership only. No targets, no models.
suppressPackageStartupMessages({ library(data.table); library(msigdbr) })
MIN_CELLS <- 20

## ---- 1. what lymphoid labels exist in the discovery cohort -------------
in_dir <- "03_original_data/GSE178341"
meta <- fread(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"))
clus <- fread(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"))
d <- merge(meta, clus, by.x = "cellID", by.y = "sampleID")
Tm <- d$SPECIMEN_TYPE == "T"
cat("\n=== clMidwayPr labels (tumour specimens) ===\n")
print(d[Tm, .N, by = clMidwayPr][order(-N)])
cat("\n=== sub-labels within T/NK/ILC ===\n")
print(d[Tm & grepl("^T|NK|ILC", clMidwayPr), .N, by = .(clMidwayPr, cl295v11SubShort)][order(-N)])

## ---- 2. donor coverage, using a provisional CD8 definition -------------
# adjust the pattern below once the label tables above are read
d[, is_cd8 := Tm & grepl("CD8|cytotox", paste(clMidwayPr, cl295v11SubShort), ignore.case = TRUE)]
cd8 <- d[is_cd8 == TRUE, .N, by = PID]
cat("\n=== discovery CD8 cells per donor ===\n")
cat("donors with any CD8:", nrow(cd8), "| donors >=", MIN_CELLS, "cells:", sum(cd8$N >= MIN_CELLS),
    "| median cells:", median(cd8$N), "\n")
print(summary(cd8$N))

## ---- 3. candidate gene sets and predictor-pool survival ----------------
get_set <- function(nm) {
  g <- tryCatch(as.data.table(msigdbr(species = "Homo sapiens"))[gs_name == nm],
                error = function(e) NULL)
  if (is.null(g) || !nrow(g)) return(character(0))
  unique(g[[grep("^(gene_symbol|human_gene_symbol)$", names(g), value = TRUE)[1]]])
}
cand <- c("GOBP_T_CELL_MEDIATED_CYTOTOXICITY",
          "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
          "HALLMARK_IL2_STAT5_SIGNALING",
          "HALLMARK_ALLOGRAFT_REJECTION")

# existing evaluation genes = genes a new panel may NOT use as predictors
pf <- list.files("04_processed_data/programs", full.names = TRUE)
cat("\nfiles in 04_processed_data/programs:\n"); print(basename(pf))
part <- readRDS(grep("partition|program", pf, value = TRUE)[1])
cat("\nstructure of the partition object:\n"); print(names(part))
evalg <- tryCatch(unique(unlist(lapply(part, function(x) x$evaluation_genes))),
                  error = function(e) character(0))
cat("existing evaluation genes:", length(evalg), "\n")

cat("\n=== candidate CD8 gene sets ===\n")
for (nm in cand) {
  g <- get_set(nm)
  if (!length(g)) { cat(sprintf("%-70s NOT FOUND in msigdbr\n", nm)); next }
  blocked <- sum(g %in% evalg)
  cat(sprintf("%-70s size %4d | blocked by existing eval genes %3d (%.0f%%) | free %4d\n",
              nm, length(g), blocked, 100 * blocked / length(g), length(g) - blocked))
}
 
