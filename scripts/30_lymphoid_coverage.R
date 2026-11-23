
# 30_lymphoid_coverage.R  (protocol v1.8, section 14.4) - coverage gate for the lymphoid programs
# Uses exactly the 'expressed' rule of script 09: CPM >= 1 in >= 25% of eligible donors.
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(msigdbr) })
MIN_CELLS <- 20; COVER_MIN <- 0.80
pb_dir <- "04_processed_data/GSE178341"

progs <- data.table(
  program  = c("LYM_CYTOTOX", "LYM_IL2_STAT5"),
  gene_set = c("GOBP_T_CELL_MEDIATED_CYTOTOXICITY", "HALLMARK_IL2_STAT5_SIGNALING"),
  description = c("Cytotoxic T-cell effector transcription", "IL2/STAT5 activation transcription"))

cache <- "02_metadata/genesets/msigdb_lymphoid.rds"
if (file.exists(cache)) gs <- readRDS(cache) else {
  h  <- msigdbr(species = "Homo sapiens", collection = "H")
  go <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
  gs <- as.data.table(rbind(h, go))[gs_name %in% progs$gene_set]
  saveRDS(gs, cache) }
gs <- as.data.table(gs)
stopifnot(all(progs$gene_set %in% gs$gs_name))
cat("msigdbr", as.character(packageVersion("msigdbr")),
    "| db_version", if ("db_version" %in% names(gs)) unique(gs$db_version)[1] else NA, "\n")

genes <- fread(file.path(pb_dir, "genes.csv"))
genes[, ens_clean := sub("\\..*$", "", ensembl_id)]
expressed_in <- function(comp) {
  pb <- readRDS(file.path(pb_dir, paste0("pseudobulk_", comp, ".rds")))
  keep <- names(pb$n_cells)[pb$n_cells >= MIN_CELLS]
  cp <- edgeR::cpm(DGEList(as.matrix(pb$counts[, keep])))
  list(genes = rownames(cp)[rowMeans(cp >= 1) >= 0.25], n_donors = length(keep)) }
e <- expressed_in("Lymphoid_CD8")
cat("eligible donors:", e$n_donors, "| expressed genes in compartment:", length(e$genes), "\n\n")

ens_col <- intersect(c("ensembl_gene", "db_ensembl_gene"), names(gs))[1]
gs[, ens_clean := if (is.na(ens_col)) NA_character_ else sub("\\..*$", "", get(ens_col))]
n_ens <- sum(unique(gs$ens_clean) %in% genes$ens_clean); n_sym <- sum(unique(gs$gene_symbol) %in% genes$symbol)
MATCH_BY <- if (n_ens >= n_sym) "ensembl" else "symbol"
cat("matching by", MATCH_BY, "( ensembl", n_ens, "| symbol", n_sym, ")\n\n")

out <- rbindlist(lapply(seq_len(nrow(progs)), function(k) {
  p <- progs[k]; g <- unique(gs[gs_name == p$gene_set, .(ens_clean, gene_symbol)])
  g <- if (MATCH_BY == "ensembl") merge(g, genes[, .(ens_clean, row_id = ensembl_id)], by = "ens_clean")
       else merge(g, unique(genes[, .(gene_symbol = symbol, row_id = ensembl_id)]), by = "gene_symbol")
  g <- unique(g, by = "row_id"); g[, expressed := row_id %in% e$genes]
  data.table(program = p$program, gene_set = p$gene_set, set_size = nrow(g),
             expressed = sum(g$expressed), coverage = round(mean(g$expressed), 3),
             retained = mean(g$expressed) >= COVER_MIN) }))
cat("=== LYMPHOID COVERAGE GATE (threshold", COVER_MIN, ") ===\n"); print(out)
fwrite(out, "06_results/programs/lymphoid_coverage.csv")
 

dir.create("06_results/sessionInfo", recursive = TRUE, showWarnings = FALSE); writeLines(capture.output(sessionInfo()), "06_results/sessionInfo/sessionInfo_30.txt")
