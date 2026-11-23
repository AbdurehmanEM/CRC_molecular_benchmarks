# 09_programs_and_gene_partition.R
# Protocol v0.6, steps 9-10: fix program gene sets, check coverage in the DISCOVERY
# pseudobulk only, and split each program's genes into non-overlapping
#   EVALUATION genes (define the target score) and
#   PREDICTOR genes (the only genes the compact panel may use).
# No validation-cohort data are read here.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(msigdbr); library(digest) })
SEED <- 20261008
MIN_CELLS <- 20
COVER_MIN <- 0.80          # program kept only if >= 80% of its genes are expressed
pb_dir  <- "04_processed_data/GSE178341"
out_dir <- "04_processed_data/programs"
res_dir <- "06_results/programs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("02_metadata/genesets", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------- 1. program definitions (pre-specified)
programs <- data.table(
  program     = c("EPI_HYPOXIA", "EPI_APOPTOSIS", "STR_TGFB", "STR_ECM",
                  "MYE_INFLAM", "MYE_IFNG", "END_ANGIO", "END_TNFA_NFKB"),
  compartment = c("Epithelial", "Epithelial", "Stromal", "Stromal",
                  "Myeloid", "Myeloid", "Endothelial", "Endothelial"),
  gene_set    = c("HALLMARK_HYPOXIA", "HALLMARK_APOPTOSIS", "HALLMARK_TGF_BETA_SIGNALING",
                  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
                  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                  "HALLMARK_ANGIOGENESIS", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  description = c("Hypoxia-associated transcription", "Apoptosis-related transcription",
                  "TGF-beta response transcription", "Extracellular-matrix organisation transcription",
                  "Inflammatory-response transcription", "Interferon-gamma response transcription",
                  "Angiogenesis-associated transcription", "TNF/NF-kB (inflammatory activation) transcription"))

# MSigDB via msigdbr (needs internet the first time; a copy is cached for offline reuse)
cache <- "02_metadata/genesets/msigdb_programs.rds"
if (file.exists(cache)) gs <- readRDS(cache) else {
  h  <- msigdbr(species = "Homo sapiens", collection = "H")
  rc <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
  gs <- as.data.table(rbind(h, rc))[gs_name %in% programs$gene_set]
  saveRDS(gs, cache)
}
gs <- as.data.table(gs)
stopifnot(all(programs$gene_set %in% gs$gs_name))
version_info <- data.table(msigdbr = as.character(packageVersion("msigdbr")),
                           db_version = if ("db_version" %in% names(gs)) unique(gs$db_version)[1] else NA,
                           retrieved = as.character(Sys.Date()))
print(version_info)
cat("Gene-set rows:", nrow(gs), "| columns:", paste(names(gs), collapse = ", "), "\n")

# ---------------------------------------------------------------- 2. coverage in discovery pseudobulk
genes <- fread(file.path(pb_dir, "genes.csv"))
expressed_in <- function(comp) {
  pb <- readRDS(file.path(pb_dir, paste0("pseudobulk_", comp, ".rds")))
  keep <- names(pb$n_cells)[pb$n_cells >= MIN_CELLS]
  y <- DGEList(as.matrix(pb$counts[, keep]))
  cp <- cpm(y)
  # expressed = CPM >= 1 in >= 25% of eligible donors
  rownames(cp)[rowMeans(cp >= 1) >= 0.25]
}
expr <- lapply(setNames(unique(programs$compartment), unique(programs$compartment)), expressed_in)

# ---- gene-ID matching (strip version suffixes; use whichever ID type matches better)
genes[, ens_clean := sub("\\..*$", "", ensembl_id)]
ens_col <- intersect(c("ensembl_gene", "db_ensembl_gene"), names(gs))[1]
gs[,    ens_clean := if (is.na(ens_col)) NA_character_ else sub("\\..*$", "", get(ens_col))]
n_ens <- sum(unique(gs$ens_clean)   %in% genes$ens_clean)
n_sym <- sum(unique(gs$gene_symbol) %in% genes$symbol)
MATCH_BY <- if (n_ens >= n_sym) "ensembl" else "symbol"
cat("\nGene-set genes found in data: by Ensembl =", n_ens, "| by symbol =", n_sym,
    "-> matching by", MATCH_BY, "\n")
print(head(genes, 3))
stopifnot(max(n_ens, n_sym) > 1000)            # stop if matching still fails

part <- list(); cov <- list()
set.seed(SEED)
for (k in seq_len(nrow(programs))) {
  p <- programs[k]
  g <- unique(gs[gs_name == p$gene_set, .(ens_clean, gene_symbol)])
  g <- if (MATCH_BY == "ensembl") {
    merge(g, genes[, .(ens_clean, row_id = ensembl_id)], by = "ens_clean")
  } else {
    merge(g, unique(genes[, .(gene_symbol = symbol, row_id = ensembl_id)]), by = "gene_symbol")
  }
  g <- unique(g, by = "row_id")
  g[, expressed := row_id %in% expr[[p$compartment]]]
  cov[[k]] <- data.table(program = p$program, set_size = nrow(g),
                         expressed = sum(g$expressed), coverage = round(mean(g$expressed), 3))
  ge <- g[expressed == TRUE]
  idx <- sample(nrow(ge))
  ge[, role := ifelse(seq_len(.N) %in% idx[seq_len(ceiling(.N / 2))], "evaluation", "predictor")]
  part[[k]] <- cbind(program = p$program, compartment = p$compartment, ge)
}
cov <- rbindlist(cov); part <- rbindlist(part)

eval_all <- unique(part[role == "evaluation", row_id])
n_conflict <- part[role == "predictor" & row_id %in% eval_all, .N]
part <- part[!(role == "predictor" & row_id %in% eval_all)]
stopifnot(!any(part[role == "predictor", row_id] %in% eval_all))

cov <- merge(cov, part[, .(n_evaluation = sum(role == "evaluation"),
                           n_predictor  = sum(role == "predictor")), by = program], by = "program")
cov[, retained := coverage >= COVER_MIN & n_evaluation >= 10 & n_predictor >= 5]
programs <- merge(programs, cov, by = "program")
cat("\n=== Program coverage in discovery (GSE178341) ===\n"); print(programs[, .(program, gene_set, set_size, expressed, coverage, n_evaluation, n_predictor, retained)])
cat("\nPredictor genes removed because they are evaluation genes in another program:", n_conflict, "\n")

# ---------------------------------------------------------------- 4. save (fixed before modelling)
fwrite(programs, file.path(out_dir, "program_definitions.csv"))
fwrite(part,     file.path(out_dir, "gene_partition.tsv"), sep = "\t")
fwrite(version_info, file.path(out_dir, "geneset_version.csv"))
h <- digest(list(programs$program, part[order(program, row_id), .(program, row_id, gene_symbol, role)]))
cat("\nPartition hash:", h, "\n")
writeLines(c(format(Sys.time()), h), file.path(out_dir, "gene_partition.sha"))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_09.txt"))
