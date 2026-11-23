
# 31_lymphoid_partition.R  (protocol v1.8, section 14.4) - additive gene partition for the lymphoid compartment
# Writes NEW files only. gene_partition.tsv and its hash are NOT modified.
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(msigdbr); library(digest) })
SEED <- 20261008; MIN_CELLS <- 20; COVER_MIN <- 0.80
pb_dir <- "04_processed_data/GSE178341"; out_dir <- "04_processed_data/programs"

progs <- data.table(program = "LYM_IL2_STAT5", compartment = "Lymphoid",
                    gene_set = "HALLMARK_IL2_STAT5_SIGNALING",
                    description = "IL2/STAT5 activation transcription")
conventional_LYM <- c("IL2RA","IL2RB","CISH","SOCS1","SOCS2","PIM1","MYC","BCL2L1")

gs <- as.data.table(readRDS("02_metadata/genesets/msigdb_lymphoid.rds"))
genes <- fread(file.path(pb_dir, "genes.csv")); genes[, ens_clean := sub("\\..*$", "", ensembl_id)]

pb <- readRDS(file.path(pb_dir, "pseudobulk_Lymphoid_CD8.rds"))
keep <- names(pb$n_cells)[pb$n_cells >= MIN_CELLS]
cp <- edgeR::cpm(DGEList(as.matrix(pb$counts[, keep])))
expressed <- rownames(cp)[rowMeans(cp >= 1) >= 0.25]
cat("eligible donors:", length(keep), "| expressed genes:", length(expressed), "\n")

ens_col <- intersect(c("ensembl_gene", "db_ensembl_gene"), names(gs))[1]
gs[, ens_clean := if (is.na(ens_col)) NA_character_ else sub("\\..*$", "", get(ens_col))]
n_ens <- sum(unique(gs$ens_clean) %in% genes$ens_clean); n_sym <- sum(unique(gs$gene_symbol) %in% genes$symbol)
MATCH_BY <- if (n_ens >= n_sym) "ensembl" else "symbol"
cat("matching by", MATCH_BY, "\n")

g <- unique(gs[gs_name == progs$gene_set, .(ens_clean, gene_symbol)])
g <- if (MATCH_BY == "ensembl") { merge(g, genes[, .(ens_clean, row_id = ensembl_id)], by = "ens_clean")
     } else { merge(g, unique(genes[, .(gene_symbol = symbol, row_id = ensembl_id)]), by = "gene_symbol") }
g <- unique(g, by = "row_id"); g[, expressed := row_id %in% expressed]
coverage <- mean(g$expressed)
if (coverage < COVER_MIN) stop("coverage below threshold - program should be dropped")

set.seed(SEED)
ge <- g[expressed == TRUE]; idx <- sample(nrow(ge))
ge[, role := ifelse(seq_len(.N) %in% idx[seq_len(ceiling(.N / 2))], "evaluation", "predictor")]
part_new <- cbind(program = progs$program, compartment = progs$compartment, ge)

old <- fread(file.path(out_dir, "gene_partition.tsv"))
eval_all <- unique(c(old[role == "evaluation", row_id], part_new[role == "evaluation", row_id]))
n_conflict <- part_new[role == "predictor" & row_id %in% eval_all, .N]
part_new <- part_new[!(role == "predictor" & row_id %in% eval_all)]
stopifnot(!any(part_new[role == "predictor", row_id] %in% eval_all))

n_eval <- part_new[role == "evaluation", .N]; n_pred <- part_new[role == "predictor", .N]
retained <- coverage >= COVER_MIN & n_eval >= 10 & n_pred >= 5
cat("\ncoverage", round(coverage, 3), "| evaluation", n_eval, "| predictor", n_pred,
    "| removed as evaluation genes elsewhere", n_conflict, "| RETAINED:", retained, "\n")
if (!retained) stop("lymphoid program fails retention - STOP")

frz <- readRDS("04_processed_data/frozen_panel/frozen_panel_v2.rds")
fg <- unlist(lapply(frz$programs, function(x) x$panel$row_id))
cat("frozen panel genes also in lymphoid evaluation pool:",
    sum(fg %in% part_new[role == "evaluation", row_id]), "of", length(fg),
    "(no leakage: different compartment, different target)\n")

defs <- cbind(progs, set_size = nrow(g), expressed = sum(g$expressed), coverage = round(coverage, 3),
              n_evaluation = n_eval, n_predictor = n_pred, retained = retained)
fwrite(defs, file.path(out_dir, "program_definitions_lymphoid.csv"))
fwrite(part_new, file.path(out_dir, "gene_partition_lymphoid.tsv"), sep = "\t")
fwrite(data.table(program = "LYM_IL2_STAT5", gene_symbol = conventional_LYM),
       file.path(out_dir, "conventional_markers_lymphoid.csv"))
h <- digest(list(defs$program, part_new[order(row_id), .(program, row_id, gene_symbol, role)]))
cat("\nLymphoid partition hash:", h, "\n")
writeLines(c(format(Sys.time()), h), file.path(out_dir, "gene_partition_lymphoid.sha"))
stopifnot(isTRUE(all.equal(fread(file.path(out_dir, "gene_partition.tsv")), old, check.attributes = FALSE)))
cat("original gene_partition.tsv verified unchanged\n")
print(defs)
 

dir.create("06_results/sessionInfo", recursive = TRUE, showWarnings = FALSE); writeLines(capture.output(sessionInfo()), "06_results/sessionInfo/sessionInfo_31.txt")
