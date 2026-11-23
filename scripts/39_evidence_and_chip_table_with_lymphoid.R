
# 39_evidence_and_chip_table_with_lymphoid.R  (protocol v1.9) - 8-program evidence + chip tables.
# Evidence logic copied from 20_evidence_and_chip_table.R. Existing 7 rows are read, not recomputed.
suppressPackageStartupMessages({ library(data.table); library(digest) })
if (digest(file = "04_processed_data/frozen_panel/frozen_panel_v2.rds", algo = "sha256") !=
    "a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb") stop("frozen_panel_v2 altered - STOP")
if (digest(file = "04_processed_data/frozen_panel/frozen_panel_LYM_v1.rds", algo = "sha256") !=
    "7dc3256d97477a796fcacd1df71021740888f7978a01ea36caf70787da7f513f") stop("lymphoid panel altered - STOP")
frzL <- readRDS("04_processed_data/frozen_panel/frozen_panel_LYM_v1.rds")

E7  <- fread("08_tables/Table_program_evidence.csv")
CH7 <- fread("08_tables/Table_chip_assay_proposal.csv")
M   <- fread("06_results/validation_lymphoid/meta_analysis_results_lymphoid.csv")
SP  <- fread("06_results/spatial_support_lymphoid_combined/spatial_support_summary_combined.csv")[analysis == "All spatial cohorts (main)"]
AP  <- fread("06_results/spatial_support_lymphoid_combined/criterionA_patients_combined.csv")

dl <- function(est, se) { w <- 1 / se^2; fe <- sum(w * est) / sum(w); k <- length(est); Q <- sum(w * (est - fe)^2)
  t2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))); wr <- 1 / (se^2 + t2)
  m <- sum(wr * est) / sum(wr); s <- sqrt(1 / sum(wr)); c(est = m, lo = m - 1.96 * s, hi = m + 1.96 * s) }
lo_ <- AP[is.finite(se) & se > 0, { r <- t(sapply(unique(patient), function(p) dl(effect[patient != p], se[patient != p])))
  .(lopo_min_effect = min(r[, "est"]), lopo_min_lo = min(r[, "lo"]), A_robust_all_LOPO = all(r[, "lo"] > 0)) }, by = program]

mf <- M[scaling == "frozen"]; mw <- M[scaling == "within"]
r  <- mf$pooled_MAE_ratio
ci <- sprintf("%.2f (%.2f-%.2f)", r, mf$ratio_lo, mf$ratio_hi)
src <- paste0("meta-analysis, ", mf$k_cohorts, " cohorts, ", mf$total_donors, " donors")
indep <- mf$meta_status == "Supported"
wsup  <- isTRUE(mw$meta_status == "Supported")
explor <- !indep && (r < 1 || wsup)
cat_ <- if (indep && SP$spatial_label == "Spatial support") "Independently supported with spatial evidence" else
        if (indep) "Independently supported in single-cell data" else if (explor) "Exploratory" else "Unsupported"
P <- frzL$programs[["LYM_IL2_STAT5"]]
Enew <- data.table(program = "LYM_IL2_STAT5", description = "IL2/STAT5 activation transcription (CD8 T cells)",
  compartment = "Lymphoid", panel = paste(P$panel$gene_symbol, collapse = ", "),
  sc_source = src, frozen_MAE_ratio_CI = ci, sc_independent = indep, within_scaling_supported = wsup,
  spatial_label = SP$spatial_label,
  spatial_effect_CI = sprintf("%.3f (%.3f-%.3f)", SP$pooled_effect, SP$lo, SP$hi), evidence_category = cat_)
Enew <- merge(Enew, lo_, by = "program", all.x = TRUE)
E <- rbind(E7, Enew, fill = TRUE)

chipL <- "None in the configuration considered (the chip has no lymphoid compartment)"
limL <- paste("Proposed for a chip configuration that does not yet exist: THP-1-derived macrophages are myeloid, not lymphoid,",
  "so no current compartment expresses this programme. Measuring it would require adding primary CD8 T cells.",
  "Spatial composition was adjusted using a pooled T/NK reference type, which is broader than the CD8 population the panel was derived from.",
  "RNA model coefficients are not transferable to qPCR delta-Ct or protein readouts without separate calibration.")
CHnew <- data.table(gene = P$panel$gene_symbol,
  coefficient_sign = ifelse(P$panel$coefficient > 0, "+", "-"),
  compartment = "Lymphoid", chip_cell_model = chipL,
  program = "IL2/STAT5 activation transcription (CD8 T cells)",
  single_cell_evidence = paste0(ci, "; ", src), spatial_evidence = SP$spatial_label, evidence_category = cat_,
  proposed_future_assay = "RT-qPCR in a CD8 T-cell compartment, if added to the chip, relative to within-experiment reference conditions",
  interpretation_limits = limL)
CH <- rbind(CH7, CHnew, fill = TRUE)

fwrite(E, "08_tables/Table_program_evidence_8.csv"); fwrite(CH, "08_tables/Table_chip_assay_proposal_8.csv")
cat("\n=== PROGRAM EVIDENCE (8 programs) ===\n")
print(E[, .(program, frozen_MAE_ratio_CI, sc_independent, within_scaling_supported, spatial_label, A_robust_all_LOPO, evidence_category)])
cat("\n=== CHIP TABLE:", nrow(CH), "genes ===\n")
print(CH[, .N, by = .(compartment, evidence_category)])
writeLines(capture.output(sessionInfo()), "06_results/validation_lymphoid/sessionInfo_39.txt")
 
