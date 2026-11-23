
# 20_evidence_and_chip_table.R  (protocol v1.3, section 11) - evidence categories, chip assay table, spatial LOPO sensitivity
suppressPackageStartupMessages({ library(data.table); library(digest) })
frz_file <- "04_processed_data/frozen_panel/frozen_panel_v2.rds"
if (digest(file = frz_file, algo = "sha256") != "a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb") stop("Frozen panel hash mismatch - STOP")
frz <- readRDS(frz_file); dir.create("08_tables", showWarnings = FALSE)
M <- fread("06_results/cross_cohort_meta/meta_analysis_results.csv")
V <- fread("06_results/validation_GSE132465/validation_results.csv"); V[scaling == "within_cohort", scaling := "within"]
SP <- fread("06_results/spatial_support_combined/spatial_support_summary_main.csv")
AP <- fread("06_results/spatial_support_combined/criterionA_patients_combined.csv")
dl <- function(est, se) { w <- 1 / se^2; fe <- sum(w * est) / sum(w); k <- length(est); Q <- sum(w * (est - fe)^2)
  t2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))); wr <- 1 / (se^2 + t2); m <- sum(wr * est) / sum(wr); s <- sqrt(1 / sum(wr))
  c(est = m, lo = m - 1.96 * s, hi = m + 1.96 * s) }
lopo <- AP[is.finite(se) & se > 0, { r <- t(sapply(unique(patient), function(p) dl(effect[patient != p], se[patient != p])))
  .(lopo_min_effect = min(r[, "est"]), lopo_min_lo = min(r[, "lo"]), A_robust_all_LOPO = all(r[, "lo"] > 0)) }, by = program]
desc <- c(EPI_HYPOXIA = "Hypoxia-associated transcription", EPI_APOPTOSIS = "Apoptosis-related transcription",
  STR_TGFB = "TGF-beta response", MYE_IFNG = "Interferon-gamma response", MYE_INFLAM = "Inflammatory response",
  END_ANGIO = "Angiogenesis-associated transcription", END_TNFA_NFKB = "TNF-alpha/NF-kB endothelial activation")
chip <- c(Epithelial = "HCT116 or SW48 (MMR-deficient); SW480 or HT-29 (MMR-proficient)", Stromal = "hMSC", Myeloid = "THP-1-derived macrophages", Endothelial = "HUVEC")
lim <- c(Epithelial = "One MMR-deficient and one MMR-proficient line are proposed because the discovery cohort contained 34 MMR-deficient and 28 MMR-proficient donors; a single microsatellite class does not represent that population. Cell lines lack patient tumour heterogeneity; transcriptional signature only, not a direct measure of oxygen tension or cell death",
  Stromal = "hMSC are not cancer-associated fibroblasts; TGF-beta response transcripts may differ in baseline and dynamic range",
  Myeloid = "THP-1 is a monocytic leukaemia line; differentiation and polarisation protocol strongly affects these transcripts",
  Endothelial = "HUVEC are venous, non-tumour endothelium; tumour endothelial states are only partly represented")
E <- rbindlist(lapply(names(frz$programs), function(pg) {
  P <- frz$programs[[pg]]; mf <- M[analysis == "main" & scaling == "frozen" & program == pg]; mw <- M[analysis == "main" & scaling == "within" & program == pg]
  if (nrow(mf) && is.finite(mf$pooled_MAE_ratio)) {
    src <- paste0("meta-analysis, ", mf$k_cohorts, " cohorts, ", mf$total_donors, " donors"); r <- mf$pooled_MAE_ratio
    ci <- sprintf("%.2f (%.2f-%.2f)", r, mf$ratio_lo, mf$ratio_hi); indep <- mf$meta_status == "Supported"; wsup <- isTRUE(mw$meta_status == "Supported")
  } else { s <- V[scaling == "frozen" & program == pg]; sw <- V[scaling == "within" & program == pg]
    src <- paste0("single evaluable cohort GSE132465, ", s$n_donors, " donors"); r <- s$MAE_ratio
    ci <- sprintf("%.2f (%.2f-%.2f)", r, s$ratio_lo, s$ratio_hi); indep <- s$status == "Supported"; wsup <- isTRUE(sw$status == "Supported") }
  sp <- SP[program == pg]; explor <- !indep && (r < 1 || wsup)
  cat_ <- if (indep && sp$spatial_label == "Spatial support") "Independently supported with spatial evidence" else
          if (indep) "Independently supported in single-cell data" else if (explor) "Exploratory" else "Unsupported"
  data.table(program = pg, description = desc[[pg]], compartment = P$compartment, panel = paste(P$panel$gene_symbol, collapse = ", "),
    sc_source = src, frozen_MAE_ratio_CI = ci, sc_independent = indep, within_scaling_supported = wsup,
    spatial_label = sp$spatial_label, spatial_effect_CI = sprintf("%.3f (%.3f-%.3f)", sp$pooled_effect, sp$lo, sp$hi),
    evidence_category = cat_)
}))
E <- merge(E, lopo, by = "program", all.x = TRUE)
CH <- rbindlist(lapply(names(frz$programs), function(pg) { P <- frz$programs[[pg]]; e <- E[program == pg]
  data.table(gene = P$panel$gene_symbol, coefficient_sign = ifelse(P$panel$coefficient > 0, "+", "-"),
    compartment = P$compartment, chip_cell_model = chip[[P$compartment]], program = desc[[pg]],
    single_cell_evidence = paste0(e$frozen_MAE_ratio_CI, "; ", e$sc_source), spatial_evidence = e$spatial_label,
    evidence_category = e$evidence_category,
    proposed_future_assay = "RT-qPCR in the corresponding chip compartment, relative to within-experiment reference conditions",
    interpretation_limits = paste0(lim[[P$compartment]], ". RNA model coefficients are not transferable to qPCR delta-Ct or protein readouts without separate calibration.")) }))
fwrite(E, "08_tables/Table_program_evidence.csv"); fwrite(CH, "08_tables/Table_chip_assay_proposal.csv")
cat("\n=== PROGRAM EVIDENCE ===\n"); print(E[, .(program, frozen_MAE_ratio_CI, sc_independent, within_scaling_supported, spatial_label, A_robust_all_LOPO, evidence_category)])
cat("\n=== CHIP ASSAY TABLE (", nrow(CH), "genes) ===\n"); print(CH[, .(gene, coefficient_sign, chip_cell_model, program, evidence_category)], nrows = 50)
 

dir.create("06_results/sessionInfo", recursive = TRUE, showWarnings = FALSE); writeLines(capture.output(sessionInfo()), "06_results/sessionInfo/sessionInfo_20.txt")
