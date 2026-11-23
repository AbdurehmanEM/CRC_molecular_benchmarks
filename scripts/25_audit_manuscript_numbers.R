# 25_audit_manuscript_numbers.R - traces every numerical claim in the manuscript to its source file
suppressPackageStartupMessages({library(data.table); library(digest)})
# 19 Oct 2026: hard-coded setwd("F:/...") removed - it stopped 00_run_all.R (stage 26)
# on any other machine, which blocked every later stage. Run from the project root.
stopifnot(dir.exists("05_analysis_scripts"))
sep <- function(x) cat("\n========== ", x, " ==========\n")
has <- function(f) if (file.exists(f)) TRUE else {cat("MISSING:", f, "\n"); FALSE}
rd  <- function(f) if (has(f)) fread(f) else NULL

sep("0. FILE INDEX")
print(data.table(file = list.files(c("06_results","08_tables","04_processed_data/programs"),
  pattern="\\.(csv|tsv)$", recursive=TRUE, full.names=TRUE)))

sep("1. FROZEN PANEL HASH + PANEL GENES")
fz <- "04_processed_data/frozen_panel/frozen_panel_v2.rds"
cat("sha256:", digest(file=fz, algo="sha256"), "\n")
cat("expected: a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb\n")
P <- readRDS(fz)$programs
print(rbindlist(lapply(names(P), function(k) data.table(program=k, compartment=P[[k]]$compartment,
  n_panel=nrow(P[[k]]$panel), genes=paste(P[[k]]$panel$gene_symbol, collapse=","),
  n_eval=nrow(P[[k]]$evaluation_genes)))))

sep("2. PROGRAM DEFINITIONS (coverage, ECM 77.4%?) + PARTITION (129 removed?)")
pd <- rd("04_processed_data/programs/program_definitions.csv"); if(!is.null(pd)) print(pd)
gp <- "04_processed_data/programs/gene_partition.tsv"
if (has(gp)) { g <- fread(gp); print(g[, .N, by=.(program, role)][order(program, role)]) }

sep("3. DISCOVERY DONOR COUNTS (expect 60/27/60/40 at >=20 cells)")
for (cmp in c("Epithelial","Stromal","Myeloid","Endothelial")) {
  f <- sprintf("04_processed_data/GSE178341/pseudobulk_%s.rds", cmp)
  if (file.exists(f)) { pb <- readRDS(f); cat(sprintf("%-12s n>=20: %d | n>=10: %d | n>=30: %d | n>=50: %d\n",
    cmp, sum(pb$n_cells>=20), sum(pb$n_cells>=10), sum(pb$n_cells>=30), sum(pb$n_cells>=50))) }
  else cat("MISSING:", f, "\n")
}

sep("4. DISCOVERY CV PERFORMANCE (MAE ratios, calibration slopes)")
cv <- rd("06_results/panel/cv_performance_discovery_v2.csv")
if (!is.null(cv)) { print(dcast(cv, program ~ method, value.var="MAE_vs_mean"))
  print(cv[method=="panel_refit", .(program, cal_slope, cal_intercept, spearman)]) }   # column is cal_slope (was looked up as calib_slope and silently skipped)
cv1 <- rd("06_results/panel/cv_performance_discovery.csv")
if (!is.null(cv1)) { cat("\n-- v1 (shrunk) slopes 1.45-5.61? --\n")
  print(cv1[method %in% c("panel","panel_refit"), .(program, method,
    slope=if("calib_slope" %in% names(cv1)) calib_slope else NA, MAE_vs_mean)]) }

sep("5. SELECTION STABILITY (15/26 genes >=75%? STAT1/IFIT3/NAMPT/MSX1/KDM6B =100%? ARID4B 26%? FOSB?)")
# explicit file: pattern "stab" + [1] picked selection_stability.csv (v1) or
# selection_stability_lymphoid.csv, whichever sorts first - never the v2 table
st <- "06_results/panel/selection_stability_v2.csv"
if (has(st)) { s <- fread(st); print(s[order(-selection_freq)][1:40]) }

sep("6. VALIDATION PER COHORT (frozen scaling)")
vf <- c(GSE132465="06_results/validation_GSE132465/validation_results.csv",
        GSE166555="06_results/validation_GSE166555/validation_results.csv",
        GSE200997="06_results/validation_GSE200997/validation_results.csv")
for (k in names(vf)) if (has(vf[k])) { cat("\n--", k, "--\n")
  x <- fread(vf[k]); print(x[scaling %in% c("frozen"), .(program, n_donors, MAE_ratio, ratio_lo, ratio_hi, spearman, status)]) }
kf <- "06_results/validation_GSE144735/validation_results_KUL_and_pooled.csv"
if (has(kf)) { cat("\n-- GSE144735 / pooled SMC+KUL (expect TNF 0.46 0.40-0.52; angio 2.22 1.49-3.35) --\n")
  x <- fread(kf); print(x[scaling=="frozen", .(analysis, program, n_donors, MAE_ratio, ratio_lo, ratio_hi, status)]) }

sep("7. META-ANALYSIS (frozen + within) & OFFSET DIAGNOSTIC")
m <- rd("06_results/cross_cohort_meta/meta_analysis_results.csv")
if (!is.null(m)) print(m[analysis=="main", .(scaling, program, k_cohorts, pooled_MAE_ratio, ratio_lo, ratio_hi,
  I2_ratio, pooled_rho, meta_status)])
o <- rd("06_results/cross_cohort_meta/offset_diagnostic_targets.csv")
if (!is.null(o)) { cat("\n-- offset: mean predicted - observed (expect -0.2..-0.7; IFN-g GSE166555 -1.08) --\n")
  print(dcast(o, program ~ cohort, value.var="mean_offset_pred_minus_obs")) }

sep("8. SPATIAL - primary cohort (6 patients) spot total (expect 15,904)")
sp <- "04_processed_data/Valdeolivas_Visium/spot_program_scores.rds"
if (file.exists(sp)) { S <- readRDS(sp); cat("unique spots:", uniqueN(S[, paste(section, barcode)]),
  "| sections:", uniqueN(S$section), "| patients:", uniqueN(S$patient), "\n") } else cat("MISSING:", sp, "\n")

sep("9. SPATIAL - coverage + marker check in added cohorts")
cv2 <- rd("06_results/spatial_support_extra/spatial_gene_coverage_extra.csv")
if (!is.null(cv2)) print(cv2[panel_present < panel_genes | eval_coverage < 0.8])
mk <- rd("06_results/spatial_support_extra/RCTD_marker_check_extra.csv")
if (!is.null(mk)) { cat("\n-- own-marker rho range per section --\n")
  print(mk[, .(min_rho=min(own_rho), max_rho=max(own_rho), endo=own_rho[type=="Endothelial"]), by=section]) }

sep("10. SPATIAL - combined 13-patient summary (all three analyses)")
sc <- rd("06_results/spatial_support_combined/spatial_support_summary_combined.csv")
if (!is.null(sc)) print(sc[, .(analysis, program, n_patients, pooled_effect, lo, hi, I2, A_positive, B_positive, spatial_label)])
# spatial_support_summary_main.csv feeds the evidence table (script 20) and figures;
# check it is identical to the main rows of the combined summary that script 23 writes
mn <- rd("06_results/spatial_support_combined/spatial_support_summary_main.csv")
if (!is.null(sc) && !is.null(mn)) {
  ref <- sc[analysis == "All spatial cohorts (main)", .(program, n_patients, pooled_effect, lo, hi, spatial_label)]
  cmp <- merge(ref, mn[, .(program, n_patients, pooled_effect, lo, hi, spatial_label)], by = "program",
               suffixes = c("_combined", "_main"), all = TRUE)
  same <- with(cmp, n_patients_combined == n_patients_main & abs(pooled_effect_combined - pooled_effect_main) < 1e-3 &
                    abs(lo_combined - lo_main) < 1e-3 & spatial_label_combined == spatial_label_main)
  cat("\n-- _main.csv vs combined summary (main rows):", if (all(same %in% TRUE)) "IDENTICAL" else "DIFFERENT - evidence table used stale spatial values", "--\n")
  if (!all(same %in% TRUE)) print(cmp[!(same %in% TRUE)])
}
cat("\n-- patients per cohort --\n")
ap <- rd("06_results/spatial_support_combined/criterionA_patients_combined.csv")
if (!is.null(ap)) print(ap[, .(patients=uniqueN(patient)), by=cohort])

sep("11. SENSITIVITY (expect 82/85 below 1; exceptions v2-only angio 1.12, TGFB 1.20; partition TGFB 1.08)")
ss <- rd("06_results/sensitivity_discovery/discovery_sensitivity_results.csv")
if (!is.null(ss)) { f <- ss[is.finite(MAE_ratio)]
  cat("fits:", nrow(f), "| below 1:", sum(f$MAE_ratio < 1), "| >=1:", sum(f$MAE_ratio >= 1), "\n")
  print(f[MAE_ratio >= 1, .(scenario, program, n_donors, MAE_ratio)])
  cat("\nfrozen genes re-selected, range:\n"); print(f[, .(range=paste(range(as.integer(sub("/.*","",frozen_genes_reselected))), collapse="-"))]) }

sep("12. EVIDENCE TABLE + CHIP TABLE (26 genes?)")
ev <- rd("08_tables/Table_program_evidence.csv"); if (!is.null(ev)) print(ev)
ch <- rd("08_tables/Table_chip_assay_proposal.csv")
if (!is.null(ch)) { cat("\nchip table rows:", nrow(ch), "\n"); print(ch[, .N, by=program]) }
cat("\n===== AUDIT COMPLETE =====\n")
