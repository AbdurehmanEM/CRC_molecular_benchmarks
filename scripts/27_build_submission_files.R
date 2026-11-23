# 27_build_submission_files.R - assemble everything needed for submission
#   1. widen the Figure 4 (validation forest) axis so no interval is clipped
#   2. rebuild Figures 1-6 from script 21
#   3. draw the two new supplementary figures (S3 marker check, S4 cohort-split spatial)
#   4. assemble Supplementary Tables S1-S10 into one workbook
#   5. copy main and supplementary figures into 10_submission/ with final names
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(openxlsx) })
SUB <- "10_submission"; dir.create(file.path(SUB, "figures"), recursive = TRUE, showWarnings = FALSE)

## ---- 1. axis fix -----------------------------------------------------------
fg <- "05_analysis_scripts/21_figures.R"; s <- readLines(fg)
n <- sum(grepl("limits = c(0.12, 8.5)", s, fixed = TRUE))
if (n) { s <- gsub("limits = c(0.12, 8.5)", "limits = c(0.10, 8.5)", s, fixed = TRUE)
         writeLines(s, fg); cat("Figure 4 axis widened to 0.10 (", n, "occurrence(s) )\n") } else
         cat("axis limits not found - skipped\n")

## ---- 2. rebuild figures ----------------------------------------------------
cat("\nRebuilding figures (this takes a few minutes)...\n"); source(fg)

## ---- 3. supplementary figures ---------------------------------------------
th <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 9))
sv <- function(p, nm, w, h) { ggsave(file.path(SUB, "figures", paste0(nm, ".pdf")), p, width = w, height = h, device = cairo_pdf)
                              ggsave(file.path(SUB, "figures", paste0(nm, ".png")), p, width = w, height = h, dpi = 300) }

# S3: RCTD marker-correlation check, every section
m1 <- fread("06_results/spatial_RCTD_marker_check.csv")
mk <- melt(m1, id.vars = "marker_set", variable.name = "column", value.name = "rho")
mk[, column := sub("^col[0-9]+_", "", column)]
p1 <- ggplot(mk, aes(column, marker_set, fill = rho)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.3) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "Spearman") +
  labs(x = "RCTD weight column", y = "Marker set",
       title = "A  Primary spatial cohort: mean marker-score vs RCTD weight (12 sections)") +
  th + theme(axis.text.x = element_text(angle = 45, hjust = 1))
m2 <- fread("06_results/spatial_support_extra/RCTD_marker_check_extra.csv")
p2 <- ggplot(m2, aes(section, type, fill = own_rho)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", own_rho)), size = 2.3) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "Spearman") +
  labs(x = NULL, y = "Cell type", title = "B  Added cohorts: own-marker correlation per section") +
  th + theme(axis.text.x = element_text(angle = 45, hjust = 1))
sv(patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 1)), "FigureS3_RCTD_marker_check", 7.2, 8)

# S4: spatial support by cohort group
S <- fread("06_results/spatial_support_combined/spatial_support_summary_combined.csv")[is.finite(pooled_effect)]
PL <- c(EPI_HYPOXIA = "Epithelial hypoxia", EPI_APOPTOSIS = "Epithelial apoptosis", STR_TGFB = "Stromal TGF-beta",
        MYE_IFNG = "Myeloid IFN-gamma", MYE_INFLAM = "Myeloid inflammatory",
        END_ANGIO = "Endothelial angiogenesis", END_TNFA_NFKB = "Endothelial TNF/NF-kB")
S[, plab := factor(PL[program], levels = rev(PL))]
S[, analysis := factor(analysis, levels = c("All spatial cohorts (main)", "Valdeolivas only", "Additional cohorts only"))]
p3 <- ggplot(S, aes(pooled_effect, plab, colour = analysis)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.3, colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.7,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 1.9, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c("#1b4965", "#5fa8d3", "#e07a5f"), name = NULL) +
  labs(x = "Composition-adjusted panel coefficient (criterion A), pooled with 95% CI", y = NULL,
       title = "Spatial support by cohort group") + th + theme(legend.position = "bottom")
sv(p3, "FigureS4_spatial_by_cohort", 7.2, 4.2)
cat("Supplementary figures S3 and S4 written\n")

## ---- 4. supplementary tables ----------------------------------------------
frz <- readRDS("04_processed_data/frozen_panel/frozen_panel_v2.rds")
panel_tbl <- rbindlist(lapply(names(frz$programs), function(pg) { P <- frz$programs[[pg]]
  data.table(program = pg, compartment = P$compartment, gene = P$panel$gene_symbol,
             ensembl = P$panel$row_id, coefficient = P$panel$coefficient,
             pred_mean = P$panel$pred_mean, pred_sd = P$panel$pred_sd, intercept = P$intercept,
             n_evaluation_genes = nrow(P$evaluation_genes)) }))
rd <- function(p) if (file.exists(p)) fread(p) else data.table(note = paste("file not found:", p))
stab <- { f <- list.files("06_results/panel", pattern = "stab", full.names = TRUE)
          if (length(f)) fread(f[1]) else data.table(note = "stability file not found") }
AU <- list.files("02_metadata", pattern = "eligibility_audit", full.names = TRUE)[1]
rdx <- function(sheet) if (!is.na(AU) && file.exists(AU) && sheet %in% openxlsx::getSheetNames(AU))
  as.data.table(openxlsx::read.xlsx(AU, sheet = sheet)) else data.table(note = paste("sheet not found:", sheet))
S1a <- rdx("Study_level"); S1b <- rdx("Sample_level"); S1c <- rdx("Paper_notes")
SHEETS <- list(
 "S1a_study_eligibility"  = S1a,
 "S1b_sample_eligibility" = S1b,
 "S1c_paper_notes"        = S1c,
 "S2_cohort_provenance"   = rd("02_metadata/cohort_provenance.csv"),
 "S3_software_versions"   = rd("06_results/software_versions.csv"),
 "S4_program_definitions" = rd("04_processed_data/programs/program_definitions.csv"),
 "S5_gene_partition"      = rd("04_processed_data/programs/gene_partition.tsv"),
 "S6_frozen_panel"        = panel_tbl,
 "S7_discovery_CV"        = rd("06_results/panel/cv_performance_discovery_v2.csv"),
 "S7b_selection_stability"= stab,
 "S8_validation_by_cohort"= rd("06_results/cross_cohort_meta/per_cohort_results_long.csv"),
 "S8b_meta_analysis"      = rd("06_results/cross_cohort_meta/meta_analysis_results.csv"),
 "S9_sensitivity"         = rd("06_results/sensitivity_discovery/discovery_sensitivity_results.csv"),
 "S10a_spatial_criterionA"= rd("06_results/spatial_support_combined/criterionA_patients_combined.csv"),
 "S10b_spatial_criterionB"= rd("06_results/spatial_support_combined/criterionB_patients_combined.csv"),
 "S10c_spatial_summary"   = rd("06_results/spatial_support_combined/spatial_support_summary_combined.csv"),
 "S10d_spatial_coverage"  = rd("06_results/spatial_support_extra/spatial_gene_coverage_extra.csv"))
wb <- createWorkbook()
for (nm in names(SHEETS)) { addWorksheet(wb, nm); writeData(wb, nm, SHEETS[[nm]]); freezePane(wb, nm, firstRow = TRUE) }
saveWorkbook(wb, file.path(SUB, "Supplementary_Tables.xlsx"), overwrite = TRUE)
cat("Supplementary_Tables.xlsx written with", length(SHEETS), "sheets\n")
print(data.table(sheet = names(SHEETS), rows = vapply(SHEETS, nrow, integer(1))))

## ---- 5. final figure names -------------------------------------------------
MAP <- c("Figure1_cohorts_coverage_updated"   = "Figure1",
         "Figure4_validation_forest_updated"  = "Figure2",
         "Figure5_spatial_support_updated"    = "Figure3",
         "Figure6_evidence_summary_updated"   = "Figure4",
         "Figure2_program_variation_updated"  = "FigureS1_program_variation",
         "Figure3_panel_baselines_updated"    = "FigureS2_panel_baselines")
for (old in names(MAP)) for (ext in c("pdf", "tiff", "png")) {
  src <- file.path("07_figures", paste0(old, ".", ext))
  if (file.exists(src)) file.copy(src, file.path(SUB, "figures", paste0(MAP[[old]], ".", ext)), overwrite = TRUE)
}
cat("\n=== 10_submission/figures ===\n"); print(list.files(file.path(SUB, "figures")))
cat("\nNOTE: submit the PDF versions (vector). Do not submit .docx or .jpg figures.\n")
