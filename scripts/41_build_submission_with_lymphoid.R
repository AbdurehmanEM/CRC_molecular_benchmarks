# 41_build_submission_with_lymphoid.R - rebuild figures, draw Figures S3-S4, copy final figure names (tables: script 44)
#   1. widen the Figure 4 (validation forest) axis so no interval is clipped
#   2. rebuild Figures 1-6 from script 21
#   3. draw the two new supplementary figures (S3 marker check, S4 cohort-split spatial)
#   4. assemble Supplementary Tables S1-S10 into one workbook
#   5. copy main and supplementary figures into 10_submission/ with final names
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(openxlsx) })
SUB <- "10_submission"; dir.create(file.path(SUB, "figures"), recursive = TRUE, showWarnings = FALSE)

fg <- "05_analysis_scripts/40_figures_with_lymphoid.R"
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
S <- rbind(fread("06_results/spatial_support_combined/spatial_support_summary_combined.csv"),
           fread("06_results/spatial_support_lymphoid_combined/spatial_support_summary_combined.csv"),
           fill = TRUE)[is.finite(pooled_effect)]
PL <- c(EPI_HYPOXIA = "Epithelial hypoxia", EPI_APOPTOSIS = "Epithelial apoptosis", STR_TGFB = "Stromal TGF-beta",
        MYE_IFNG = "Myeloid IFN-gamma", MYE_INFLAM = "Myeloid inflammatory",
        END_ANGIO = "Endothelial angiogenesis", END_TNFA_NFKB = "Endothelial TNF/NF-kB",
        LYM_IL2_STAT5 = "CD8 IL2/STAT5")
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
