# 45_revision_figures.R - revised publication figures (main 1-4, supplementary S1-S5)
# Implements the figure comments of the pre-submission review (reviewer figures 1-9).
# Based on 43_publication_figures.R (same visual system) and 41_build_submission_with_lymphoid.R (S3, S4).
# Reads saved results only; no new analysis. Writes NEW files; nothing existing is overwritten:
#   07_figures/revision_2026-11/<FigureX>_rev.{pdf,tiff,png}
# Run from the project root:  source("05_analysis_scripts/revision/45_revision_figures.R")
# ASCII-only file: non-ASCII characters are written as \u escapes.
#
# Reviewer figure -> output
#   1 coverage         -> Figure1_study_design_coverage_rev   (n/N in every tile, n.d. tiles, CD8 marked as amendment)
#   2 validation       -> Figure2_validation_forest_rev        (larger text, symbol key, endothelial n = 11, one rounding rule)
#   3 spatial          -> Figure3_spatial_support_rev          (scale bar, orientation, section ID, cohort legend, criteria defined)
#   4 scorecard        -> Figure4_evidence_scorecard_rev       (category legend, separate k / n column, 0.50 limit, hypoxia qualifier)
#   5 donor variation  -> FigureS1_program_variation_rev
#   6 panel baselines  -> FigureS2_panel_baselines_rev
#   7 RCTD markers     -> FigureS3_RCTD_marker_check_rev
#   8 spatial by group -> FigureS4_spatial_by_cohort_rev
#   9 reference sens.  -> FigureS5_reference_sensitivity_rev

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(scales) })

# ---------------------------------------------------------------- settings you may need to check
# Total patients per single-cell cohort (denominators in Figure 1B; from Table 1 / Figure 1 legend)
N_TOTAL <- c(GSE178341 = 62, GSE132465 = 23, GSE166555 = 12, GSE200997 = 16, GSE144735 = 6)
# Programs whose single-cell support is qualified as "not distinguishable from expression-matched
# random panels" (hypoxia option (a); script 53). The script prints the script-53 p-values as a check.
QUALIFIED <- c("EPI_HYPOXIA")
SEC <- "SN048_A416371_Rep1"; SHOW <- c("MYE_IFNG", "END_TNFA_NFKB")   # example section in Figure 3A
NULL53 <- "06_results/revision_2026-11/53_null_random_panels_summary.csv"

# ---------------------------------------------------------------- inputs
P <- function(path) {
  if (file.exists(path)) return(path)
  f <- file.path("figure_data", gsub("/", "__", path)); if (file.exists(f)) return(f)
  f <- file.path("figure_data", basename(path));        if (file.exists(f)) return(f)
  stop("Input not found: ", path)
}
rd <- function(path) fread(P(path))
OUT <- "07_figures/revision_2026-11"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------- visual system (as script 43)
FONT <- if (.Platform$OS.type == "windows") "sans" else {
  ff <- tryCatch(system("fc-list : family", intern = TRUE), error = function(e) character(0))
  if (any(grepl("Liberation Sans", ff))) "Liberation Sans" else "sans" }
INK <- "#1d1d1b"; INK2 <- "#55544f"; INK3 <- "#97968f"; GRID <- "#ecebe6"; BAND <- "#f6f5f1"; HIL <- "#fbeccd"
COMP_COL <- c(Epithelial = "#2a78d6", Stromal = "#eb6834", Myeloid = "#1baf7a", Endothelial = "#eda100", Lymphoid = "#6250d6")
EVID <- c("Independently supported with spatial evidence", "Independently supported in single-cell data", "Exploratory", "Unsupported")
EVID_SHORT <- setNames(c("Validated + spatial", "Validated (single-cell)", "Exploratory", "Unsupported"), EVID)
EVID_COL <- setNames(c("#0e5a34", "#3f9864", "#a8d5b5", "#dcdbd5"), EVID)
EVID_TXT <- setNames(c("white", "white", INK, INK), EVID)
DIV <- c("#1c5aa6", "#6f9fd8", "#e9e8e4", "#ee9a72", "#c0441a")
PLAB <- c(EPI_HYPOXIA = "Epithelial hypoxia", EPI_APOPTOSIS = "Epithelial apoptosis", STR_TGFB = "Stromal TGF-\u03b2",
          MYE_IFNG = "Myeloid IFN-\u03b3", MYE_INFLAM = "Myeloid inflammatory", END_ANGIO = "Endothelial angiogenesis",
          END_TNFA_NFKB = "Endothelial TNF/NF-\u03baB", LYM_IL2_STAT5 = "CD8 IL2/STAT5")
PORD <- names(PLAB)
COMP_OF <- c(EPI_HYPOXIA = "Epithelial", EPI_APOPTOSIS = "Epithelial", STR_TGFB = "Stromal", MYE_IFNG = "Myeloid",
             MYE_INFLAM = "Myeloid", END_ANGIO = "Endothelial", END_TNFA_NFKB = "Endothelial", LYM_IL2_STAT5 = "Lymphoid")
SP_COH <- c("Valdeolivas", "GSE280318", "GSE226997")
SP_PLAT <- c(Valdeolivas = "Visium", GSE280318 = "Visium HD", GSE226997 = "FFPE Visium")
SP_COL <- c(Valdeolivas = "#1b4965", GSE280318 = "#5fa8d3", GSE226997 = "#e07a5f")   # cohort colours (script 41)
SP_SHP <- c(Valdeolivas = 21, GSE280318 = 22, GSE226997 = 24)
W2 <- 7.0
EN <- "\u2013"; DAG <- "\u2020"; MID <- "\u00b7"

theme_pub <- function(base = 7.5) theme_minimal(base_size = base, base_family = FONT) + theme(
  text = element_text(colour = INK), axis.text = element_text(colour = INK, size = base - 0.5),
  axis.title = element_text(colour = INK, size = base - 0.2), panel.grid.minor = element_blank(),
  panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
  axis.line.x = element_line(colour = INK3, linewidth = 0.3), axis.ticks.x = element_line(colour = INK3, linewidth = 0.3),
  axis.ticks.length = unit(1.6, "pt"), strip.text = element_text(face = "bold", size = base - 0.2, hjust = 0, colour = INK),
  plot.tag = element_text(face = "bold", size = base + 3, colour = INK), plot.title = element_text(face = "bold", size = base + 0.2, colour = INK),
  plot.title.position = "plot", legend.position = "bottom", legend.text = element_text(size = base - 0.3, colour = INK),
  legend.title = element_text(size = base - 0.3, colour = INK), legend.key.size = unit(9, "pt"),
  legend.margin = margin(0, 0, 0, 0), plot.margin = margin(3, 5, 3, 3), panel.spacing = unit(3, "pt"))
CAP <- function(size = 6.3) theme(plot.caption = element_text(size = size, colour = INK, hjust = 0, family = FONT, lineheight = 1.1),
                                  plot.caption.position = "plot")
blank_y <- theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
blank_x <- theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.line.x = element_blank(), panel.grid.major.x = element_blank())
no_strip <- theme(strip.text = element_blank())
half_up <- function(v, d = 2) sign(v) * floor(abs(v) * 10^d + 0.5 + 1e-9) / 10^d
# One rounding rule everywhere (reviewer: "reconcile CI rounding"): values are first rounded to the
# 3 decimals shown in the Supplementary Tables, then half-up to 2 decimals (0.495 -> 0.50, as in the text).
r2 <- function(v) half_up(half_up(v, 3), 2)
fmt_ci <- function(e, l, h) ifelse(is.finite(e), sprintf("%.2f (%.2f%s%.2f)", r2(e), r2(l), EN, r2(h)), NA_character_)

SUFFIX <- "_rev"
save_fig <- function(p, name, w, h) {
  name <- paste0(name, SUFFIX)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".tiff")), p, width = w, height = h, dpi = 600, compression = "lzw", bg = "white")
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  message("saved ", name)
}

# ---------------------------------------------------------------- shared data
ev   <- rd("08_tables/Table_program_evidence_8.csv")
spA  <- rbind(rd("06_results/spatial_support_combined/criterionA_patients_combined.csv"),
              rd("06_results/spatial_support_lymphoid_combined/criterionA_patients_combined.csv"), fill = TRUE)
spB  <- rbind(rd("06_results/spatial_support_combined/criterionB_patients_combined.csv"),
              rd("06_results/spatial_support_lymphoid_combined/criterionB_patients_combined.csv"), fill = TRUE)
spAll <- rbind(rd("06_results/spatial_support_combined/spatial_support_summary_combined.csv"),
               rd("06_results/spatial_support_lymphoid_combined/spatial_support_summary_combined.csv"), fill = TRUE)
spS  <- spAll[analysis == "All spatial cohorts (main)"]
coh_long <- rbind(rd("06_results/cross_cohort_meta/per_cohort_results_long.csv"),
                  rd("06_results/validation_lymphoid/per_cohort_results_long_lymphoid.csv"), fill = TRUE)
coh_long[cohort == "GSE144735_KUL", cohort := "GSE144735"]
meta <- rbind(rd("06_results/cross_cohort_meta/meta_analysis_results.csv")[analysis == "main"],
              rd("06_results/validation_lymphoid/meta_analysis_results_lymphoid.csv"), fill = TRUE)
n_spatial_pat <- uniqueN(spA$patient)
SP_N   <- setNames(spA[, .(n = uniqueN(patient)), by = cohort][match(SP_COH, cohort), n], SP_COH)
SP_LAB <- setNames(sprintf("%s (%s; n = %d)", SP_COH, SP_PLAT[SP_COH], SP_N), SP_COH)
# patients evaluable per program (criterion A); counts shown as "positive / evaluable"
nA <- spA[is.finite(effect), .(nA = uniqueN(patient)), by = program]
nB <- spB[is.finite(diff),   .(nB = uniqueN(patient)), by = program]
cnt <- function(x, n) { x <- as.character(x); ifelse(grepl("/", x), x, ifelse(is.na(x) | x == "NA", "\u2014", paste0(x, "/", n))) }
spS <- merge(merge(spS, nA, by = "program", all.x = TRUE), nB, by = "program", all.x = TRUE)
spS[, `:=`(A_lab = cnt(A_positive, nA), B_lab = cnt(B_positive, nB))]

frz <- NULL
if (file.exists("04_processed_data/frozen_panel/frozen_panel_v2.rds")) {
  frz <- readRDS("04_processed_data/frozen_panel/frozen_panel_v2.rds")
  frz$programs <- c(frz$programs, readRDS("04_processed_data/frozen_panel/frozen_panel_LYM_v1.rds")$programs)
}
n_disc_elig <- if (!is.null(frz)) length(unique(unlist(lapply(frz$programs, `[[`, "donors")))) else 60L

# =========================================================================================
# FIGURE 1 - study workflow (A) and cohort x compartment coverage (B)
# Comments: discovery vs validation clearer; post hoc CD8 in the flow; n.d. distinct from 0;
# denominator for every tile.
# =========================================================================================
n_prog <- nrow(ev); n_genes <- sum(lengths(strsplit(ev$panel, ",\\s*")))
steps <- data.table(
  x = 1:6, grp = c("Discovery", rep("Discovery", 2), "Independent", "Independent", "Synthesis"),
  head = c("Discovery", "Programs", "Panels", "Validation", "Spatial", "Evidence"),
  body = c(sprintf("GSE178341\n%d of %d patients\n4 compartments\n+ CD8 T cells%s", n_disc_elig, N_TOTAL[["GSE178341"]], DAG),
           sprintf("%d programs retained\npanel and target\ngenes disjoint", n_prog),
           sprintf("elastic net, 1%s4 genes\n(%d genes in total)\nfrozen before\nany testing", EN, n_genes),
           "4 external\nsingle-cell cohorts\nfrozen panels,\nno refitting",
           sprintf("%d patients\n3 cohorts\ncomposition-\nadjusted", n_spatial_pat),
           "4 evidence\ncategories;\nchip assays\nproposed only"))
GRP_COL <- c(Discovery = "#3d3c38", Independent = "#1b4965", Synthesis = "#7a7972")
pA <- ggplot(steps) +
  annotate("rect", xmin = 0.55, xmax = 3.5, ymin = 1.0, ymax = 1.16, fill = GRP_COL[["Discovery"]], alpha = 0.12) +
  annotate("rect", xmin = 3.5, xmax = 5.5, ymin = 1.0, ymax = 1.16, fill = GRP_COL[["Independent"]], alpha = 0.15) +
  annotate("text", x = 2.02, y = 1.08, label = "Discovery data only (GSE178341)", size = 2.3, fontface = "bold", colour = INK, family = FONT) +
  annotate("text", x = 4.5, y = 1.08, label = "Independent patients, panels frozen", size = 2.3, fontface = "bold", colour = INK, family = FONT) +
  geom_segment(data = steps[x < 6], aes(x = x + 0.43, xend = x + 0.57, y = 0.5, yend = 0.5), colour = INK2, linewidth = 0.4,
               arrow = arrow(length = unit(3, "pt"), type = "closed")) +
  geom_rect(aes(xmin = x - 0.42, xmax = x + 0.42, ymin = 0.03, ymax = 0.95), fill = BAND, colour = NA) +
  geom_rect(aes(xmin = x - 0.42, xmax = x + 0.42, ymin = 0.79, ymax = 0.95, fill = grp), colour = NA) +
  scale_fill_manual(values = GRP_COL, guide = "none") +
  geom_text(aes(x = x, y = 0.87, label = head), colour = "white", fontface = "bold", size = 2.4, family = FONT) +
  geom_text(aes(x = x, y = 0.42, label = body), colour = INK, size = 2.15, lineheight = 1.0, family = FONT) +
  scale_x_continuous(limits = c(0.55, 6.45), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1.17), expand = c(0, 0)) +
  labs(tag = "A") + theme_void(base_family = FONT) + theme(plot.tag = element_text(face = "bold", size = 10.5), plot.margin = margin(3, 3, 6, 3))

el  <- rd("06_results/GSE178341/pseudobulk_eligibility.csv")[min_cells == 20]
elL <- rd("06_results/GSE178341/pseudobulk_eligibility_lymphoid.csv")[min_cells == 20]
xtra <- rd("06_results/validation_lymphoid/validation_coverage_extra_cohorts.csv")
cov_n <- function(path, col = "n_donors") { x <- rd(path); x[, compartment := COMP_OF[program]]; x[, .(n = max(get(col), na.rm = TRUE)), by = compartment] }
cv1 <- rbind(
  data.table(cohort = "GSE178341", role = "Discovery", compartment = names(COMP_COL),
             n = as.numeric(c(unlist(el[, setdiff(names(COMP_COL), "Lymphoid"), with = FALSE]), elL$Lymphoid_CD8))),
  cbind(cohort = "GSE132465", role = "Validation", rbind(cov_n("06_results/validation_GSE132465/validation_gene_coverage.csv"),
        data.table(compartment = "Lymphoid", n = rd("06_results/validation_GSE132465/validation_gene_coverage_lymphoid.csv")$n_donors[1]))),
  cbind(cohort = "GSE166555", role = "Validation", rbind(cov_n("06_results/validation_GSE166555/validation_gene_coverage.csv"),
        data.table(compartment = "Lymphoid", n = xtra[cohort == "GSE166555", n_donors][1]))),
  cbind(cohort = "GSE200997", role = "Validation", rbind(cov_n("06_results/validation_GSE200997/validation_gene_coverage.csv"),
        data.table(compartment = "Lymphoid", n = NA_real_))),
  cbind(cohort = "GSE144735", role = "Descriptive", rbind(cov_n("06_results/validation_GSE144735/KUL_gene_coverage.csv", "kul_donors"),
        data.table(compartment = "Lymphoid", n = xtra[cohort == "GSE144735_KUL", n_donors][1]))),
  use.names = TRUE)
cv1[, N := N_TOTAL[cohort]]
spc <- copy(spA)[is.finite(effect)][, compartment := COMP_OF[program]][, .(n = uniqueN(patient)), by = .(cohort, compartment)]
spc <- CJ(cohort = SP_COH, compartment = names(COMP_COL))[spc, n := as.numeric(i.n), on = .(cohort, compartment)]
spc[, N := SP_N[cohort]]
cv1 <- rbind(cv1, spc[, .(cohort, role = "Spatial (patients)", compartment, n, N)], use.names = TRUE)
ROLE_LAB <- c(Discovery = "Discovery", Validation = "Independent\nvalidation", Descriptive = "Descriptive\nonly", `Spatial (patients)` = "Spatial\n(patients)")
cv1[, `:=`(role = factor(ROLE_LAB[as.character(role)], levels = ROLE_LAB), compartment = factor(compartment, levels = names(COMP_COL)))]
PLAT <- setNames(sprintf("%s (%s)", SP_COH, SP_PLAT[SP_COH]), SP_COH)
cv1[cohort %in% names(PLAT), cohort := PLAT[cohort]]
cv1[, cohort := factor(cohort, levels = rev(c("GSE178341", "GSE132465", "GSE166555", "GSE200997", "GSE144735", PLAT)))]
cv1[, `:=`(nd = is.na(n), alpha = ifelse(is.na(n), 0, 0.15 + 0.85 * pmin(n, 60) / 60),
           lab = ifelse(is.na(n), "n.d.", sprintf("%d/%d", as.integer(n), as.integer(N))),
           low = !is.na(n) & n < 10 & role %in% ROLE_LAB[c("Validation", "Descriptive")])]
cv1[, txt := ifelse(!is.na(n) & n >= 30, "white", INK)]
XLAB <- c(Epithelial = "Epithelial", Stromal = "Stromal", Myeloid = "Myeloid", Endothelial = "Endothelial", Lymphoid = paste0("CD8 T cells", DAG))
pB <- ggplot(cv1, aes(compartment, cohort)) +
  geom_tile(aes(fill = compartment, alpha = alpha), colour = "white", linewidth = 0.9, width = 0.97, height = 0.92) +
  geom_tile(data = cv1[nd == TRUE], fill = "white", colour = INK3, linewidth = 0.3, width = 0.9, height = 0.8) +
  geom_tile(data = cv1[low == TRUE], fill = NA, colour = INK2, linewidth = 0.35, linetype = "22", width = 0.86, height = 0.78) +
  geom_text(aes(label = lab, colour = txt, fontface = ifelse(nd, "italic", "bold")), size = 2.35, family = FONT) +
  geom_vline(xintercept = 4.5, colour = INK2, linewidth = 0.35, linetype = "22") +
  scale_fill_manual(values = COMP_COL, guide = "none") + scale_alpha_identity() + scale_colour_identity() +
  scale_x_discrete(position = "top", expand = c(0, 0), labels = XLAB) + scale_y_discrete(expand = c(0, 0)) +
  facet_grid(role ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL, tag = "B",
       caption = paste0("Tiles: eligible donors / patients in the cohort (\u226520 cells per donor \u00d7 compartment, tumour specimens); spatial rows: patients evaluable / patients in the cohort.\n",
                        "Dashed outline: <10 donors, not evaluable for validation. n.d. (white tile): not definable because GSE200997 annotates T and NK cells jointly; this is not an absence of cells.\n",
                        DAG, " CD8 T-cell compartment added by protocol amendment after the four-compartment analysis (right of the dashed line).")) +
  theme_pub() + CAP(6.2) +
  theme(panel.grid = element_blank(), panel.grid.major.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(face = "bold", colour = INK, size = 7), axis.text.y = element_text(colour = INK, size = 6.8),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", colour = INK, size = 6.6, lineheight = 0.95),
        panel.spacing.y = unit(7, "pt"), strip.background = element_rect(fill = BAND, colour = NA))
save_fig(wrap_elements(full = pA) / pB + plot_layout(heights = c(1.05, 2.6)), "Figure1_study_design_coverage", W2, 5.2)

# =========================================================================================
# FIGURE 2 - frozen-panel validation forest
# Comments: too dense (enlarged, taller, larger text); distinguish pooled diamonds from single
# cohorts (explicit key); label n = 11 endothelial directly; reconcile CI rounding (r2()).
# The within-cohort (secondary) panel is also saved alone as Figure S6 in case the journal
# prefers frozen-only in the main figure (set MAIN_BOTH <- FALSE).
# =========================================================================================
MAIN_BOTH <- TRUE
COH <- c("GSE132465", "GSE166555", "GSE200997", "GSE144735")
fr <- CJ(program = PORD, scaling = c("frozen", "within"), cohort = c(COH, "Pooled"), sorted = FALSE)
cl <- coh_long[, .(program, scaling, cohort, n_donors, MAE_ratio, ratio_lo, ratio_hi, status)]
mt <- meta[, .(program, scaling, cohort = "Pooled", n_donors = total_donors, MAE_ratio = pooled_MAE_ratio, ratio_lo, ratio_hi,
               status = ifelse(is.finite(pooled_MAE_ratio), "pooled", "not pooled"))]
fr <- merge(fr, rbind(cl, mt, fill = TRUE), by = c("program", "scaling", "cohort"), all.x = TRUE, sort = FALSE)
fr[, y := 6 - match(cohort, c(COH, "Pooled"))]
fr[, `:=`(plab = factor(PLAB[program], levels = PLAB), comp = COMP_OF[program], ok = is.finite(MAE_ratio) & MAE_ratio > 0)]
sing <- fr[ok & !cohort %in% c("Pooled", "GSE144735") & n_donors >= 10, .(nk = .N, s_coh = cohort[1], s_n = n_donors[1]), by = .(program, scaling)]
fr <- merge(fr, sing, by = c("program", "scaling"), all.x = TRUE, sort = FALSE)
fr[, note := fifelse(ok, NA_character_,
             fifelse(grepl("^Excluded", status), "excluded (T/NK not separable)",
             fifelse(cohort == "Pooled" & !is.na(nk) & nk == 1, sprintf("not pooled: single\ncohort, n = %d", as.integer(s_n)),
             fifelse(cohort == "Pooled", "not pooled (<2 cohorts)", "not evaluable"))))]
fr[, est_txt := fifelse(ok, fmt_ci(MAE_ratio, ratio_lo, ratio_hi), "\u2014")]
fr[, desc := cohort == "GSE144735"]
fr[, rowlab := fifelse(cohort == "Pooled", "Pooled (RE)", sprintf("%s  %s", cohort, fifelse(is.finite(n_donors) & n_donors > 0, paste0("n = ", n_donors), "")))]
fr[, single_hl := ok & !is.na(nk) & nk == 1 & cohort == s_coh & !cohort %in% c("Pooled", "GSE144735")]
XL <- c(0.1, 10); TS <- 2.3
forest <- function(sc, title) {
  d <- fr[scaling == sc]; dp <- d[ok == TRUE & cohort == "Pooled"]; dc <- d[ok == TRUE & cohort != "Pooled"]
  dia <- dp[, .(x = c(ratio_lo, MAE_ratio, ratio_hi, MAE_ratio), yy = c(y, y + 0.36, y, y - 0.36)), by = .(plab, program, comp)]
  ggplot(d) +
    annotate("rect", xmin = XL[1], xmax = 1, ymin = -Inf, ymax = Inf, fill = BAND) +
    geom_vline(xintercept = 1, colour = INK2, linewidth = 0.3, linetype = "22") +
    geom_linerange(data = dc, aes(y = y, xmin = pmax(ratio_lo, XL[1]), xmax = pmin(ratio_hi, XL[2]), colour = comp, alpha = desc), linewidth = 0.6) +
    geom_point(data = dc, aes(y = y, x = MAE_ratio, fill = comp, alpha = desc), shape = 21, colour = "white", stroke = 0.3, size = 2.0) +
    geom_polygon(data = dia, aes(x = x, y = yy, group = program), fill = INK, colour = INK, linewidth = 0.2) +
    geom_text(data = d[ok == FALSE], aes(y = y, x = 1, label = note), colour = INK2, size = 1.95, fontface = "italic", family = FONT, lineheight = 0.9) +
    scale_colour_manual(values = COMP_COL, guide = "none") + scale_fill_manual(values = COMP_COL, guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.38), guide = "none") +
    scale_x_log10(limits = XL, breaks = c(0.25, 1, 4), labels = c("0.25", "1", "4"), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.4, 5.6), expand = c(0, 0)) +
    facet_grid(plab ~ ., switch = "y") + labs(x = "MAE ratio vs discovery mean (log scale)", y = NULL, title = title) +
    theme_pub() + blank_y + no_strip + theme(plot.title = element_text(hjust = 0, size = 7.3), panel.spacing = unit(5, "pt"))
}
txtcol <- function(sc, title) ggplot(fr[scaling == sc]) +
  geom_text(aes(x = 0, y = y, label = est_txt, fontface = ifelse(cohort == "Pooled" | single_hl, "bold", "plain"),
                colour = ifelse(desc, INK3, INK)), hjust = 0, size = TS, family = FONT) +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = c(0.4, 5.6), expand = c(0, 0)) +
  facet_grid(plab ~ .) + labs(x = " ", y = NULL, title = title) + theme_pub() + blank_y + no_strip +
  theme(axis.text.x = element_text(colour = "white"), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(), panel.grid.major.x = element_blank(), panel.spacing = unit(5, "pt"),
        plot.title = element_text(size = 7.0, colour = INK, face = "bold", hjust = 0))
labcol <- ggplot(fr[scaling == "frozen"]) +
  geom_text(aes(x = 1, y = y, label = rowlab, fontface = ifelse(cohort == "Pooled", "bold", "plain"),
                colour = ifelse(desc, INK3, INK)), hjust = 1, size = TS, family = FONT) +
  scale_colour_identity() + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = c(0.4, 5.6), expand = c(0, 0)) +
  facet_grid(plab ~ ., switch = "y") + labs(x = " ", y = NULL, title = " \n ") + theme_pub() + blank_y +
  theme(plot.title = element_text(size = 7.3), axis.text.x = element_text(colour = "white"), axis.line.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(),
        panel.grid.major.x = element_blank(), panel.spacing = unit(5, "pt"),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 7.2, colour = INK))
# symbol key (distinguishes single-cohort estimates, descriptive cohort and pooled diamonds)
key <- ggplot() +
  annotate("segment", x = 0.10, xend = 0.16, y = 0.5, yend = 0.5, colour = INK2, linewidth = 0.6) +
  annotate("point", x = 0.13, y = 0.5, shape = 21, fill = INK2, colour = "white", size = 2.2) +
  annotate("text", x = 0.17, y = 0.5, label = "Single cohort (95% CI)", hjust = 0, size = 2.1, family = FONT, colour = INK) +
  annotate("segment", x = 0.38, xend = 0.44, y = 0.5, yend = 0.5, colour = INK2, alpha = 0.38, linewidth = 0.6) +
  annotate("point", x = 0.41, y = 0.5, shape = 21, fill = INK2, colour = "white", alpha = 0.38, size = 2.2) +
  annotate("text", x = 0.45, y = 0.5, label = "<10 donors, descriptive", hjust = 0, size = 2.1, family = FONT, colour = INK) +
  annotate("polygon", x = c(0.655, 0.67, 0.685, 0.67), y = c(0.5, 0.95, 0.5, 0.05), fill = INK) +
  annotate("text", x = 0.695, y = 0.5, label = "Pooled RE estimate (95% CI)", hjust = 0, size = 2.1, family = FONT, colour = INK) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + theme_void()
cap2 <- paste0("MAE, mean absolute error. Points: single cohorts with 95% percentile CI (2,000 donor bootstrap resamples). Diamonds: DerSimonian", EN, "Laird\n",
               "random-effects (RE) pooled estimate, computed only where \u22652 cohorts had \u226510 donors. Both endothelial programs are single-cohort estimates\n",
               "(GSE132465, n = 11; bold) and were not pooled. Shaded: ratio <1 (panel better than the discovery-mean baseline). Not evaluable: a panel gene missing,\n",
               "<80% evaluation-gene coverage or <10 donors. Excluded: T and NK cells not separable in GSE200997. Values rounded to 3, then 2 decimals (half-up), as in the text.")
if (MAIN_BOTH) {
  top <- wrap_plots(labcol, forest("frozen", "Frozen discovery scaling\n(primary)"), txtcol("frozen", "MAE ratio\n(95% CI)"),
                    forest("within", "Within-cohort scaling\n(secondary)"), txtcol("within", "MAE ratio\n(95% CI)"),
                    nrow = 1, widths = c(1.45, 1.7, 1.15, 1.7, 1.15))
  p2 <- (top / key) + plot_layout(heights = c(1, 0.02)) + plot_annotation(caption = cap2, theme = CAP(6.0))
  save_fig(p2, "Figure2_validation_forest", W2, 9.6)
} else {
  top <- wrap_plots(labcol, forest("frozen", "Frozen discovery scaling\n(primary)"), txtcol("frozen", "MAE ratio\n(95% CI)"),
                    nrow = 1, widths = c(1.45, 2.4, 1.25))
  p2 <- (top / key) + plot_layout(heights = c(1, 0.02)) + plot_annotation(caption = cap2, theme = CAP(6.0))
  save_fig(p2, "Figure2_validation_forest", W2, 9.6)
}
s6 <- (wrap_plots(labcol, forest("within", "Within-cohort scaling\n(secondary)"), txtcol("within", "MAE ratio\n(95% CI)"),
                  nrow = 1, widths = c(1.45, 2.4, 1.25)) / key) + plot_layout(heights = c(1, 0.02)) +
  plot_annotation(caption = cap2, theme = CAP(6.0))
save_fig(s6, "FigureS6_validation_within_cohort", W2, 9.6)

# =========================================================================================
# FIGURE 3 - spatial support
# Comments: A needs orientation, scale bar and section ID; B/C need an exact legend for shapes,
# dots and diamonds, criterion definitions and the patient-level bootstrap; avoid functional reading.
# =========================================================================================
spot_csv <- file.path("figure_data", paste0("spot_scores_", SEC, ".csv"))
sp <- if (file.exists(spot_csv)) fread(spot_csv) else readRDS("04_processed_data/Valdeolivas_Visium/spot_program_scores.rds")[section == SEC]
sp <- sp[program %in% SHOW]
rho <- sp[, .(r = cor(target, panel, method = "spearman")), by = program]
# scale bar: Visium spot centres are 100 um apart, so the median nearest-neighbour distance in pixels = 100 um
xy <- unique(sp[, .(pxl_col, pxl_row)]); set.seed(1); s_i <- sample(nrow(xy), min(nrow(xy), 400))
nn <- vapply(s_i, function(i) { dd <- sqrt((xy$pxl_col - xy$pxl_col[i])^2 + (xy$pxl_row - xy$pxl_row[i])^2); min(dd[dd > 0]) }, numeric(1))
px_per_um <- median(nn) / 100; bar_px <- 1000 * px_per_um
message(sprintf("Figure 3A scale: %.3f px per um (median spot spacing %.1f px = 100 um)", px_per_um, median(nn)))
m <- melt(sp, id.vars = c("program", "pxl_row", "pxl_col"), measure.vars = c("target", "panel"), variable.name = "score")
m[, value := ave(value, program, score, FUN = function(v) as.numeric(scale(v)))]
lev3 <- sprintf("%s\nSpearman \u03c1 = %.2f (spots)", PLAB[SHOW], rho$r[match(SHOW, rho$program)])
m[, plab := factor(sprintf("%s\nSpearman \u03c1 = %.2f (spots)", PLAB[program], rho$r[match(program, rho$program)]), levels = lev3)]
m[, score := factor(ifelse(score == "target", "Target (evaluation genes)", "Frozen panel score"), levels = c("Target (evaluation genes)", "Frozen panel score"))]
xr <- range(sp$pxl_col); yr <- range(-sp$pxl_row)
sb <- data.table(score = factor("Frozen panel score", levels = levels(m$score)), plab = factor(lev3[1], levels = lev3),
                 x0 = xr[1], x1 = xr[1] + bar_px, y = yr[1] - 0.05 * diff(yr))
pat <- sub("^SN[0-9]+_(A[0-9]+)_.*$", "\\1", SEC)
p3a <- ggplot(m, aes(pxl_col, -pxl_row, colour = pmin(pmax(value, -2.5), 2.5))) + geom_point(size = 0.42, stroke = 0) +
  geom_segment(data = sb, aes(x = x0, xend = x1, y = y, yend = y), inherit.aes = FALSE, colour = INK, linewidth = 0.8) +
  geom_text(data = sb, aes(x = (x0 + x1) / 2, y = y - 0.05 * diff(yr), label = "1 mm"), inherit.aes = FALSE, size = 2.2, colour = INK, family = FONT) +
  facet_grid(score ~ plab) + coord_equal(clip = "off") +
  scale_colour_gradientn(colours = DIV, limits = c(-2.5, 2.5), breaks = c(-2, 0, 2), name = "z-score\nwithin section") +
  labs(x = NULL, y = NULL, tag = "A",
       title = sprintf("Valdeolivas Visium section %s (patient %s; one of two sections)", SEC, pat),
       subtitle = "Spots drawn at their image coordinates, in the orientation of the H&E image (top of panel = top of image). Scale bar 1 mm.") +
  theme_pub() + theme(axis.text = element_blank(), axis.line.x = element_blank(), axis.ticks = element_blank(), axis.ticks.x = element_blank(),
                      panel.grid = element_blank(), panel.grid.major.x = element_blank(), plot.subtitle = element_text(size = 6.4, colour = INK),
                      legend.position = "right", legend.key.height = unit(16, "pt"), legend.key.width = unit(6, "pt"),
                      strip.text.x = element_text(hjust = 0.5, lineheight = 1.05), strip.text.y = element_text(angle = -90, hjust = 0.5, face = "bold", colour = INK))
POOL <- "Pooled random-effects estimate, 95% CI"
GRP_LAB <- c(SP_LAB, Pooled = POOL)
crit_plot <- function(dat, vcol, summ, xlab, title, tag, cnt_col, add_ci = TRUE) {
  dat <- copy(dat)[is.finite(get(vcol))][, `:=`(plab = factor(PLAB[program], levels = rev(PLAB)), grp = factor(cohort, levels = c(SP_COH, "Pooled")))]
  summ <- copy(summ)[, `:=`(plab = factor(PLAB[program], levels = rev(PLAB)), grp = factor("Pooled", levels = c(SP_COH, "Pooled")))]
  g <- ggplot(dat, aes(y = plab)) + geom_vline(xintercept = 0, colour = INK2, linewidth = 0.35, linetype = "22") +
    geom_point(aes(x = .data[[vcol]], fill = grp, shape = grp, colour = grp), stroke = 0.3, size = 1.7, alpha = 0.85,
               position = position_jitter(height = 0.17, width = 0, seed = 3))
  if (add_ci) g <- g + geom_linerange(data = summ, aes(xmin = lo, xmax = hi), colour = INK, linewidth = 0.8) +
    geom_point(data = summ, aes(x = pooled_effect, fill = grp, shape = grp, colour = grp), stroke = 0.6, size = 2.5)
  g + geom_text(data = summ, aes(x = Inf, label = .data[[cnt_col]]), hjust = 1.05, size = 2.2, colour = INK, family = FONT) +
    scale_fill_manual(values = c(SP_COL, Pooled = "white"), labels = GRP_LAB, name = NULL, drop = FALSE) +
    scale_colour_manual(values = c(setNames(rep("white", 3), SP_COH), Pooled = INK), labels = GRP_LAB, name = NULL, drop = FALSE) +
    scale_shape_manual(values = c(SP_SHP, Pooled = 23), labels = GRP_LAB, name = NULL, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.2))) +
    labs(x = xlab, y = NULL, tag = tag, title = title) + theme_pub() +
    theme(plot.title = element_text(size = 7, face = "bold")) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.4, alpha = 1, stroke = 0.6)))
}
p3b <- crit_plot(spA[, .(cohort, patient, program, effect)], "effect", spS,
                 "Composition-adjusted panel coefficient", "Criterion A (right: patients positive / evaluable)", "B", "A_lab") +
  theme(axis.text.y = element_text(colour = INK, size = 6.8))
p3c <- crit_plot(spB[, .(cohort, patient, program, diff)], "diff", spS,
                 "Target difference, expected-high minus comparison region", "Criterion B (right: patients positive / evaluable)", "C", "B_lab", add_ci = FALSE) +
  theme(axis.text.y = element_blank()) + guides(fill = "none", shape = "none", colour = "none")
cap3 <- paste0("B, C: one point per patient (colour and shape = cohort); ", n_spatial_pat, " patients in three cohorts, fewer for angiogenesis and IFN-\u03b3 (probe-panel coverage).\n",
               "Criterion A: spot-level regression of the program target on the frozen panel score, adjusted for RCTD-estimated proportions of nine cell types, averaged over\n",
               "a patient's sections; patient-level uncertainty from a spatial block bootstrap (8 \u00d7 8-spot blocks, 1,000 resamples). Open diamond: random-effects estimate\n",
               "pooled over patients with 95% CI. Criterion B: target in expected-high minus comparison regions (pathologist regions in Valdeolivas; upper vs lower quartile of\n",
               "compartment proportion elsewhere and for CD8 in all cohorts). Spatial association shows co-location of transcriptional programs; it is not evidence of function,\n",
               "cell communication or causality.")
p3 <- p3a / ((p3b | p3c) + plot_layout(guides = "collect") & theme(legend.position = "bottom")) + plot_layout(heights = c(1.25, 1)) +
  plot_annotation(caption = cap3, theme = CAP(6.0))
save_fig(p3, "Figure3_spatial_support", W2, 8.6)

# =========================================================================================
# FIGURE 4 - evidence scorecard
# Comments: evidence categories with a clear legend; separate k / n column; standardise the 0.44
# interval (0.38-0.50, not 0.49); remove caption text about symbol shapes that do not exist;
# hypoxia qualifier (option (a)).
# =========================================================================================
sc_rows <- rbindlist(lapply(PORD, function(pg) {
  mm <- meta[scaling == "frozen" & program == pg & is.finite(pooled_MAE_ratio)]
  cc <- coh_long[scaling == "frozen" & program == pg & is.finite(MAE_ratio) & MAE_ratio > 0 & n_donors >= 10 & cohort != "GSE144735" & !grepl("^Excluded", status)]
  if (nrow(mm)) {
    k <- if ("k_cohorts" %in% names(mm) && is.finite(mm$k_cohorts[1])) mm$k_cohorts[1] else nrow(cc)
    n <- if ("total_donors" %in% names(mm) && is.finite(mm$total_donors[1])) mm$total_donors[1] else sum(cc$n_donors)
    data.table(program = pg, sc_est = mm$pooled_MAE_ratio[1], sc_lo = mm$ratio_lo[1], sc_hi = mm$ratio_hi[1], k = k, n = n, sc_src = "pooled")
  } else if (nrow(cc) == 1) {
    data.table(program = pg, sc_est = cc$MAE_ratio, sc_lo = cc$ratio_lo, sc_hi = cc$ratio_hi, k = 1L, n = cc$n_donors, sc_src = cc$cohort)
  } else data.table(program = pg, sc_est = NA_real_, sc_lo = NA_real_, sc_hi = NA_real_, k = nrow(cc), n = sum(cc$n_donors), sc_src = "none")
}))
e4 <- merge(copy(ev), sc_rows, by = "program", all.x = TRUE)
e4[, sc_txt := fmt_ci(sc_est, sc_lo, sc_hi)]
chk <- e4[, .(program, table_value = frozen_MAE_ratio_CI, figure_value = gsub(EN, "-", sc_txt))]
cat("\nFigure 4 single-cell values (table string vs value recomputed with the common rounding rule):\n"); print(chk)
e4 <- merge(e4, spS[, .(program, sp_est = pooled_effect, sp_lo = lo, sp_hi = hi, A_lab, B_lab)], by = "program", all.x = TRUE)
# hypoxia qualifier: fixed list (manuscript decision) checked against script 53
if (file.exists(NULL53)) {
  q53 <- fread(NULL53)[null == "A fixed random panels" & (setting == "mean log ratio over cohorts" | grepl("^validation", setting))]
  cat("\nScript 53 (fixed random panels) p-values, for the qualifier check:\n"); print(q53[, .(program, setting, observed, null_median, p_empirical)])
  pm <- q53[setting == "mean log ratio over cohorts", .(program, p = p_empirical)]
  flag <- intersect(pm[p >= 0.05, program], e4[sc_independent == TRUE, program])
  if (!setequal(flag, QUALIFIED)) message("NOTE: script 53 pooled p >= 0.05 for supported programs: ", paste(flag, collapse = ", "),
                                          "; QUALIFIED is set to ", paste(QUALIFIED, collapse = ", "), ". Check before submission.")
} else message("NOTE: ", NULL53, " not found; using QUALIFIED = ", paste(QUALIFIED, collapse = ", "))
e4[, qual := program %in% QUALIFIED]
e4[, `:=`(evidence_category = factor(evidence_category, levels = EVID), comp = factor(COMP_OF[program], levels = names(COMP_COL)))]
setorder(e4, evidence_category, comp, program)
e4[, y := rev(seq_len(.N))]
e4[, kn := fifelse(is.na(sc_est), "\u2014", sprintf("%d / %d", as.integer(k), as.integer(n)))]
e4[, genes := vapply(strsplit(panel, ",\\s*"), function(g) if (length(g) > 2) paste(paste(g[1:2], collapse = "  "), paste(g[-(1:2)], collapse = "  "), sep = "\n") else paste(g, collapse = "  "), "")]
EVID_2L <- setNames(c("Validated\n+ spatial", "Validated\n(single-cell)", "Exploratory", "Unsupported"), EVID)
e4[, ev_lab := ifelse(qual, paste0(EVID_2L[as.character(evidence_category)], "*"), EVID_2L[as.character(evidence_category)])]
YL <- c(0.5, nrow(e4) + 0.5); FS <- 2.15
base4 <- function() list(scale_y_continuous(limits = YL, expand = c(0, 0)), theme_pub(), blank_y,
                         theme(plot.title = element_text(size = 6.8, hjust = 0, lineheight = 0.95)))
txt_only <- theme(axis.text.x = element_text(colour = "white"), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                  panel.grid = element_blank(), panel.grid.major.x = element_blank())
c1 <- ggplot(e4) + geom_rect(aes(xmin = 0, xmax = 0.04, ymin = y - 0.38, ymax = y + 0.38, fill = comp)) +
  geom_text(aes(x = 0.06, y = y, label = PLAB[program]), hjust = 0, size = 2.0, fontface = "bold", colour = INK, family = FONT) +
  scale_fill_manual(values = COMP_COL, name = "Compartment", drop = FALSE) + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Program", x = " ") + base4() + txt_only
c2 <- ggplot(e4, aes(y = y)) + annotate("rect", xmin = 0.1, xmax = 1, ymin = -Inf, ymax = Inf, fill = BAND) +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.3, linetype = "22") +
  geom_linerange(aes(xmin = pmax(sc_lo, 0.1), xmax = pmin(sc_hi, 10), colour = comp), linewidth = 0.8) +
  geom_point(aes(x = sc_est, colour = comp), fill = "white", shape = 21, size = 2.0, stroke = 0.7) +
  scale_colour_manual(values = COMP_COL, guide = "none") +
  scale_x_log10(limits = c(0.1, 10), breaks = c(0.25, 1, 4), labels = c("0.25", "1", "4"), expand = c(0, 0)) +
  labs(title = "Single-cell\n(frozen scaling)", x = "MAE ratio") + base4()
c2t <- ggplot(e4) + geom_text(aes(x = 0, y = y, label = ifelse(is.na(sc_txt), "\u2014", sc_txt)), hjust = 0, size = FS, colour = INK, family = FONT) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + labs(title = "MAE ratio\n(95% CI)", x = " ") + base4() + txt_only
c2k <- ggplot(e4) + geom_text(aes(x = 0.5, y = y, label = kn), size = FS, colour = INK, family = FONT) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + labs(title = "k / n", x = " ") + base4() + txt_only
cW <- ggplot(e4) + geom_text(aes(x = 0.5, y = y, label = ifelse(as.logical(within_scaling_supported), "Supported", "\u2014")), size = FS, colour = INK, family = FONT) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + labs(title = "Within-\ncohort", x = " ") + base4() + txt_only
c3 <- ggplot(e4, aes(y = y)) + geom_vline(xintercept = 0, colour = INK2, linewidth = 0.3, linetype = "22") +
  geom_linerange(aes(xmin = sp_lo, xmax = sp_hi, colour = comp), linewidth = 0.8) +
  geom_point(aes(x = sp_est, colour = comp), fill = "white", shape = 23, size = 2.1, stroke = 0.7) +
  scale_colour_manual(values = COMP_COL, guide = "none") + scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  labs(title = "Spatial\n(criterion A)", x = "Adjusted coefficient") + base4()
c3t <- ggplot(e4) + geom_text(aes(x = 0.5, y = y, label = paste0(A_lab, "\n", B_lab)), size = FS, colour = INK, family = FONT, lineheight = 0.95) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + labs(title = "Patients\nA / B", x = " ") + base4() + txt_only
c4 <- ggplot(e4) + geom_tile(aes(x = 0.5, y = y, fill = evidence_category), width = 0.96, height = 0.78) +
  geom_text(aes(x = 0.5, y = y, label = ev_lab, colour = EVID_TXT[as.character(evidence_category)]), size = 1.95, fontface = "bold", family = FONT, lineheight = 0.9) +
  scale_fill_manual(values = EVID_COL, labels = EVID, name = "Evidence category", drop = FALSE) +
  scale_colour_identity() + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Evidence\ncategory", x = " ") + base4() + txt_only
c5 <- ggplot(e4) + geom_text(aes(x = 0, y = y, label = genes), hjust = 0, size = 2.0, fontface = "italic", colour = INK, family = FONT, lineheight = 0.95) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + labs(title = "Panel genes", x = " ") + base4() + txt_only
p4 <- (c1 + c2 + c2t + c2k + cW + c3 + c3t + c4 + c5 + plot_layout(nrow = 1, widths = c(1.45, 0.9, 0.92, 0.5, 0.62, 0.85, 0.5, 0.85, 1.1), guides = "collect")) &
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(2, "pt"))
p4 <- p4 + plot_annotation(caption = paste0(
  "Single-cell: random-effects pooled MAE ratio (frozen scaling) over k external cohorts with \u226510 donors (n = patients);\n",
  "endothelial programs: single eligible cohort GSE132465 (k = 1, n = 11). Within-cohort: supported under the secondary within-cohort scaling.\n",
  "Spatial: pooled composition-adjusted coefficient (criterion A). Patients A / B: patients positive / evaluable for criteria A and B.\n",
  "* Not distinguishable from expression-matched random panels (Supplementary Table S15). Categories follow the pre-specified rule;\n",
  "spatial evidence cannot promote a program that failed single-cell validation. Individual genes are not separately validated markers."), theme = CAP(6.0))
save_fig(p4, "Figure4_evidence_scorecard", W2, 4.4)

# =========================================================================================
# FIGURE S1 - discovery donor variation
# Comments: explain target scaling and the evaluation / panel gene counts; point transparency; readable axes.
# =========================================================================================
tg_csv <- file.path("figure_data", "discovery_targets_per_donor.csv")
tg <- if (file.exists(tg_csv)) fread(tg_csv) else {
  suppressPackageStartupMessages({ library(Matrix); library(edgeR) })
  rbindlist(lapply(names(frz$programs), function(pg) { Pp <- frz$programs[[pg]]
    pbf <- if (Pp$compartment == "Lymphoid") "pseudobulk_Lymphoid_CD8.rds" else paste0("pseudobulk_", Pp$compartment, ".rds")
    pb <- readRDS(file.path("04_processed_data/GSE178341", pbf)); L <- edgeR::cpm(as.matrix(pb$counts[, Pp$donors]), log = TRUE, prior.count = 2)
    e <- Pp$evaluation_genes; y <- rowMeans(sweep(sweep(t(L[e$row_id, , drop = FALSE]), 2, e$eval_mean, "-"), 2, e$eval_sd, "/"))
    data.table(program = pg, compartment = Pp$compartment, donor = names(y), target = y) })) }
gc <- if (!is.null(frz)) data.table(program = names(frz$programs), n_eval = vapply(frz$programs, function(x) nrow(x$evaluation_genes), 1L),
                                    n_panel = vapply(frz$programs, function(x) nrow(x$panel), 1L)) else
      data.table(program = ev$program, n_eval = NA_integer_, n_panel = lengths(strsplit(ev$panel, ",\\s*")))
tg[, `:=`(plab = factor(PLAB[program], levels = rev(PLAB)), comp = COMP_OF[program])]
ns <- merge(tg[, .(n = .N), by = .(program, plab, comp)], gc, by = "program", all.x = TRUE)
ns[, lab := ifelse(is.na(n_eval), sprintf("n = %d donors; %d panel genes", n, n_panel),
                   sprintf("n = %d donors\n%d evaluation / %d panel genes", n, n_eval, n_panel))]
xmax1 <- max(tg$target, na.rm = TRUE)
ps1 <- ggplot(tg, aes(target, plab)) + geom_vline(xintercept = 0, colour = INK2, linewidth = 0.35, linetype = "22") +
  geom_point(aes(fill = comp), shape = 21, colour = "white", stroke = 0.15, size = 1.15, alpha = 0.5, position = position_jitter(height = 0.15, seed = 3)) +
  geom_boxplot(aes(colour = comp), fill = NA, outlier.shape = NA, width = 0.55, linewidth = 0.45) +
  geom_text(data = ns, aes(x = xmax1 + 0.1, label = lab), hjust = 0, size = 2.15, colour = INK, family = FONT, lineheight = 0.95) +
  scale_colour_manual(values = COMP_COL, guide = "none") + scale_fill_manual(values = COMP_COL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.5))) +
  labs(x = "Program target (mean z-score over evaluation genes; 0 = discovery average, units = discovery SD)", y = NULL,
       caption = paste0("Target: each evaluation gene's log2-CPM is centred and scaled with its discovery mean and SD, then averaged over the evaluation genes of the program.\n",
                        "Evaluation genes define the target; panel genes are the frozen predictors and are disjoint from them. One point per eligible discovery donor (GSE178341).")) +
  theme_pub() + CAP(6.0) + theme(axis.text.y = element_text(colour = INK, size = 7))
save_fig(ps1, "FigureS1_program_variation", W2, 3.6)

# =========================================================================================
# FIGURE S2 - panel versus pre-specified alternatives (A) and selection stability (B)
# Comments: typography; expose comparison models and their definition; explain frequencies from
# 50 fits; flag low-stability genes.
# =========================================================================================
cvp <- rbind(rd("06_results/panel/cv_performance_discovery_v2.csv"), rd("06_results/panel/cv_performance_discovery_lymphoid.csv"), fill = TRUE)
MLAB <- c(panel_refit = "Compact panel (frozen)", equal_weight = "Equal-weight top-4", conventional = "Conventional markers",
          panel_shrunk = "Panel, shrunk coefficients", larger_model = "Unrestricted elastic net")
cvp[, base := MAE[method == "mean"], by = program]
cvp <- cvp[method %in% names(MLAB)][, `:=`(lo = MAE_rep_min / base, hi = MAE_rep_max / base,
          mlab = factor(MLAB[method], levels = rev(MLAB)), plab = factor(PLAB[program], levels = PLAB), comp = COMP_OF[program], hl = method == "panel_refit")]
s2a <- ggplot(cvp, aes(y = mlab)) + annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf, fill = BAND) +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.3, linetype = "22") +
  geom_linerange(aes(xmin = lo, xmax = hi, colour = ifelse(hl, comp, "grey")), linewidth = 0.5) +
  geom_point(aes(x = MAE_vs_mean, fill = ifelse(hl, comp, "grey"), size = hl), shape = 21, colour = "white", stroke = 0.25) +
  scale_colour_manual(values = c(COMP_COL, grey = INK3), guide = "none") + scale_fill_manual(values = c(COMP_COL, grey = INK3), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 2.3, `FALSE` = 1.6), guide = "none") +
  scale_x_continuous(breaks = scales::breaks_extended(3)) + facet_wrap(~plab, ncol = 4) +
  labs(x = "Cross-validated MAE ratio vs discovery-mean prediction (point: mean over 10 repeats; bar: range across repeats)", y = NULL, tag = "A") +
  theme_pub(7.3) + theme(axis.text.y = element_text(colour = INK, size = 6.8), panel.border = element_rect(fill = NA, colour = GRID, linewidth = 0.3),
                         panel.spacing = unit(6, "pt"), strip.clip = "off")
stb <- rbind(rd("06_results/panel/selection_stability_v2.csv"), rd("06_results/panel/selection_stability_lymphoid.csv"), fill = TRUE)
pan <- rbindlist(lapply(seq_len(nrow(ev)), function(i) data.table(program = ev$program[i], gene_symbol = trimws(strsplit(ev$panel[i], ",")[[1]]))))
pan <- merge(pan, stb[, .(program, gene_symbol, selection_freq)], by = c("program", "gene_symbol"), all.x = TRUE)
pan[, `:=`(plab = factor(PLAB[program], levels = PLAB), comp = COMP_OF[program])]
pad <- pan[, .(k = seq_len(max(0, 4 - .N))), by = .(program, plab, comp)][, gene_symbol := paste0(".pad", k)][, k := NULL]
pan <- rbind(pan, pad[, selection_freq := 0], fill = TRUE)
pan[, is_pad := grepl("^\\.pad", gene_symbol)]; setorder(pan, program, -is_pad, selection_freq)
pan[, glab := factor(paste(program, gene_symbol), levels = unique(paste(program, gene_symbol)))]
pan[, low := !is_pad & selection_freq < 0.75]
s2b <- ggplot(pan, aes(selection_freq, glab, fill = comp)) +
  geom_col(aes(alpha = low), width = 0.66, na.rm = TRUE) + geom_vline(xintercept = 0.75, colour = INK2, linewidth = 0.35, linetype = "22") +
  geom_text(data = pan[is_pad == FALSE], aes(label = ifelse(low, sprintf("%.2f low", selection_freq), sprintf("%.2f", selection_freq)),
            fontface = ifelse(low, "bold", "plain")), hjust = -0.12, size = 2.1, colour = INK, family = FONT) +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.35), guide = "none") +
  scale_y_discrete(labels = function(z) ifelse(grepl("\\.pad", z), "", sub("^\\S+ ", "", z))) +
  scale_fill_manual(values = COMP_COL, guide = "none") + scale_x_continuous(limits = c(0, 1.42), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  facet_wrap(~plab, scales = "free_y", ncol = 4) +
  labs(x = "Selection frequency across 50 fits (5-fold cross-validation \u00d7 10 repeats); dashed line 0.75", y = NULL, tag = "B") +
  theme_pub(7.3) + theme(axis.text.y = element_text(face = "italic", colour = INK, size = 6.8), panel.border = element_rect(fill = NA, colour = GRID, linewidth = 0.3),
                         panel.spacing = unit(6, "pt"), strip.clip = "off")
cap_s2 <- paste0("A: all models fitted and compared within the same donor-grouped cross-validation folds; MAE ratio <1 (shaded) = better than predicting the discovery mean.\n",
                 "Compact panel (frozen) = relaxed elastic net, 1", EN, "4 genes, unpenalised refit; Panel, shrunk coefficients = same genes with elastic-net coefficients;\n",
                 "Equal-weight top-4 = equal-weight score of the 4 predictor genes most correlated with the target in training; Conventional markers = pre-specified literature\n",
                 "markers in the predictor pool (both calibrated in training); Unrestricted elastic net = elastic net at lambda.min without the gene cap (cost of compression).\n",
                 "B: share of the 50 fits in which each frozen panel gene was selected. Faded bars marked 'low': selected in <75% of fits (gene-level choice unstable).")
save_fig((s2a / s2b + plot_layout(heights = c(1.1, 1))) + plot_annotation(caption = cap_s2, theme = CAP(6.0)), "FigureS2_panel_baselines", W2, 8.0)

# =========================================================================================
# FIGURE S3 - RCTD marker check (script 41)
# Comments: separate own-type (diagonal) from cross-type (off-diagonal) correlations; readable axis
# text; present as a circular consistency check, not a validation.
# =========================================================================================
norm <- function(z) gsub("[^a-z0-9]", "", tolower(z))
m1 <- rd("06_results/spatial_RCTD_marker_check.csv")
mk <- melt(m1, id.vars = "marker_set", variable.name = "column", value.name = "rho")
mk[, column := sub("^col[0-9]+_", "", column)]
mk[, diag := norm(marker_set) == norm(column)]
if (!any(mk$diag)) message("NOTE: no marker set matched an RCTD column name; diagonal not outlined in Figure S3A")
lv <- sort(unique(mk$column)); mk[, `:=`(column = factor(column, levels = lv), marker_set = factor(marker_set, levels = rev(intersect(lv, unique(marker_set))) ))]
if (anyNA(mk$marker_set)) mk[, marker_set := factor(as.character(marker_set))]
pretty_lab <- function(z) gsub("_", " ", z)
d_med <- median(mk[diag == TRUE, rho]); o_med <- median(mk[diag == FALSE, rho])
s3a <- ggplot(mk, aes(column, marker_set, fill = rho)) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_tile(data = mk[diag == TRUE], fill = NA, colour = INK, linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", rho), fontface = ifelse(diag, "bold", "plain")), size = 2.3, colour = INK) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "Spearman \u03c1") +
  scale_x_discrete(labels = pretty_lab) + scale_y_discrete(labels = pretty_lab) +
  labs(x = "RCTD cell-type weight", y = "Marker set (mean marker score)", tag = "A",
       title = sprintf("Valdeolivas cohort (12 sections): outlined = own type (median \u03c1 %.2f); other cells = cross-type (median \u03c1 %.2f)", d_med, o_med)) +
  theme_pub(7.3) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6.8), axis.text.y = element_text(size = 6.8),
                         panel.grid.major.x = element_blank(), axis.line.x = element_blank(), legend.position = "right")
m2 <- rd("06_results/spatial_support_extra/RCTD_marker_check_extra.csv")
ordA <- rev(levels(mk$marker_set)); ordB <- c(intersect(ordA, unique(m2$type)), setdiff(sort(unique(m2$type)), ordA))
m2[, type := factor(type, levels = rev(ordB))]   # same row order as panel A
s3b <- ggplot(m2, aes(section, type, fill = own_rho)) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", own_rho)), size = 2.2, colour = INK) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "Spearman \u03c1") +
  scale_y_discrete(labels = pretty_lab) +
  labs(x = "Section", y = "Cell type", tag = "B", title = "GSE280318 and GSE226997: own-type correlation per section (diagonal only)") +
  theme_pub(7.3) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5), axis.text.y = element_text(size = 6.8),
                         panel.grid.major.x = element_blank(), axis.line.x = element_blank(), legend.position = "right")
cap_s3 <- paste0("Circular consistency check, not an independent validation: marker sets and the RCTD reference come from the same single-cell data, so high own-type correlation is expected\n",
                 "if deconvolution behaved as intended. Off-diagonal values show where cell types co-occur or share markers in tissue (informative about spot mixing).")
save_fig((s3a / s3b + plot_layout(heights = c(1.15, 1))) + plot_annotation(caption = cap_s3, theme = CAP(6.0)), "FigureS3_RCTD_marker_check", W2, 8.6)

# =========================================================================================
# FIGURE S4 - spatial support by cohort group (script 41)
# Comments: platforms and numbers for all groups; label the null line; CIs are not tests of function.
# =========================================================================================
S <- copy(spAll)[is.finite(pooled_effect)]
grp_n <- function(coh) uniqueN(spA[cohort %in% coh, patient])
GL <- c("All spatial cohorts (main)" = sprintf("All cohorts: %d patients (Visium, Visium HD, FFPE Visium)", grp_n(SP_COH)),
        "Valdeolivas only" = sprintf("Valdeolivas: %d patients (Visium)", grp_n("Valdeolivas")),
        "Additional cohorts only" = sprintf("GSE280318 + GSE226997: %d patients (Visium HD, FFPE Visium)", grp_n(c("GSE280318", "GSE226997"))))
S <- S[analysis %in% names(GL)]
S[, `:=`(plab = factor(PLAB[program], levels = rev(PLAB)), analysis = factor(GL[analysis], levels = GL))]
s4 <- ggplot(S, aes(pooled_effect, plab, colour = analysis, shape = analysis)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.4, colour = INK2) +
  annotate("text", x = 0, y = length(PLAB) + 0.65, label = "0 = no association", hjust = -0.05, size = 2.2, colour = INK, family = FONT) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.7, position = position_dodge(width = 0.65)) +
  geom_point(size = 2.1, fill = "white", stroke = 0.7, position = position_dodge(width = 0.65)) +
  scale_colour_manual(values = c("#1d1d1b", "#1b4965", "#e07a5f"), name = NULL) +
  scale_shape_manual(values = c(23, 21, 22), name = NULL) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.0))) +
  guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
  labs(x = "Composition-adjusted panel coefficient (criterion A), pooled across patients, 95% CI", y = NULL,
       caption = paste0("Patient numbers are the maximum per group; fewer patients were evaluable for angiogenesis and IFN-\u03b3 (probe-panel coverage).\n",
                        "Intervals describe the panel-target association in tissue; they are not tests of angiogenesis or of any other function.")) +
  theme_pub(7.3) + CAP(6.0) + theme(axis.text.y = element_text(colour = INK, size = 7), legend.text = element_text(size = 7))
save_fig(s4, "FigureS4_spatial_by_cohort", W2, 4.8)

# =========================================================================================
# FIGURE S5 - spatial reference sensitivity
# Comments: highlight the apoptosis 9/13 -> 8/13 change and reference sampling as a caveat;
# enlarge the symbol legend; no small grey explanatory text.
# =========================================================================================
sens_f <- "06_results/sens_randomref/comparison_published_vs_randomref.csv"
if (file.exists(sens_f) || file.exists(file.path("figure_data", gsub("/", "__", sens_f)))) {
  sv <- rd(sens_f)
  s5 <- rbind(sv[, .(program, ref = "Published reference (primary)", est = pooled_effect_published, lo = lo_published, hi = hi_published)],
              sv[, .(program, ref = "Random-selection reference", est = pooled_effect_randomref, lo = lo_randomref, hi = hi_randomref)])
  s5[, `:=`(plab = factor(PLAB[program], levels = rev(PLAB)), ref = factor(ref, levels = c("Published reference (primary)", "Random-selection reference")))]
  chg <- sv[label_changed == TRUE, .(program, plab = factor(PLAB[program], levels = rev(PLAB)),
           note = sprintf("Label changed: %s \u2192\n%s (criterion B\n%s \u2192 %s patients; threshold 9)", spatial_label_published, spatial_label_randomref,
                          B_positive_published, B_positive_randomref))]
  xr5 <- range(c(s5$lo, s5$hi, 0), na.rm = TRUE)
  ps5 <- ggplot(s5, aes(y = plab)) +
    geom_tile(data = chg, aes(x = mean(xr5), y = plab, width = Inf), height = 0.9, fill = HIL, inherit.aes = FALSE) +
    geom_vline(xintercept = 0, colour = INK2, linewidth = 0.35, linetype = "22") +
    geom_linerange(aes(xmin = lo, xmax = hi, colour = ref), linewidth = 0.8, position = position_dodge(width = 0.55)) +
    geom_point(aes(x = est, shape = ref, colour = ref), fill = "white", size = 2.4, stroke = 0.8, position = position_dodge(width = 0.55)) +
    geom_text(data = chg, aes(x = xr5[2], y = plab, label = note), hjust = 0, nudge_x = 0.01 * diff(xr5), size = 2.2, colour = INK,
              fontface = "bold", family = FONT, lineheight = 0.95, inherit.aes = FALSE) +
    scale_colour_manual(values = c("Published reference (primary)" = INK, "Random-selection reference" = "#c0441a"), name = NULL) +
    scale_shape_manual(values = c("Published reference (primary)" = 21, "Random-selection reference" = 23), name = NULL) +
    scale_y_discrete(limits = unname(rev(PLAB)[rev(PLAB) %in% as.character(s5$plab)])) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.7))) + coord_cartesian(clip = "off") +
    guides(colour = guide_legend(override.aes = list(size = 3)), shape = guide_legend(override.aes = list(size = 3))) +
    labs(x = "Pooled composition-adjusted panel coefficient (criterion A, 95% CI)", y = NULL,
         caption = paste0("Caveat: spatial results depend on how reference cells are sampled for RCTD. The random-selection reference draws cells in random order instead of\n",
                          "source-file order; all other steps unchanged. ",
                          if (nrow(chg)) "The only label change is highlighted; no evidence category changed (apoptosis is exploratory on its single-cell result)." else "No spatial label changed.")) +
    theme_pub(7.3) + CAP(6.2) + theme(axis.text.y = element_text(colour = INK, size = 7), legend.text = element_text(size = 7.3), legend.key.size = unit(11, "pt"))
  save_fig(ps5, "FigureS5_reference_sensitivity", W2, 4.0)
} else message("NOTE: ", sens_f, " not found; Figure S5 skipped")

writeLines(capture.output(sessionInfo()), file.path(OUT, paste0("sessionInfo_45", SUFFIX, ".txt")))
cat("\nRevised figures written to", OUT, "\n"); print(list.files(OUT, pattern = "\\.pdf$"))
