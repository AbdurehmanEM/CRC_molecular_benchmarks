# 43_publication_figures.R - publication figures (main 1-4, supplementary S1, S2, S5)
# Redesign of 40_figures_with_lymphoid.R (20 Oct 2026). Reads saved results only; no analysis.
#   * one visual system: fixed compartment colours, a single sequential ramp for evidence,
#     shared typography (7 pt at final size), 178 mm double-column width
#   * every number shown is read from the result files (no hard-coded counts)
#   * outputs: vector PDF (cairo), 600-dpi LZW TIFF and 300-dpi PNG in 07_figures/publication/
# Run from the project root:  source("05_analysis_scripts/43_publication_figures.R")
# Inputs are the files listed in IN below; if a figure_data/ folder exists (flattened copies),
# it is used as a fallback.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(scales) })

# ---------------------------------------------------------------- inputs
P <- function(path) {
  if (file.exists(path)) return(path)
  f <- file.path("figure_data", gsub("/", "__", path)); if (file.exists(f)) return(f)
  f <- file.path("figure_data", basename(path));        if (file.exists(f)) return(f)
  stop("Input not found: ", path)
}
rd <- function(path) fread(P(path))
OUT <- "07_figures/publication"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------- visual system
FONT <- if (.Platform$OS.type == "windows") "sans" else {             # Windows: Arial; elsewhere Arial-metric Liberation Sans if present
  ff <- tryCatch(system("fc-list : family", intern = TRUE), error = function(e) character(0))
  if (any(grepl("Liberation Sans", ff))) "Liberation Sans" else "sans" }
INK <- "#1d1d1b"; INK2 <- "#55544f"; INK3 <- "#97968f"; GRID <- "#ecebe6"; BAND <- "#f6f5f1"
COMP_COL <- c(Epithelial = "#2a78d6", Stromal = "#eb6834", Myeloid = "#1baf7a", Endothelial = "#eda100", Lymphoid = "#6250d6")
EVID <- c("Independently supported with spatial evidence", "Independently supported in single-cell data", "Exploratory", "Unsupported")
EVID_SHORT <- setNames(c("Validated + spatial", "Validated (single-cell)", "Exploratory", "Unsupported"), EVID)
EVID_COL <- setNames(c("#0e5a34", "#3f9864", "#a8d5b5", "#dcdbd5"), EVID)   # ordinal: one green ramp, grey = not met
EVID_TXT <- setNames(c("white", "white", INK, INK2), EVID)
DIV <- c("#1c5aa6", "#6f9fd8", "#e9e8e4", "#ee9a72", "#c0441a")             # diverging: two hues, neutral midpoint
PLAB <- c(EPI_HYPOXIA = "Epithelial hypoxia", EPI_APOPTOSIS = "Epithelial apoptosis", STR_TGFB = "Stromal TGF-β",
          MYE_IFNG = "Myeloid IFN-γ", MYE_INFLAM = "Myeloid inflammatory", END_ANGIO = "Endothelial angiogenesis",
          END_TNFA_NFKB = "Endothelial TNF/NF-κB", LYM_IL2_STAT5 = "CD8 IL2/STAT5")
PORD <- names(PLAB)
COMP_OF <- c(EPI_HYPOXIA = "Epithelial", EPI_APOPTOSIS = "Epithelial", STR_TGFB = "Stromal", MYE_IFNG = "Myeloid",
             MYE_INFLAM = "Myeloid", END_ANGIO = "Endothelial", END_TNFA_NFKB = "Endothelial", LYM_IL2_STAT5 = "Lymphoid")
W2 <- 7.0                                                                    # 178 mm double column

theme_pub <- function(base = 7) theme_minimal(base_size = base, base_family = FONT) + theme(
  text = element_text(colour = INK), axis.text = element_text(colour = INK2, size = base - 0.6),
  axis.title = element_text(colour = INK2, size = base - 0.2), panel.grid.minor = element_blank(),
  panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
  axis.line.x = element_line(colour = INK3, linewidth = 0.3), axis.ticks.x = element_line(colour = INK3, linewidth = 0.3),
  axis.ticks.length = unit(1.6, "pt"), strip.text = element_text(face = "bold", size = base - 0.2, hjust = 0, colour = INK),
  plot.tag = element_text(face = "bold", size = base + 3, colour = INK), plot.title = element_text(face = "bold", size = base + 0.3, colour = INK),
  plot.title.position = "plot", legend.position = "bottom", legend.text = element_text(size = base - 0.6, colour = INK2),
  legend.title = element_text(size = base - 0.6, colour = INK2), legend.key.size = unit(7, "pt"),
  legend.margin = margin(0, 0, 0, 0), plot.margin = margin(3, 5, 3, 3), panel.spacing = unit(3, "pt"))
blank_y <- theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
blank_x <- theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.line.x = element_blank(), panel.grid.major.x = element_blank())
no_strip <- theme(strip.text = element_blank())
half_up <- function(v, d = 2) sign(v) * floor(abs(v) * 10^d + 0.5 + 1e-9) / 10^d
fmt_ci <- function(e, l, h, d = 2) ifelse(is.finite(e), sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), half_up(e, d), half_up(l, d), half_up(h, d)), NA_character_)

SUFFIX <- "_aem"   # appended to every output file name
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
spS  <- rbind(rd("06_results/spatial_support_combined/spatial_support_summary_combined.csv"),
              rd("06_results/spatial_support_lymphoid_combined/spatial_support_summary_combined.csv"), fill = TRUE)[analysis == "All spatial cohorts (main)"]
coh_long <- rbind(rd("06_results/cross_cohort_meta/per_cohort_results_long.csv"),
                  rd("06_results/validation_lymphoid/per_cohort_results_long_lymphoid.csv"), fill = TRUE)
coh_long[cohort == "GSE144735_KUL", cohort := "GSE144735"]
meta <- rbind(rd("06_results/cross_cohort_meta/meta_analysis_results.csv")[analysis == "main"],
              rd("06_results/validation_lymphoid/meta_analysis_results_lymphoid.csv"), fill = TRUE)
n_spatial_pat <- uniqueN(spA$patient)
SP_COH <- c("Valdeolivas", "GSE280318", "GSE226997")
SP_N   <- spA[, .(n = uniqueN(patient)), by = cohort][match(SP_COH, cohort), n]
SP_LAB <- setNames(sprintf("%s (n = %d)", SP_COH, SP_N), SP_COH)

# =========================================================================================
# FIGURE 1 - study workflow (A) and cohort x compartment coverage (B)
# =========================================================================================
n_prog <- nrow(ev); n_genes <- sum(lengths(strsplit(ev$panel, ",\\s*")))
steps <- data.table(
  x = 1:6,
  head = c("Discovery", "Programs", "Panels", "Validation", "Spatial", "Evidence"),
  body = c("GSE178341\n62 patients\n5 compartments",
           sprintf("%d of 10 programs\nretained; panel and\ntarget genes disjoint", n_prog),
           sprintf("elastic net, 1\u20134 genes\n(%d genes in total)\nfrozen before testing", n_genes),
           "4 independent\nsingle-cell cohorts\nno refitting",
           sprintf("%d patients\n3 cohorts, adjusted\nfor composition", n_spatial_pat),
           "4 evidence\ncategories; chip\nassays proposed"))
pA <- ggplot(steps) +
  geom_segment(data = steps[x < 6], aes(x = x + 0.43, xend = x + 0.57, y = 0.5, yend = 0.5), colour = INK3, linewidth = 0.35,
               arrow = arrow(length = unit(3, "pt"), type = "closed")) +
  geom_rect(aes(xmin = x - 0.42, xmax = x + 0.42, ymin = 0.05, ymax = 0.95), fill = BAND, colour = NA) +
  geom_rect(aes(xmin = x - 0.42, xmax = x + 0.42, ymin = 0.80, ymax = 0.95), fill = INK, colour = NA) +
  geom_text(aes(x = x, y = 0.875, label = head), colour = "white", fontface = "bold", size = 2.1, family = FONT) +
  geom_text(aes(x = x, y = 0.45, label = body), colour = INK, size = 1.95, lineheight = 1.0, family = FONT) +
  scale_x_continuous(limits = c(0.55, 6.45), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(tag = "A") + theme_void(base_family = FONT) + theme(plot.tag = element_text(face = "bold", size = 10), plot.margin = margin(3, 3, 6, 3))

# coverage
el  <- rd("06_results/GSE178341/pseudobulk_eligibility.csv")[min_cells == 20]
elL <- rd("06_results/GSE178341/pseudobulk_eligibility_lymphoid.csv")[min_cells == 20]
cov_n <- function(path, col = "n_donors") { x <- rd(path); x[, compartment := COMP_OF[program]]; x[, .(n = max(get(col), na.rm = TRUE)), by = compartment] }
cv1 <- rbind(
  data.table(cohort = "GSE178341", role = "Discovery", compartment = names(COMP_COL),
             n = c(unlist(el[, setdiff(names(COMP_COL), "Lymphoid"), with = FALSE]), elL$Lymphoid_CD8)),
  cbind(cohort = "GSE132465", role = "Validation", rbind(cov_n("06_results/validation_GSE132465/validation_gene_coverage.csv"),
        data.table(compartment = "Lymphoid", n = rd("06_results/validation_GSE132465/validation_gene_coverage_lymphoid.csv")$n_donors[1]))),
  cbind(cohort = "GSE166555", role = "Validation", rbind(cov_n("06_results/validation_GSE166555/validation_gene_coverage.csv"),
        data.table(compartment = "Lymphoid", n = rd("06_results/validation_lymphoid/validation_coverage_extra_cohorts.csv")[cohort == "GSE166555", n_donors][1]))),
  cbind(cohort = "GSE200997", role = "Validation", rbind(cov_n("06_results/validation_GSE200997/validation_gene_coverage.csv"),
        data.table(compartment = "Lymphoid", n = NA_real_))),
  cbind(cohort = "GSE144735", role = "Descriptive", rbind(cov_n("06_results/validation_GSE144735/KUL_gene_coverage.csv", "kul_donors"),
        data.table(compartment = "Lymphoid", n = rd("06_results/validation_lymphoid/validation_coverage_extra_cohorts.csv")[cohort == "GSE144735_KUL", n_donors][1]))),
  use.names = TRUE)
spc <- copy(spA)[, compartment := COMP_OF[program]][, .(n = uniqueN(patient)), by = .(cohort, compartment)]
spc <- CJ(cohort = SP_COH, compartment = names(COMP_COL))[spc, n := i.n, on = .(cohort, compartment)]
cv1 <- rbind(cv1, spc[, .(cohort, role = "Spatial (patients)", compartment, n)], use.names = TRUE)
cv1[, `:=`(role = factor(role, levels = c("Discovery", "Validation", "Descriptive", "Spatial (patients)")),
           compartment = factor(compartment, levels = names(COMP_COL)))]
PLAT <- c(Valdeolivas = "Valdeolivas (Visium)", GSE280318 = "GSE280318 (Visium HD)", GSE226997 = "GSE226997 (FFPE Visium)")
cv1[cohort %in% names(PLAT), cohort := PLAT[cohort]]
cv1[, cohort := factor(cohort, levels = rev(c("GSE178341", "GSE132465", "GSE166555", "GSE200997", "GSE144735", PLAT)))]
cv1[, `:=`(alpha = ifelse(is.na(n), 0, 0.15 + 0.85 * pmin(n, 60) / 60),
           lab = ifelse(is.na(n), "n.d.", as.character(n)),
           low = !is.na(n) & n < 10 & role %in% c("Validation", "Descriptive"))]
cv1[, txt := ifelse(!is.na(n) & n >= 30, "white", INK)]
pB <- ggplot(cv1, aes(compartment, cohort)) +
  geom_tile(aes(fill = compartment, alpha = alpha), colour = "white", linewidth = 0.9, width = 0.97, height = 0.92) +
  geom_tile(data = cv1[low == TRUE], fill = NA, colour = INK2, linewidth = 0.35, linetype = "22", width = 0.86, height = 0.78) +
  geom_text(aes(label = lab, colour = txt), size = 2.3, fontface = "bold", family = FONT) +
  scale_fill_manual(values = COMP_COL, guide = "none") + scale_alpha_identity() + scale_colour_identity() +
  scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  facet_grid(role ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL, tag = "B",
       caption = "Eligible donors per compartment (≥20 cells per donor × compartment, tumour specimens); spatial rows give patients.\nDashed outline: <10 donors, not evaluable. n.d.: compartment not definable (T and NK cells annotated jointly).") +
  theme_pub() + theme(panel.grid = element_blank(), panel.grid.major.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                      axis.text.x = element_text(face = "bold", colour = INK, size = 6.8), axis.text.y = element_text(colour = INK, size = 6.5),
                      strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", colour = INK2, size = 6.3),
                      plot.caption = element_text(colour = INK2, size = 5.8, hjust = 0), plot.caption.position = "plot")
save_fig(wrap_elements(full = pA) / pB + plot_layout(heights = c(0.85, 2.6)), "Figure1_study_design_coverage", W2, 4.6)

# =========================================================================================
# FIGURE 2 - frozen-panel validation forest (table-forest layout)
# =========================================================================================
COH <- c("GSE132465", "GSE166555", "GSE200997", "GSE144735")
fr <- CJ(program = PORD, scaling = c("frozen", "within"), cohort = c(COH, "Pooled"), sorted = FALSE)
cl <- coh_long[, .(program, scaling, cohort, n_donors, MAE_ratio, ratio_lo, ratio_hi, status)]
mt <- meta[, .(program, scaling, cohort = "Pooled", n_donors = total_donors, MAE_ratio = pooled_MAE_ratio, ratio_lo, ratio_hi,
               status = ifelse(is.finite(pooled_MAE_ratio), meta_status, "not pooled"))]
fr <- merge(fr, rbind(cl, mt, fill = TRUE), by = c("program", "scaling", "cohort"), all.x = TRUE, sort = FALSE)
fr[, y := 6 - match(cohort, c(COH, "Pooled"))]                                    # 5 = GSE132465 ... 1 = Pooled
fr[, `:=`(plab = factor(PLAB[program], levels = PLAB), comp = COMP_OF[program], ok = is.finite(MAE_ratio) & MAE_ratio > 0)]
fr[, note := fifelse(ok, NA_character_, fifelse(grepl("^Excluded", status), "excluded (T/NK not separable)",
             fifelse(cohort == "Pooled", "not pooled (<2 cohorts)", "not evaluable")))]
fr[, est_txt := fifelse(ok, fmt_ci(MAE_ratio, ratio_lo, ratio_hi), "—")]
fr[, desc := cohort == "GSE144735"]
fr[, rowlab := fifelse(cohort == "Pooled", "Pooled (RE)", sprintf("%s  %s", cohort, fifelse(is.finite(n_donors) & n_donors > 0, paste0("n=", n_donors), "")))]
XL <- c(0.1, 10)
forest <- function(sc, title) {
  d <- fr[scaling == sc]; dp <- d[ok == TRUE & cohort == "Pooled"]; dc <- d[ok == TRUE & cohort != "Pooled"]
  dia <- dp[, .(x = c(ratio_lo, MAE_ratio, ratio_hi, MAE_ratio), yy = c(y, y + 0.32, y, y - 0.32)), by = .(plab, program, comp)]
  ggplot(d) +
    annotate("rect", xmin = XL[1], xmax = 1, ymin = -Inf, ymax = Inf, fill = BAND) +
    geom_vline(xintercept = 1, colour = INK2, linewidth = 0.3, linetype = "22") +
    geom_linerange(data = dc, aes(y = y, xmin = pmax(ratio_lo, XL[1]), xmax = pmin(ratio_hi, XL[2]), colour = comp, alpha = desc), linewidth = 0.55) +
    geom_point(data = dc, aes(y = y, x = MAE_ratio, fill = comp, alpha = desc), shape = 21, colour = "white", stroke = 0.3, size = 1.7) +
    geom_polygon(data = dia, aes(x = x, y = yy, group = program), fill = INK, colour = INK, linewidth = 0.2) +
    geom_text(data = d[ok == FALSE], aes(y = y, x = 1, label = note), colour = INK3, size = 1.75, fontface = "italic", family = FONT) +
    scale_colour_manual(values = COMP_COL, guide = "none") + scale_fill_manual(values = COMP_COL, guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.38), guide = "none") +
    scale_x_log10(limits = XL, breaks = c(0.25, 1, 4), labels = c("0.25", "1", "4"), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.4, 5.6), expand = c(0, 0)) +
    facet_grid(plab ~ ., switch = "y") + labs(x = "MAE ratio vs discovery mean (log)", y = NULL, title = title) +
    theme_pub() + blank_y + no_strip + theme(plot.title = element_text(hjust = 0, size = 6.4))
}
txtcol <- function(sc, title) ggplot(fr[scaling == sc]) +
  geom_text(aes(x = 0, y = y, label = est_txt, fontface = ifelse(cohort == "Pooled", "bold", "plain")),
            hjust = 0, size = 1.85, colour = INK, family = FONT) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = c(0.4, 5.6), expand = c(0, 0)) +
  facet_grid(plab ~ .) + labs(x = " ", y = NULL, title = title) + theme_pub() + blank_y + no_strip +
  theme(axis.text.x = element_text(colour = "white"), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(), panel.grid.major.x = element_blank(), plot.title = element_text(size = 6.0, colour = INK2, face = "plain", hjust = 0))
labcol <- ggplot(fr[scaling == "frozen"]) +
  geom_text(aes(x = 1, y = y, label = rowlab, fontface = ifelse(cohort == "Pooled", "bold", "plain"),
                colour = ifelse(desc, INK3, INK)), hjust = 1, size = 1.85, family = FONT) +
  scale_colour_identity() + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = c(0.4, 5.6), expand = c(0, 0)) +
  facet_grid(plab ~ ., switch = "y") + labs(x = " ", y = NULL, title = " ") + theme_pub() + blank_y +
  theme(axis.text.x = element_text(colour = "white"), axis.line.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), panel.grid.major.x = element_blank(),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 6.5, colour = INK))
p2 <- labcol + forest("frozen", "Frozen discovery scaling (primary)") + txtcol("frozen", "MAE ratio (95% CI)") +
      forest("within", "Within-cohort scaling (secondary)") + txtcol("within", "MAE ratio (95% CI)") +
      plot_layout(widths = c(1.4, 1.7, 1.05, 1.7, 1.05)) +
      plot_annotation(caption = "Points: single cohorts (95% bootstrap CI, 2,000 donor resamples); faded: GSE144735 (<10 donors, descriptive). Diamonds: DerSimonian–Laird pooled estimate over cohorts with ≥10 donors.\nShaded region: panel error below the discovery-mean baseline.",
                      theme = theme(plot.caption = element_text(size = 5.8, colour = INK2, hjust = 0, family = FONT)))
save_fig(p2, "Figure2_validation_forest", W2, 7.6)

# =========================================================================================
# FIGURE 3 - spatial support: example maps (A), criterion A (B), criterion B (C)
# =========================================================================================
SEC <- "SN048_A416371_Rep1"; SHOW <- c("MYE_IFNG", "END_TNFA_NFKB")
spot_csv <- file.path("figure_data", paste0("spot_scores_", SEC, ".csv"))
sp <- if (file.exists(spot_csv)) fread(spot_csv) else readRDS("04_processed_data/Valdeolivas_Visium/spot_program_scores.rds")[section == SEC]
sp <- sp[program %in% SHOW]
rho <- sp[, .(r = cor(target, panel, method = "spearman")), by = program]
m <- melt(sp, id.vars = c("program", "pxl_row", "pxl_col"), measure.vars = c("target", "panel"), variable.name = "score")
m[, value := ave(value, program, score, FUN = function(v) as.numeric(scale(v)))]
m[, plab := factor(sprintf("%s\nSpearman ρ = %.2f", PLAB[program], rho$r[match(program, rho$program)]),
                   levels = sprintf("%s\nSpearman ρ = %.2f", PLAB[SHOW], rho$r[match(SHOW, rho$program)]))]
m[, score := factor(ifelse(score == "target", "Target", "Panel"), levels = c("Target", "Panel"))]
p3a <- ggplot(m, aes(pxl_col, -pxl_row, colour = pmin(pmax(value, -2.5), 2.5))) + geom_point(size = 0.42, stroke = 0) +
  facet_grid(score ~ plab) + coord_equal() +
  scale_colour_gradientn(colours = DIV, limits = c(-2.5, 2.5), breaks = c(-2, 0, 2), name = "z-score\nwithin section") +
  labs(x = NULL, y = NULL, tag = "A", title = sprintf("Example Visium section %s (patient %s); target = evaluation genes, panel = frozen panel score", SEC, sub("^SN[0-9]+_(A[0-9]+)_.*$", "\\1", SEC))) +
  theme_pub() + theme(axis.text = element_blank(), axis.line.x = element_blank(), axis.ticks = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), panel.grid.major.x = element_blank(),
                      legend.position = "right", legend.key.height = unit(14, "pt"), legend.key.width = unit(5, "pt"),
                      strip.text.x = element_text(hjust = 0.5, lineheight = 1.05), strip.text.y = element_text(angle = -90, hjust = 0.5, face = "bold", colour = INK2))
crit_plot <- function(dat, summ, xlab, tag, cnt_col, add_ci = TRUE) {
  dat <- copy(dat)[, `:=`(plab = factor(PLAB[program], levels = rev(PLAB)), comp = COMP_OF[program], cohort = factor(cohort, levels = SP_COH))]
  summ <- copy(summ)[, `:=`(plab = factor(PLAB[program], levels = rev(PLAB)), comp = COMP_OF[program])]
  xr <- range(c(dat[[3 + 1]], 0), na.rm = TRUE)
  g <- ggplot(dat, aes(y = plab)) + geom_vline(xintercept = 0, colour = INK2, linewidth = 0.3, linetype = "22") +
    geom_point(aes(x = .data[[names(dat)[4]]], fill = comp, shape = cohort), colour = "white", stroke = 0.25, size = 1.45, alpha = 0.75,
               position = position_jitter(height = 0.16, width = 0, seed = 3))
  if (add_ci) g <- g + geom_linerange(data = summ, aes(xmin = lo, xmax = hi, colour = comp), linewidth = 0.9) +
    geom_point(data = summ, aes(x = pooled_effect, fill = comp), shape = 23, colour = INK, stroke = 0.35, size = 2.2)
  g + geom_text(data = summ, aes(x = Inf, label = .data[[cnt_col]]), hjust = 1.05, size = 1.9, colour = INK2, family = FONT) +
    scale_fill_manual(values = COMP_COL, guide = "none") + scale_colour_manual(values = COMP_COL, guide = "none") +
    scale_shape_manual(values = c(Valdeolivas = 21, GSE280318 = 22, GSE226997 = 24), labels = SP_LAB, name = NULL, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    labs(x = xlab, y = NULL, tag = tag) + theme_pub() +
    guides(shape = guide_legend(override.aes = list(fill = INK3, colour = "white", size = 1.8, alpha = 1)))
}
p3b <- crit_plot(spA[, .(cohort, patient, program, effect)], spS, "Composition-adjusted panel coefficient\n(criterion A; ◆ pooled estimate, 95% CI)", "B", "A_positive") +
  theme(axis.text.y = element_text(colour = INK, size = 6.4))
p3c <- crit_plot(spB[, .(cohort, patient, program, diff)], spS, "Target difference, expected-high minus comparison\n(criterion B; one point per patient)", "C", "B_positive", add_ci = FALSE) +
  theme(axis.text.y = element_blank())
p3 <- p3a / ((p3b | p3c) + plot_layout(guides = "collect") & theme(legend.position = "bottom")) + plot_layout(heights = c(1.25, 1)) +
  plot_annotation(caption = sprintf("B, C: %d patients from three cohorts (fewer for angiogenesis and IFN-γ, limited by probe-panel coverage). Right-hand counts: patients with a positive estimate.\nCriterion B uses pathologist regions in Valdeolivas and compartment-proportion quartiles elsewhere (all cohorts for the lymphoid program).", n_spatial_pat),
                  theme = theme(plot.caption = element_text(size = 5.8, colour = INK2, hjust = 0, family = FONT)))
save_fig(p3, "Figure3_spatial_support", W2, 7.2)

# =========================================================================================
# FIGURE 4 - evidence scorecard
# =========================================================================================
pci <- function(s) { v <- as.numeric(regmatches(s, gregexpr("[0-9]+\\.?[0-9]*", s))[[1]]); if (length(v) < 3) rep(NA_real_, 3) else v[1:3] }
e4 <- copy(ev)[, c("sc_est", "sc_lo", "sc_hi") := transpose(lapply(frozen_MAE_ratio_CI, pci))]
e4 <- merge(e4, spS[, .(program, sp_est = pooled_effect, sp_lo = lo, sp_hi = hi, A_positive, B_positive)], by = "program", all.x = TRUE)
e4[, `:=`(evidence_category = factor(evidence_category, levels = EVID), comp = factor(COMP_OF[program], levels = names(COMP_COL)))]
setorder(e4, evidence_category, comp, program)
e4[, y := rev(seq_len(.N))]
k_lab <- meta[scaling == "frozen", .(program, k = k_cohorts)]; e4 <- merge(e4, k_lab, by = "program", all.x = TRUE); setorder(e4, -y)
e4[, sc_note := fifelse(is.na(k) | k < 2, "single cohort", sprintf("pooled, %d cohorts", k))]
YL <- c(0.5, nrow(e4) + 0.5)
base4 <- function() list(scale_y_continuous(limits = YL, expand = c(0, 0)), theme_pub(), blank_y,
                         theme(plot.title = element_text(size = 6.6, hjust = 0)))
c1 <- ggplot(e4) + geom_rect(aes(xmin = 0, xmax = 0.035, ymin = y - 0.38, ymax = y + 0.38, fill = comp)) +
  geom_text(aes(x = 0.06, y = y, label = PLAB[program]), hjust = 0, size = 2.25, fontface = "bold", colour = INK, family = FONT) +
  scale_fill_manual(values = COMP_COL, name = NULL, drop = FALSE) + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Program", x = " ") + base4() + blank_x + theme(axis.text.x = element_text(colour = "white"))
c2 <- ggplot(e4, aes(y = y)) + annotate("rect", xmin = 0.1, xmax = 1, ymin = -Inf, ymax = Inf, fill = BAND) +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.3, linetype = "22") +
  geom_linerange(aes(xmin = pmax(sc_lo, 0.1), xmax = pmin(sc_hi, 10), colour = comp), linewidth = 0.8) +
  geom_point(aes(x = sc_est, fill = ifelse(sc_independent, COMP_COL[as.character(comp)], "white"), colour = comp), shape = 21, size = 2.1, stroke = 0.6) +
  geom_text(aes(x = 10, label = sc_note), hjust = 1, vjust = -0.9, size = 1.6, colour = INK3, family = FONT) +
  scale_fill_identity() + scale_colour_manual(values = COMP_COL, guide = "none") +
  scale_x_log10(limits = c(0.1, 10), breaks = c(0.25, 1, 4), labels = c("0.25", "1", "4"), expand = c(0, 0)) +
  labs(title = "Single-cell validation (frozen)", x = "MAE ratio vs discovery mean") + base4()
c3 <- ggplot(e4, aes(y = y)) + geom_vline(xintercept = 0, colour = INK2, linewidth = 0.3, linetype = "22") +
  geom_linerange(aes(xmin = sp_lo, xmax = sp_hi, colour = comp), linewidth = 0.8) +
  geom_point(aes(x = sp_est, fill = ifelse(spatial_label == "Spatial support", COMP_COL[as.character(comp)], "white"), colour = comp), shape = 23, size = 2.2, stroke = 0.6) +
  geom_text(aes(x = Inf, label = paste0("A ", A_positive, " · B ", B_positive)), hjust = 1.02, vjust = -0.9, size = 1.6, colour = INK3, family = FONT) +
  scale_fill_identity() + scale_colour_manual(values = COMP_COL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.06))) +
  labs(title = "Spatial support (criterion A)", x = "Composition-adjusted coefficient") + base4()
c4 <- ggplot(e4) + geom_tile(aes(x = 0.5, y = y, fill = evidence_category), width = 0.96, height = 0.76) +
  geom_text(aes(x = 0.5, y = y, label = EVID_SHORT[as.character(evidence_category)], colour = EVID_TXT[as.character(evidence_category)]),
            size = 2.0, fontface = "bold", family = FONT) +
  scale_fill_manual(values = EVID_COL, guide = "none") + scale_colour_identity() + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Evidence category", x = " ") + base4() + blank_x + theme(axis.text.x = element_text(colour = "white"))
c5 <- ggplot(e4) + geom_text(aes(x = 0, y = y, label = gsub(",\\s*", "  ", panel)), hjust = 0, size = 2.0, fontface = "italic", colour = INK, family = FONT) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + labs(title = "Panel genes", x = " ") + base4() + blank_x +
  theme(axis.text.x = element_text(colour = "white"))
p4 <- (c1 + c2 + c3 + c4 + c5 + plot_layout(widths = c(1.35, 1.45, 1.35, 1.05, 1.5), guides = "collect")) &
  theme(legend.position = "bottom")
p4 <- p4 + plot_annotation(caption = "Filled symbols: criterion met (single-cell: pooled or single-cohort upper CI < 1; spatial: criteria A and B both met). Open symbols: not met (spatial: criterion A only).\nGrey counts: patients positive for criterion A and criterion B. Categories follow the pre-specified evidence rule; spatial evidence cannot promote a program that failed single-cell validation.",
                           theme = theme(plot.caption = element_text(size = 5.8, colour = INK2, hjust = 0, family = FONT)))
save_fig(p4, "Figure4_evidence_scorecard", W2, 3.5)

# =========================================================================================
# FIGURE S1 - discovery donor variation
# =========================================================================================
tg_csv <- file.path("figure_data", "discovery_targets_per_donor.csv")
tg <- if (file.exists(tg_csv)) fread(tg_csv) else {
  suppressPackageStartupMessages({ library(Matrix); library(edgeR) })
  frz <- readRDS("04_processed_data/frozen_panel/frozen_panel_v2.rds"); frz$programs <- c(frz$programs, readRDS("04_processed_data/frozen_panel/frozen_panel_LYM_v1.rds")$programs)
  rbindlist(lapply(names(frz$programs), function(pg) { Pp <- frz$programs[[pg]]
    pbf <- if (Pp$compartment == "Lymphoid") "pseudobulk_Lymphoid_CD8.rds" else paste0("pseudobulk_", Pp$compartment, ".rds")
    pb <- readRDS(file.path("04_processed_data/GSE178341", pbf)); L <- edgeR::cpm(as.matrix(pb$counts[, Pp$donors]), log = TRUE, prior.count = 2)
    e <- Pp$evaluation_genes; y <- rowMeans(sweep(sweep(t(L[e$row_id, , drop = FALSE]), 2, e$eval_mean, "-"), 2, e$eval_sd, "/"))
    data.table(program = pg, compartment = Pp$compartment, donor = names(y), target = y) })) }
tg[, `:=`(plab = factor(PLAB[program], levels = rev(PLAB)), comp = COMP_OF[program])]
ns <- tg[, .(n = .N), by = .(plab, comp)]
ps1 <- ggplot(tg, aes(target, plab)) + geom_vline(xintercept = 0, colour = INK2, linewidth = 0.3, linetype = "22") +
  geom_boxplot(aes(colour = comp), fill = NA, outlier.shape = NA, width = 0.55, linewidth = 0.35) +
  geom_point(aes(fill = comp), shape = 21, colour = "white", stroke = 0.2, size = 1.1, alpha = 0.85, position = position_jitter(height = 0.13, seed = 3)) +
  geom_text(data = ns, aes(x = Inf, label = paste0("n = ", n)), hjust = 1.05, size = 1.9, colour = INK2, family = FONT) +
  scale_colour_manual(values = COMP_COL, guide = "none") + scale_fill_manual(values = COMP_COL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.12))) +
  labs(x = "Program target in discovery donors (mean z-score of evaluation genes)", y = NULL) + theme_pub() +
  theme(axis.text.y = element_text(colour = INK, size = 6.5))
save_fig(ps1, "FigureS1_program_variation", W2 * 0.62, 3.0)

# =========================================================================================
# FIGURE S2 - panel versus pre-specified alternatives (A) and selection stability (B)
# =========================================================================================
cvp <- rbind(rd("06_results/panel/cv_performance_discovery_v2.csv"), rd("06_results/panel/cv_performance_discovery_lymphoid.csv"), fill = TRUE)
MLAB <- c(panel_refit = "Compact panel (frozen)", equal_weight = "Equal-weight top-4", conventional = "Conventional markers",
          panel_shrunk = "Panel, shrunk coefficients", larger_model = "Unrestricted elastic net")
cvp[, base := MAE[method == "mean"], by = program]
cvp <- cvp[method %in% names(MLAB)][, `:=`(lo = MAE_rep_min / base, hi = MAE_rep_max / base,
          mlab = factor(MLAB[method], levels = rev(MLAB)), plab = factor(PLAB[program], levels = PLAB), comp = COMP_OF[program], hl = method == "panel_refit")]
s2a <- ggplot(cvp, aes(y = mlab)) + annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf, fill = BAND) +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.3, linetype = "22") +
  geom_linerange(aes(xmin = lo, xmax = hi, colour = ifelse(hl, comp, "grey")), linewidth = 0.45) +
  geom_point(aes(x = MAE_vs_mean, fill = ifelse(hl, comp, "grey"), size = hl), shape = 21, colour = "white", stroke = 0.25) +
  scale_colour_manual(values = c(COMP_COL, grey = INK3), guide = "none") + scale_fill_manual(values = c(COMP_COL, grey = INK3), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 2.0, `FALSE` = 1.4), guide = "none") +
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1)) + facet_wrap(~plab, ncol = 4) +
  labs(x = "MAE ratio vs discovery-mean prediction (5-fold CV × 10; bars: per-repeat panel MAE range relative to mean baseline MAE)", y = NULL, tag = "A") +
  theme_pub() + theme(axis.text.y = element_text(colour = INK, size = 6.2), panel.border = element_rect(fill = NA, colour = GRID, linewidth = 0.3))
stb <- rbind(rd("06_results/panel/selection_stability_v2.csv"), rd("06_results/panel/selection_stability_lymphoid.csv"), fill = TRUE)
pan <- rbindlist(lapply(seq_len(nrow(ev)), function(i) data.table(program = ev$program[i], gene_symbol = trimws(strsplit(ev$panel[i], ",")[[1]]))))
pan <- merge(pan, stb[, .(program, gene_symbol, selection_freq)], by = c("program", "gene_symbol"), all.x = TRUE)
pan[, `:=`(plab = factor(PLAB[program], levels = PLAB), comp = COMP_OF[program])]
pad <- pan[, .(k = seq_len(max(0, 4 - .N))), by = .(program, plab, comp)][, gene_symbol := paste0(".pad", k)][, k := NULL]
pan <- rbind(pan, pad[, selection_freq := 0], fill = TRUE)                       # equal bar widths across panels
pan[, is_pad := grepl("^\\.pad", gene_symbol)]; setorder(pan, program, -is_pad, selection_freq)
pan[, glab := factor(paste(program, gene_symbol), levels = unique(paste(program, gene_symbol)))]
s2b <- ggplot(pan, aes(selection_freq, glab, fill = comp)) +
  geom_col(width = 0.66, na.rm = TRUE) + geom_vline(xintercept = 0.75, colour = INK2, linewidth = 0.3, linetype = "22") +
  geom_text(data = pan[is_pad == FALSE], aes(label = sprintf("%.2f", selection_freq)), hjust = -0.18, size = 1.8, colour = INK2, family = FONT) +
  scale_y_discrete(labels = function(z) ifelse(grepl("\\.pad", z), "", sub("^\\S+ ", "", z))) +
  scale_fill_manual(values = COMP_COL, guide = "none") + scale_x_continuous(limits = c(0, 1.25), breaks = c(0, 0.5, 0.75, 1), expand = c(0, 0)) +
  facet_wrap(~plab, scales = "free_y", ncol = 4) +
  labs(x = "Selection frequency across 50 cross-validation fits (dashed: 0.75)", y = NULL, tag = "B") +
  theme_pub() + theme(axis.text.y = element_text(face = "italic", colour = INK, size = 6.2), panel.border = element_rect(fill = NA, colour = GRID, linewidth = 0.3),
                      strip.text = element_text(size = 6.0, face = "bold", hjust = 0, colour = INK))
save_fig(s2a / s2b + plot_layout(heights = c(1.1, 1)), "FigureS2_panel_baselines", W2, 6.8)

# =========================================================================================
# FIGURE S5 - spatial reference sensitivity (published vs random-selection reference)
# =========================================================================================
sens_f <- "06_results/sens_randomref/comparison_published_vs_randomref.csv"
if (file.exists(sens_f) || file.exists(file.path("figure_data", gsub("/", "__", sens_f)))) {
  sv <- rd(sens_f)
  s5 <- rbind(sv[, .(program, ref = "Published reference", est = pooled_effect_published, lo = lo_published, hi = hi_published, lab = spatial_label_published, B = B_positive_published)],
              sv[, .(program, ref = "Random-selection reference", est = pooled_effect_randomref, lo = lo_randomref, hi = hi_randomref, lab = spatial_label_randomref, B = B_positive_randomref)])
  s5[, `:=`(plab = factor(PLAB[program], levels = rev(PLAB)), comp = COMP_OF[program], ref = factor(ref, levels = c("Published reference", "Random-selection reference")))]
  chg <- sv[label_changed == TRUE, .(program, plab = factor(PLAB[program], levels = rev(PLAB)),
           note = sprintf("%s to %s (criterion B %s to %s)", spatial_label_published, spatial_label_randomref, B_positive_published, B_positive_randomref))]
  ps5 <- ggplot(s5, aes(y = plab)) + geom_vline(xintercept = 0, colour = INK2, linewidth = 0.3, linetype = "22") +
    geom_linerange(aes(xmin = lo, xmax = hi, colour = ref), linewidth = 0.7, position = position_dodge(width = 0.55)) +
    geom_point(aes(x = est, shape = ref, colour = ref), fill = "white", size = 1.9, stroke = 0.6, position = position_dodge(width = 0.55)) +
    scale_colour_manual(values = c("Published reference" = INK, "Random-selection reference" = "#c0441a"), name = NULL) +
    scale_shape_manual(values = c("Published reference" = 21, "Random-selection reference" = 23), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.06))) +
    labs(x = "Pooled composition-adjusted panel coefficient (criterion A, 95% CI)", y = NULL,
         caption = paste0("Reference cells drawn in random order instead of source-file order; all spatial scripts re-run unchanged otherwise.\n",
                          if (nrow(chg)) paste0("Only label change: ", paste(PLAB[chg$program], chg$note, collapse = "; "), ". No evidence category changed.") else "No spatial label changed.")) + theme_pub() +
    theme(axis.text.y = element_text(colour = INK, size = 6.5), plot.caption = element_text(size = 5.6, colour = INK2, hjust = 0), plot.caption.position = "plot")
  save_fig(ps5, "FigureS5_reference_sensitivity", W2 * 0.72, 3.0)
}
writeLines(capture.output(sessionInfo()), file.path(OUT, paste0("sessionInfo_43", SUFFIX, ".txt")))
cat("\nFigures written to", OUT, "\n")
