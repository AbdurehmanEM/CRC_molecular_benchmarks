
# 21_figures.R - publication Figures 1-6 from saved results (no new analysis)
# Double-column width 7.2 in (183 mm); vector PDF + 600-dpi TIFF (LZW) + 300-dpi PNG preview
for (p in c("ggplot2", "patchwork")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(ggplot2); library(patchwork); library(edgeR) })
frz <- readRDS("04_processed_data/frozen_panel/frozen_panel_v2.rds")
dir.create("07_figures", showWarnings = FALSE)
W <- 7.2; INK <- "#222220"; INK2 <- "#5f5e58"; GRID <- "#e8e7e1"
COMP_COL <- c(Epithelial = "#2a78d6", Stromal = "#eb6834", Myeloid = "#1baf7a", Endothelial = "#e0a100")
COH_COL <- c(GSE132465 = "#2a78d6", GSE144735 = "#eb6834", GSE166555 = "#1baf7a", GSE200997 = "#d6338a", Pooled = INK)
DIV <- c("#184f95", "#5598e7", "#cde2fb", "#f4f3ef", "#fbd3c4", "#e8693d", "#a8331a")
CAT <- c("Independently supported with spatial evidence", "Independently supported in single-cell data", "Exploratory", "Unsupported")
CATCOL <- setNames(c("#0f7a43", "#58b37a", "#f0b93b", "#d6d5cf"), CAT)
PLAB <- c(EPI_HYPOXIA = "Epithelial hypoxia", EPI_APOPTOSIS = "Epithelial apoptosis", STR_TGFB = "Stromal TGF-\u03b2",
  MYE_IFNG = "Myeloid IFN-\u03b3", MYE_INFLAM = "Myeloid inflammatory", END_ANGIO = "Endothelial angiogenesis",
  END_TNFA_NFKB = "Endothelial TNF/NF-\u03baB")
PORD <- names(PLAB); comp_of <- sapply(frz$programs, `[[`, "compartment")
tint <- function(col, f) { m <- col2rgb(col) / 255; rgb(t(1 - f * (1 - m))) }
th <- theme_classic(base_size = 8) + theme(text = element_text(colour = INK), axis.text = element_text(colour = INK2, size = 7),
  axis.line = element_line(colour = INK2, linewidth = 0.3), axis.ticks = element_line(colour = INK2, linewidth = 0.3),
  panel.grid.major.x = element_line(colour = GRID, linewidth = 0.25), strip.background = element_blank(),
  strip.text = element_text(face = "bold", hjust = 0, size = 7.5), plot.title = element_text(face = "bold", size = 8.5),
  plot.tag = element_text(face = "bold", size = 11), legend.position = "bottom", legend.text = element_text(size = 7),
  legend.key.size = unit(0.35, "cm"), strip.clip = "off", plot.caption = element_text(colour = INK2, size = 6.5, hjust = 0))
save_fig <- function(p, name, w, h) {
  ggsave(file.path("07_figures", paste0(name, "_updated.pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path("07_figures", paste0(name, "_updated.tiff")), p, width = w, height = h, dpi = 600, compression = "lzw", bg = "white")
  ggsave(file.path("07_figures", paste0(name, "_updated.png")), p, width = w, height = h, dpi = 300, bg = "white"); message("saved ", name, "_updated") }
el <- fread("06_results/GSE178341/pseudobulk_eligibility.csv")[min_cells == 20]
f1 <- data.table(cohort = "GSE178341 (discovery)", compartment = names(COMP_COL), n = unlist(el[, names(COMP_COL), with = FALSE]))
cov_n <- function(path, coh, col = "n_donors") { x <- fread(path); x[, compartment := comp_of[program]]
  x[, .(n = max(get(col), na.rm = TRUE)), by = compartment][, cohort := coh][] }
f1 <- rbind(f1, cov_n("06_results/validation_GSE132465/validation_gene_coverage.csv", "GSE132465 (validation)"),
  cov_n("06_results/validation_GSE144735/KUL_gene_coverage.csv", "GSE144735 (descriptive)", "kul_donors"),
  cov_n("06_results/validation_GSE166555/validation_gene_coverage.csv", "GSE166555 (validation)"),
  cov_n("06_results/validation_GSE200997/validation_gene_coverage.csv", "GSE200997 (validation; own annotation)"),
  data.table(cohort = "Valdeolivas Visium (spatial; patients)", compartment = names(COMP_COL), n = 6L), use.names = TRUE)
f1[, cohort := factor(cohort, levels = rev(unique(cohort)))]; f1[, compartment := factor(compartment, levels = names(COMP_COL))]
f1[, evaluable := n >= 10 | grepl("spatial", cohort)]
f1[, fillcol := mapply(function(cc, nn) tint(COMP_COL[[cc]], 0.18 + 0.82 * min(nn, 60) / 60), as.character(compartment), n)]
f1[, txtcol := ifelse(n >= 30, "white", INK)]
p1 <- ggplot(f1, aes(compartment, cohort)) +
  geom_tile(aes(fill = fillcol), colour = "white", linewidth = 1.2) + scale_fill_identity() +
  geom_tile(data = f1[evaluable == FALSE], fill = NA, colour = INK2, linewidth = 0.45, linetype = "22", width = 0.9, height = 0.84) +
  geom_text(aes(label = n, colour = txtcol), size = 3, fontface = "bold") + scale_colour_identity() +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, caption = "Numbers: eligible donors (\u226520 cells per donor \u00d7 compartment; tumour specimens);\ncolour intensity scales with donor number.\nDashed outline: <10 donors (not evaluable). Spatial row: patients with primary CRC (6 patients, 12 sections).") +
  th + theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", size = 8, colour = INK), axis.text.y = element_text(size = 7.5, colour = INK))
save_fig(p1, "Figure1_cohorts_coverage", W, 3.2)
f2 <- rbindlist(lapply(PORD, function(pg) { P <- frz$programs[[pg]]
  pb <- readRDS(file.path("04_processed_data/GSE178341", paste0("pseudobulk_", P$compartment, ".rds")))
  L <- edgeR::cpm(as.matrix(pb$counts[, P$donors]), log = TRUE, prior.count = 2); ev <- P$evaluation_genes
  y <- rowMeans(sweep(sweep(t(L[ev$row_id, , drop = FALSE]), 2, ev$eval_mean, "-"), 2, ev$eval_sd, "/"))
  data.table(program = pg, donor = names(y), target = y, compartment = P$compartment, n_eval = nrow(ev), n_panel = nrow(P$panel)) }))
f2[, `:=`(plab = factor(PLAB[program], levels = rev(PLAB)), compartment = factor(compartment, levels = names(COMP_COL)))]
lab2 <- f2[, .(n = .N, n_eval = n_eval[1], n_panel = n_panel[1]), by = .(plab, compartment)]
p2 <- ggplot(f2, aes(target, plab)) +
  geom_vline(xintercept = 0, colour = INK2, linewidth = 0.3, linetype = "22") +
  geom_boxplot(aes(fill = compartment), outlier.shape = NA, width = 0.6, colour = INK2, linewidth = 0.3, alpha = 0.25) +
  geom_point(aes(colour = compartment), position = position_jitter(height = 0.14, seed = 3), size = 0.9, alpha = 0.85) +
  geom_text(data = lab2, aes(x = max(f2$target) + 0.08, label = sprintf("n = %d donors\n%d evaluation / %d panel genes", n, n_eval, n_panel)),
    hjust = 0, lineheight = 0.9, size = 2.2, colour = INK2) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.32))) +
  scale_colour_manual(values = COMP_COL, breaks = names(COMP_COL), name = NULL) +
  scale_fill_manual(values = COMP_COL, breaks = names(COMP_COL), guide = "none") +
  labs(x = "Program target in discovery donors (mean z-score of evaluation genes)", y = NULL) + th
save_fig(p2, "Figure2_program_variation", W, 4)
cvp <- fread("06_results/panel/cv_performance_discovery_v2.csv")
MLAB <- c(panel_refit = "Compact panel (frozen)", panel_shrunk = "Panel, shrunk coefficients", equal_weight = "Equal-weight top-4",
  conventional = "Conventional markers", larger_model = "Larger elastic net", mean = "Discovery mean")
cvp[, `:=`(ratio_lo = MAE_rep_min / MAE[method == "mean"], ratio_hi = MAE_rep_max / MAE[method == "mean"]), by = program]
cvp <- cvp[method %in% names(MLAB) & method != "mean"]
cvp[, `:=`(mlab = factor(MLAB[method], levels = rev(MLAB)), plab = factor(PLAB[program], levels = PLAB),
  pcol = ifelse(method == "panel_refit", COMP_COL[comp_of[program]], "#9a998f"))]
p3a <- ggplot(cvp, aes(MAE_vs_mean, mlab)) +
  annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf, fill = "#eef5ee", alpha = 0.7) +
  geom_vline(xintercept = 1, colour = INK2, linetype = "22", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ratio_lo, xmax = ratio_hi, colour = pcol), height = 0, linewidth = 0.5) +
  geom_point(aes(colour = pcol, size = method == "panel_refit")) + scale_colour_identity() +
  scale_size_manual(values = c(`TRUE` = 2.3, `FALSE` = 1.5), guide = "none") +
  facet_wrap(~plab, ncol = 4) + scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1)) +
  labs(x = "MAE ratio vs discovery-mean prediction \n(5-fold CV \u00d7 10; bars = range across repeats; shaded = better than baseline)",
    y = NULL, tag = "A") + th
st <- fread("06_results/panel/selection_stability_v2.csv")
pan <- rbindlist(lapply(PORD, function(pg) data.table(program = pg, gene_symbol = frz$programs[[pg]]$panel$gene_symbol,
  coefficient = frz$programs[[pg]]$panel$coefficient)))
pan <- merge(pan, st[, .(program, gene_symbol, selection_freq)], by = c("program", "gene_symbol"), all.x = TRUE)
pan[, `:=`(plab = factor(PLAB[program], levels = PLAB), compartment = factor(comp_of[program], levels = names(COMP_COL)),
  glab = paste0(gene_symbol, ifelse(coefficient < 0, " (\u2212)", "")))]
p3b <- ggplot(pan, aes(selection_freq, reorder(glab, selection_freq), fill = compartment)) +
  geom_col(width = 0.68) + geom_text(aes(label = sprintf("%.2f", selection_freq)), hjust = -0.15, size = 2.1, colour = INK2) +
  facet_wrap(~plab, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = COMP_COL, breaks = names(COMP_COL), name = NULL) +
  scale_x_continuous(limits = c(0, 1.32), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  labs(x = "Selection frequency of frozen panel genes across 50 cross-validation fits; (\u2212) = negative coefficient", y = NULL, tag = "B") + th
save_fig(p3a / p3b + plot_layout(heights = c(1.15, 1)), "Figure3_panel_baselines", W, 7.6)
allc <- fread("06_results/cross_cohort_meta/per_cohort_results_long.csv")
pc <- allc[is.finite(MAE_ratio) & MAE_ratio > 0]
mm <- fread("06_results/cross_cohort_meta/meta_analysis_results.csv")[analysis == "main" & is.finite(pooled_MAE_ratio)]
f4 <- rbind(pc[, .(cohort, program, scaling, n = n_donors, est = MAE_ratio, lo = ratio_lo, hi = ratio_hi)],
  mm[, .(cohort = "Pooled", program, scaling, n = total_donors, est = pooled_MAE_ratio, lo = ratio_lo, hi = ratio_hi)])
ne <- allc[!(is.finite(MAE_ratio) & MAE_ratio > 0), .(cohort, program, scaling, lab = "not evaluable")]
ne <- rbind(ne, CJ(program = PORD, scaling = c("frozen", "within"))[!mm, on = .(program, scaling)][, .(cohort = "Pooled", program, scaling, lab = "not pooled (<2 eligible cohorts)")])
SCL <- c(frozen = "Primary: frozen discovery scaling", within = "Secondary: within-cohort scaling")
prep4 <- function(x) x[, `:=`(plab = factor(PLAB[program], levels = PLAB), cohort = factor(cohort, levels = rev(names(COH_COL))),
  scal = factor(SCL[scaling], levels = SCL))][]
f4 <- prep4(f4); ne <- prep4(ne); f4[, small := n < 10 & cohort != "Pooled"]
f4[, lbl := ifelse(cohort == "Pooled", sprintf("%.2f", est), "")]
p4 <- ggplot(f4, aes(est, cohort, colour = cohort)) +
  annotate("rect", xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf, fill = "#eef5ee", alpha = 0.7) +
  geom_vline(xintercept = 1, colour = INK2, linetype = "22", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = pmax(lo, 0.12), xmax = pmin(hi, 8.5), alpha = small), height = 0, linewidth = 0.55) +
  geom_point(aes(shape = cohort == "Pooled", alpha = small, size = cohort == "Pooled")) +
  geom_text(aes(label = lbl), nudge_y = 0.42, size = 2, colour = INK, fontface = "bold") +
  geom_text(data = ne, aes(x = 1, y = cohort, label = lab), inherit.aes = FALSE, colour = "#9a998f", size = 2, fontface = "italic") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.7, `TRUE` = 3), guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.4), guide = "none") +
  scale_colour_manual(values = COH_COL, breaks = setdiff(names(COH_COL), "Pooled"), name = NULL) +
  scale_x_log10(limits = c(0.10, 8.5), breaks = c(0.25, 0.5, 1, 2, 4), labels = c("0.25", "0.5", "1", "2", "4")) +
  facet_grid(plab ~ scal, switch = "y") +
  labs(x = "MAE ratio vs discovery-mean baseline (log scale; 95% bootstrap CI; shaded = panel better)", y = NULL,
    caption = "Pooled (\u25c6) = DerSimonian\u2013Laird random effects across cohorts with \u226510 donors.\nFaded: GSE144735 (<10 donors, descriptive only).") +
  th + theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7), strip.placement = "outside",
    axis.text.y = element_text(size = 6.2), panel.spacing.y = unit(0.25, "lines"))
save_fig(p4, "Figure4_validation_forest", W, 9.2)
sp <- readRDS("04_processed_data/Valdeolivas_Visium/spot_program_scores.rds")
sec <- "SN048_A416371_Rep1"; m <- sp[section == sec & program %in% c("MYE_IFNG", "END_TNFA_NFKB")]
m <- melt(m, id.vars = c("program", "pxl_row", "pxl_col"), measure.vars = c("target", "panel"), variable.name = "score")
m[, value := ave(value, program, score, FUN = function(v) as.numeric(scale(v)))]
m[, `:=`(plab = factor(PLAB[program], levels = PLAB[c("MYE_IFNG", "END_TNFA_NFKB")]),
  score = factor(ifelse(score == "target", "Target (evaluation genes)", "Frozen panel score"), levels = c("Target (evaluation genes)", "Frozen panel score")))]
p5a <- ggplot(m, aes(pxl_col, -pxl_row, colour = pmin(pmax(value, -2.5), 2.5))) + geom_point(size = 0.62, stroke = 0) +
  facet_grid(plab ~ score) + coord_equal() + scale_colour_gradientn(colours = DIV, limits = c(-2.5, 2.5), name = "z-score\n(within section)") +
  labs(x = NULL, y = NULL, tag = "A", title = paste0("Example section ", sec, " (patient A416371)")) +
  th + theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), panel.grid.major.x = element_blank(),
    legend.position = "right", legend.key.height = unit(0.6, "cm"), legend.key.width = unit(0.25, "cm"),
    strip.text.y = element_text(angle = -90, hjust = 0.5))
ap <- fread("06_results/spatial_support_combined/criterionA_patients_combined.csv")
ss <- fread("06_results/spatial_support_combined/spatial_support_summary_main.csv")
COH <- c(Valdeolivas = "Valdeolivas (n=6)", GSE226997 = "GSE226997 (n=4)", GSE280318 = "GSE280318 (n=3)")
ap[, `:=`(plab = factor(PLAB[program], levels = rev(PLAB)), compartment = factor(comp_of[program], levels = names(COMP_COL)))]
ss[, `:=`(plab = factor(PLAB[program], levels = rev(PLAB)), compartment = factor(comp_of[program], levels = names(COMP_COL)))]
p5b <- ggplot() + geom_vline(xintercept = 0, colour = INK2, linetype = "22", linewidth = 0.3) +
  geom_point(data = ap, aes(effect, plab, colour = compartment), alpha = 0.45, size = 1.2, position = position_jitter(height = 0.12, seed = 1)) +
  geom_errorbarh(data = ss, aes(xmin = lo, xmax = hi, y = plab, colour = compartment), height = 0, linewidth = 0.9) +
  geom_point(data = ss, aes(pooled_effect, plab, fill = compartment), shape = 23, size = 2.6, colour = INK, stroke = 0.3) +
  scale_colour_manual(values = COMP_COL, guide = "none") + scale_fill_manual(values = COMP_COL, guide = "none") +
  labs(x = "Composition-adjusted panel coefficient\n(criterion A; \u25c7 pooled, 95% CI; dots = 13 patients)", y = NULL, tag = "B") + th
bp <- fread("06_results/spatial_support_combined/criterionB_patients_combined.csv")[is.finite(diff)]
bp[, `:=`(plab = factor(PLAB[program], levels = rev(PLAB)), compartment = factor(comp_of[program], levels = names(COMP_COL)))]
bp[, cohort := factor(COH[cohort], levels = COH)]

p5c <- ggplot(bp, aes(diff, plab)) + geom_vline(xintercept = 0, colour = INK2, linetype = "22", linewidth = 0.3) +
  geom_point(aes(fill = compartment, shape = cohort), size = 1.8, colour = INK, stroke = 0.25, alpha = 0.9, position = position_jitter(height = 0.12, seed = 2)) +
  scale_fill_manual(values = COMP_COL, guide = "none") + scale_shape_manual(values = c(21, 24, 22), name = NULL, drop = FALSE) +
  guides(shape = guide_legend(nrow = 1, override.aes = list(fill = "#9a998f"))) +
  labs(x = "Target difference, expected-high minus comparison\n(criterion B; dots = patients)", y = NULL, tag = "C") +
  th + theme(axis.text.y = element_blank(), legend.position = "bottom")
save_fig(p5a / (p5b | p5c) + plot_layout(heights = c(1.15, 1)), "Figure5_spatial_support", W, 7.8)
ev <- fread("08_tables/Table_program_evidence.csv")
ordp <- PORD[order(match(ev$evidence_category[match(PORD, ev$program)], CAT))]
ev[, `:=`(plab = factor(PLAB[program], levels = rev(PLAB[ordp])), compartment = factor(comp_of[program], levels = names(COMP_COL)))]
ev[, catcol := CATCOL[evidence_category]]
long <- rbind(ev[, .(plab, col = "Panel genes", txt = gsub(", ", " \u00b7 ", panel))],
  ev[, .(plab, col = "Single-cell (frozen)", txt = frozen_MAE_ratio_CI)],
  ev[, .(plab, col = "Within-cohort", txt = ifelse(within_scaling_supported, "Supported", "\u2014"))],
  ev[, .(plab, col = "Spatial", txt = spatial_label)])
XPOS <- c(`Panel genes` = 1.0, `Single-cell (frozen)` = 3.2, `Within-cohort` = 4.4, Spatial = 5.2)
long[, x := XPOS[col]]
p6 <- ggplot() +
  geom_tile(data = ev, aes(x = -0.12, y = plab, fill = compartment), width = 0.12, height = 0.8) +
  scale_fill_manual(values = COMP_COL, breaks = names(COMP_COL), name = NULL) +
  geom_tile(data = ev, aes(x = 0.45, y = plab), fill = ev$catcol, width = 0.85, height = 0.8) +
  geom_text(data = long, aes(x = x, y = plab, label = txt), hjust = 0, size = 2.15, colour = INK) +
  annotate("text", x = c(0.45, XPOS), y = length(PLAB) + 0.75, hjust = c(0.5, rep(0, 4)), size = 2.4, fontface = "bold", colour = INK,
    label = c("Evidence", names(XPOS))) +
  scale_x_continuous(limits = c(-0.2, 6.45), expand = c(0, 0)) + coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL, caption = paste0("Evidence colour: dark green = ", CAT[1], ";\nlight green = ", CAT[2], "; amber = ", CAT[3], "; grey = ", CAT[4],
    ".\nSingle-cell: frozen-scaling MAE ratio vs discovery mean (random-effects pooled; endothelial programs: single eligible cohort GSE132465).\nLeft stripe: compartment.")) +
  th + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), panel.grid = element_blank(), panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 7.5, colour = INK, face = "bold"), plot.margin = margin(14, 6, 4, 4))
save_fig(p6, "Figure6_evidence_summary", W, 3.9)
cat("\nAll figures written to 07_figures/ (PDF vector, TIFF 600 dpi, PNG 300 dpi)\n")
 
