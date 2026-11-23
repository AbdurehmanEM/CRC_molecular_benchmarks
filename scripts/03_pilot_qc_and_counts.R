# 03_pilot_qc_and_counts.R
# Pilot QC plots, sample-aware outlier flags, compartment mapping and the
# donor x compartment cell-count table that decides go/no-go.
# Run twice: the output hash must be identical (reproducibility check).

suppressPackageStartupMessages({ library(Matrix); library(data.table); library(ggplot2); library(digest) })
set.seed(20261008)

tag     <- "full"            # "full" after the pilot passes
in_dir  <- "04_processed_data/GSE132465"
res_dir <- if (tag == "pilot") "06_results/pilot" else "06_results/GSE132465"
fig_dir <- if (tag == "pilot") "07_figures/pilot" else "07_figures/GSE132465"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(in_dir, paste0("GSE132465_", tag, "_counts_sparse.rds")))
meta   <- as.data.table(readRDS(file.path(in_dir, paste0("GSE132465_", tag, "_cell_meta.rds"))))
stopifnot(identical(meta$Index, colnames(counts)))

# ---------------------------------------------------------------- QC metrics
mt <- grepl("^MT-", rownames(counts))
meta[, nCount   := Matrix::colSums(counts)]
meta[, nFeature := Matrix::colSums(counts > 0)]
meta[, pct_mt   := 100 * Matrix::colSums(counts[mt, , drop = FALSE]) / nCount]

# Sample-aware flags: 3 MADs from each SAMPLE's median (flags for inspection)
mad_flag <- function(x, n = 3, side = c("both", "high")) {
  side <- match.arg(side); m <- median(x); d <- mad(x)
  if (side == "high") x > m + n * d else x < m - n * d | x > m + n * d
}
meta[, `:=`(
  flag_lowcount = mad_flag(log10(nCount + 1)),
  flag_features = mad_flag(log10(nFeature + 1)),
  flag_mt       = mad_flag(pct_mt, side = "high")
), by = Sample]
meta[, qc_flag := flag_lowcount | flag_features | flag_mt]

qc_summary <- meta[, .(cells = .N, median_counts = as.numeric(median(nCount)),
                       median_features = as.numeric(median(nFeature)),
                       median_pct_mt = round(median(pct_mt), 2),
                       max_pct_mt = round(max(pct_mt), 2),
                       flagged = sum(qc_flag)), by = .(Patient, Sample, Class)]
fwrite(qc_summary, file.path(res_dir, "qc_summary_per_sample.csv"))
print(qc_summary)

# ---------------------------------------------------------------- QC plots
p1 <- ggplot(meta, aes(Sample, nCount, fill = Class)) + geom_violin(scale = "width") +
  scale_y_log10() + theme_bw(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "UMIs per cell (log10)", x = NULL)
p2 <- ggplot(meta, aes(Sample, nFeature, fill = Class)) + geom_violin(scale = "width") +
  scale_y_log10() + theme_bw(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Genes per cell (log10)", x = NULL)
p3 <- ggplot(meta, aes(Sample, pct_mt, fill = Class)) + geom_violin(scale = "width") +
  theme_bw(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "% mitochondrial UMIs", x = NULL)
p4 <- ggplot(meta, aes(nCount, nFeature, colour = pct_mt)) + geom_point(size = 0.2) +
  scale_x_log10() + scale_y_log10() + facet_wrap(~Sample) + theme_bw(9) +
  scale_colour_viridis_c() + labs(x = "UMIs", y = "Genes")
ggsave(file.path(fig_dir, "qc_umis_per_sample.pdf"), p1, width = 8, height = 4)
ggsave(file.path(fig_dir, "qc_genes_per_sample.pdf"), p2, width = 8, height = 4)
ggsave(file.path(fig_dir, "qc_mito_per_sample.pdf"), p3, width = 8, height = 4)
ggsave(file.path(fig_dir, "qc_umis_vs_genes.png"),   p4, width = 10, height = 7, dpi = 200)

# ---------------------------------------------------------------- compartments
type_tab <- meta[, .N, by = .(Cell_type, Cell_subtype)][order(Cell_type, -N)]
fwrite(type_tab, file.path(res_dir, "published_celltype_subtype_table.csv"))
print(type_tab, nrows = 200)

# PROVISIONAL mapping - verify against type_tab and the marker plot
meta[, compartment := fcase(
  Cell_subtype %in% c("Tip-like ECs", "Stalk-like ECs", "Proliferative ECs"), "Endothelial",
  Cell_subtype == "Lymphatic ECs",                                       "Lymphatic_EC",
  Cell_subtype %in% c("Pericytes", "Smooth muscle cells"),              "Mural_other",
  Cell_subtype == "Enteric glial cells",                                "Glia",
  Cell_subtype %in% c("Stromal 1", "Stromal 2", "Stromal 3", "Myofibroblasts"), "Stromal",
  Cell_type == "Myeloids",                                              "Myeloid",
  Cell_type == "Epithelial cells",                                      "Epithelial",
  default = "Other_context"
)]
meta[compartment == "Epithelial", epi_status := fcase(
  Class != "Tumor", "adjacent_mucosa_epithelium",
  grepl("CMS", Cell_subtype), "putative_malignant",
  default = "uncertain_in_tumor")]

# Marker check (dot plot of canonical markers by compartment)
markers <- c("EPCAM","KRT8","KRT20","CEACAM5","COL1A1","DCN","PDGFRA","ACTA2","RGS5",
             "PECAM1","VWF","CDH5","PLVAP","LYZ","CD68","CD14","C1QA","PTPRC","CD3E","MS4A1")
markers <- intersect(markers, rownames(counts))
lib <- Matrix::colSums(counts)
lognorm <- log1p(t(t(counts[markers, , drop = FALSE]) / lib * 1e4))
dot <- rbindlist(lapply(split(seq_len(ncol(counts)), meta$compartment), function(ix) {
  data.table(gene = markers,
             mean_expr = Matrix::rowMeans(lognorm[, ix, drop = FALSE]),
             pct_expr  = 100 * Matrix::rowMeans(lognorm[, ix, drop = FALSE] > 0))
}), idcol = "compartment")
p5 <- ggplot(dot, aes(factor(gene, levels = markers), compartment, size = pct_expr, colour = mean_expr)) +
  geom_point() + scale_colour_viridis_c() + theme_bw(10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = NULL, y = NULL)
ggsave(file.path(fig_dir, "marker_dotplot_by_compartment.pdf"), p5, width = 9, height = 3.5)

# ---------------------------------------------------------------- donor x compartment counts
retained <- meta[qc_flag == FALSE & Class == "Tumor"]
cnt <- dcast(retained[, .N, by = .(Patient, Sample, compartment)],
             Patient + Sample ~ compartment, value.var = "N", fill = 0)
fwrite(cnt, file.path(res_dir, "cell_counts_donor_by_compartment_tumor.csv"))
print(cnt)
epi_cnt <- dcast(retained[compartment == "Epithelial", .N, by = .(Patient, epi_status)],
                 Patient ~ epi_status, value.var = "N", fill = 0)
fwrite(epi_cnt, file.path(res_dir, "epithelial_status_by_donor_tumor.csv"))
print(epi_cnt)

MIN_CELLS <- 30
core <- intersect(c("Epithelial", "Stromal", "Myeloid", "Endothelial"), names(cnt))
elig <- cnt[, lapply(.SD, function(x) x >= MIN_CELLS), .SDcols = core]
cat("\nDonors with >=", MIN_CELLS, "cells per compartment (tumour, post-QC flags):\n")
print(colSums(elig))
fwrite(cbind(cnt[, .(Patient, Sample)], elig), file.path(res_dir, "pseudobulk_eligibility_min30.csv"))

fwrite(meta, file.path(res_dir, "cell_meta_with_qc_and_compartment.csv.gz"))

# ---------------------------------------------------------------- reproducibility
h <- digest(list(qc_summary, cnt, type_tab))
cat("\nOutput hash:", h, "\n")
cat(format(Sys.time()), h, "\n", file = file.path(res_dir, "reproducibility_hashes.txt"), append = TRUE)
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_03.txt"))
