# 17b_check_reference_donor_balance.R  (diagnostic, read-only; added 19 Oct 2026)
# Does the RCTD reference built by script 17 over-represent donors stored early in the
# GSE178341 h5 file? Script 17 samples 1,200 cells per type at random, re-sorts them by
# h5 position, then keeps the first 800 that pass QC. If h5 order follows patient or
# channel, later donors are under-sampled.
#
# For every reference type this compares each donor's share of the SELECTED cells with
# its share of the ELIGIBLE pool (all tumour cells of that type), and correlates a
# donor's selection rate with the donor's median position in the h5 file.
# Under random selection the rate is roughly constant across donors and the correlation
# is near 0. Reads metadata, h5 barcodes and the saved reference only (a few minutes).
suppressPackageStartupMessages({ library(data.table); library(hdf5r) })
in_dir <- "03_original_data/GSE178341"
ref <- readRDS("04_processed_data/spatial_reference/GSE178341_reference_counts.rds")
meta <- fread(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"))
clus <- fread(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"))
d <- merge(meta, clus, by.x = "cellID", by.y = "sampleID")
Tm <- d$SPECIMEN_TYPE == "T"
d[, ref_type := fcase(                                  # verbatim from script 17
  clMidwayPr == "Epi" & Tm & cl295v11SubShort %in% c("cE01", "cE02", "cE03"), "Malignant_epithelial",
  clMidwayPr == "Epi" & SPECIMEN_TYPE == "N", "Normal_epithelial",
  clMidwayPr == "Fibro" & Tm, "Fibroblast",
  clMidwayPr %in% c("Peri", "SmoothMuscle") & Tm, "Pericyte_SMC",
  clMidwayPr == "Endo" & Tm, "Endothelial",
  clMidwayPr %in% c("Mono", "Macro", "DC", "Granulo") & Tm, "Myeloid",
  grepl("Mast", clMidwayPr) & Tm, "Mast",
  clMidwayPr %in% c("B", "Plasma") & Tm, "B_Plasma",
  grepl("^T|NK|ILC", clMidwayPr) & Tm, "T_NK", default = NA_character_)]
h5 <- H5File$new(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit.h5"), mode = "r")
bc <- h5[["matrix/barcodes"]]$read(); h5$close_all()
pool <- d[!is.na(ref_type), .(cellID, PID, ref_type)][, h5_pos := match(cellID, bc)]
sel  <- as.data.table(ref$meta)[, .(cellID, PID, ref_type)]

out <- rbindlist(lapply(sort(unique(pool$ref_type)), function(tp) {
  p <- pool[ref_type == tp]; s <- sel[ref_type == tp]
  dn <- p[, .(n_pool = .N, med_pos = median(h5_pos)), by = PID]
  dn <- merge(dn, s[, .(n_sel = .N), by = PID], by = "PID", all.x = TRUE); dn[is.na(n_sel), n_sel := 0L]
  dn[, rate := n_sel / n_pool]
  tvd <- 0.5 * sum(abs(dn$n_sel / sum(dn$n_sel) - dn$n_pool / sum(dn$n_pool)))
  data.table(ref_type = tp, pool_cells = nrow(p), selected = nrow(s),
             donors_pool = nrow(dn), donors_selected = sum(dn$n_sel > 0),
             tvd_vs_pool = round(tvd, 3),
             rho_rate_vs_h5_position = round(suppressWarnings(cor(dn$rate, dn$med_pos, method = "spearman")), 2),
             affected = nrow(p) > 800)                # truncation to 800 passing cells only bites above this
}))
dir.create("06_results/spatial_reference_check", recursive = TRUE, showWarnings = FALSE)
fwrite(out, "06_results/spatial_reference_check/reference_donor_balance.csv")
cat("\n=== RCTD reference: donor balance vs eligible pool ===\n"); print(out)
cat("\nInterpretation: for 'affected' types a clearly negative rho (e.g. < -0.3) together with\n",
    "fewer donors_selected than donors_pool indicates the h5-order bias. If present, rebuild the\n",
    "reference with FIX_SELECTION_ORDER <- TRUE in script 17 and repeat 18/18b/19/22/23 as a\n",
    "sensitivity analysis.\n")

dir.create("06_results/sessionInfo", recursive = TRUE, showWarnings = FALSE); writeLines(capture.output(sessionInfo()), "06_results/sessionInfo/sessionInfo_17b.txt")
