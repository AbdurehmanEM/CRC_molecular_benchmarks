# 07_GSE178341_compartment_counts.R
# Metadata only: map clusters to study compartments; count cells per donor x
# compartment x processing type, to fix specimen and sorting rules first.

suppressPackageStartupMessages(library(data.table))
in_dir  <- "03_original_data/GSE178341"
res_dir <- "06_results/GSE178341_audit"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(res_dir, "compartment_counts_log.txt"), split = TRUE)

meta <- fread(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"))
clus <- fread(file.path(in_dir, "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"))
d <- merge(meta, clus, by.x = "cellID", by.y = "sampleID")
stopifnot(nrow(d) == nrow(meta))

d[, compartment := fcase(
  clMidwayPr == "Endo" & grepl("lymphatic", cl295v11SubFull, ignore.case = TRUE), "Lymphatic_EC",
  clMidwayPr == "Endo",                              "Endothelial",
  clMidwayPr == "Fibro",                             "Stromal",
  clMidwayPr %in% c("Peri", "SmoothMuscle"),         "Mural_other",
  clMidwayPr == "Schwann",                           "Glia",
  clMidwayPr %in% c("Mono", "Macro"),                "Myeloid",
  clMidwayPr %in% c("DC", "Granulo"),                "Myeloid_other",
  clMidwayPr == "Epi",                               "Epithelial",
  default = "Other_context")]
fwrite(d[, .N, by = .(clMidwayPr, cl295v11SubFull, compartment)][order(compartment, -N)],
       file.path(res_dir, "cluster_to_compartment_map.csv"))

tum  <- d[SPECIMEN_TYPE == "T"]
core <- c("Epithelial", "Stromal", "Myeloid", "Endothelial")

cat("\n=== Tumour specimens by processing types present ===\n")
print(tum[, .(types = paste(sort(unique(PROCESSING_TYPE)), collapse = " + ")), by = PatientTypeID][, .N, by = types])

cat("\n=== Core-compartment tumour cells by processing type ===\n")
print(dcast(tum[compartment %in% core, .N, by = .(compartment, PROCESSING_TYPE)],
            compartment ~ PROCESSING_TYPE, value.var = "N", fill = 0))

count_tab <- function(x) dcast(x[compartment %in% core, .N, by = .(PID, compartment)],
                               PID ~ compartment, value.var = "N", fill = 0)
ruleA <- count_tab(tum)
ruleB <- count_tab(tum[PROCESSING_TYPE %in% c("unsorted", "LiveMACS")])
fwrite(ruleA, file.path(res_dir, "donor_compartment_counts_ruleA_all_processing.csv"))
fwrite(ruleB, file.path(res_dir, "donor_compartment_counts_ruleB_unsorted_only.csv"))

elig <- function(tab) rbindlist(lapply(c(10, 20, 30, 50), function(k)
  data.table(min_cells = k, t(colSums(as.matrix(tab[, ..core]) >= k)),
             all_four = sum(rowSums(as.matrix(tab[, ..core]) >= k) == 4))))
cat("\n=== Eligible donors (of", uniqueN(tum$PID), ") - Rule A: all processing pooled ===\n"); print(elig(ruleA))
cat("\n=== Eligible donors - Rule B: unsorted + LiveMACS only ===\n");                     print(elig(ruleB))

batch <- unique(tum[, .(PID, SOURCE_HOSPITAL, SINGLECELL_TYPE, MMRStatus, TissueSiteSimple, MetastasisStatus)])
cat("\n=== Donors with >1 chemistry across tumour channels ===\n"); print(batch[, .N, by = PID][N > 1])
fwrite(batch, file.path(res_dir, "donor_batch_clinical.csv"))

epi <- dcast(d[compartment == "Epithelial", .N, by = .(cl295v11SubFull, SPECIMEN_TYPE)],
             cl295v11SubFull ~ SPECIMEN_TYPE, value.var = "N", fill = 0)
epi[, pct_in_tumour := round(100 * T / (T + N), 1)]
cat("\n=== Epithelial subclusters: tumour (T) vs normal (N) specimens ===\n")
print(epi[order(-pct_in_tumour)])
fwrite(epi, file.path(res_dir, "epithelial_subclusters_T_vs_N.csv"))
sink()