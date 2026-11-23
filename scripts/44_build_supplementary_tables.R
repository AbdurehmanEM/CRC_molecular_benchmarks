# 44_build_supplementary_tables.R - submission-ready supplementary workbook
#
# Replaces section 4 of 41_build_submission_with_lymphoid.R. Reads saved results only (no analysis
# is re-run) and writes 10_submission/Supplementary_Tables_aem.xlsx with:
#   * a Contents sheet (title, table list, legends, abbreviations)
#   * Supplementary Tables S1-S13, each with a title row and legend above the data
#   * a References sheet at the end
# Changes against the 19-sheet workbook from script 41:
#   - S1 rebuilt as the final eligibility table (the working audit sheets still carried VERIFY/Pending
#     placeholders and wrong roles); added per-compartment coverage (S1b) and the validation patient
#     roster with distinct-patient counts (S1c)
#   - S2 now gives recruitment site, submitter and the discovery-donor clinical characteristics
#   - S3 built from the per-stage sessionInfo files, so spacexr, SpatialExperiment etc. are included
#   - S10d now covers all spatial sections (the old sheet held the added cohorts only)
#   - S13 added: GSE200997 cluster-to-compartment assignment (cited for the marker-rule annotation)
#   - readable column names, program labels, 3-decimal rounding
# Run from the project root:  source("05_analysis_scripts/44_build_supplementary_tables.R")

suppressPackageStartupMessages({ library(data.table); library(openxlsx) })
OUT   <- "10_submission/Supplementary_Tables_aem.xlsx"
TITLE <- "Evidence-Graded Molecular Benchmarks for Colorectal Cancer-on-a-Chip Models Derived from Single-Cell and Spatial Transcriptomics"
AUTH  <- "Abdurehman Eshete Mohammed"
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

rd <- function(p) { if (!file.exists(p)) stop("missing input: ", p); fread(p, na.strings = c("", "NA")) }
rd2 <- function(a, b) rbind(rd(a), rd(b), fill = TRUE)

PL <- c(EPI_HYPOXIA = "Epithelial hypoxia", EPI_APOPTOSIS = "Epithelial apoptosis",
        STR_TGFB = "Stromal TGF-β response", MYE_IFNG = "Myeloid IFN-γ response",
        MYE_INFLAM = "Myeloid inflammatory response", END_ANGIO = "Endothelial angiogenesis",
        END_TNFA_NFKB = "Endothelial TNF/NF-κB activation", LYM_IL2_STAT5 = "CD8 T-cell IL2/STAT5", STR_ECM = "Stromal ECM organisation (not retained)", LYM_CYTOTOX = "CD8 T-cell cytotoxicity (not retained)")
ORD <- names(PL)
add_label <- function(d) {                       # program label next to the program code, rows in fixed order
  if (!"program" %in% names(d)) return(d)
  d <- copy(d)[, Program := ifelse(program %in% names(PL), PL[program], program)]
  setcolorder(d, c("program", "Program"))
  d[, .o := match(program, ORD)]; setorder(d, .o, na.last = TRUE); d[, .o := NULL]
  setnames(d, "program", "Program code"); d
}
tidy <- function(d, digits = 3) {                # round non-integer numerics; logicals to Yes/No
  d <- copy(d)
  for (j in names(d)) {
    v <- d[[j]]
    if (is.logical(v)) set(d, j = j, value = ifelse(v, "Yes", "No"))
    else if (is.numeric(v) && any(abs(v - round(v)) > 1e-9, na.rm = TRUE)) set(d, j = j, value = round(v, digits))
  }
  d
}
ren <- function(d, map) { m <- map[names(map) %in% names(d)]; setnames(d, names(m), unname(m)); d }

## ---- S1a: cohort eligibility (final) ------------------------------------------------------------
pred <- list(
  GSE132465 = rbind(rd("06_results/validation_GSE132465/validation_predictions.csv")[scaling == "frozen", .(program, donor = Patient)],
                    rd("06_results/validation_GSE132465/validation_predictions_lymphoid.csv")[scaling == "frozen", .(program, donor = Patient)]),
  GSE166555 = rbind(rd("06_results/validation_GSE166555/validation_predictions.csv")[, .(program, donor)],
                    rd("06_results/validation_lymphoid/validation_predictions_extra_cohorts.csv")[cohort == "GSE166555", .(program = "LYM_IL2_STAT5", donor)]),
  GSE200997 = rd("06_results/validation_GSE200997/validation_predictions.csv")[, .(program, donor)])
roster <- rbindlist(lapply(names(pred), function(k) unique(pred[[k]])[, cohort := k]))
n_dist <- roster[, .(n = uniqueN(donor)), by = cohort]; nd <- setNames(n_dist$n, n_dist$cohort)
spA <- rbind(rd("06_results/spatial_support_combined/criterionA_patients_combined.csv"),
             rd("06_results/spatial_support_lymphoid_combined/criterionA_patients_combined.csv"))
sp_n <- spA[, .(n = uniqueN(patient)), by = cohort]; spn <- setNames(sp_n$n, sp_n$cohort)
kul <- rd("06_results/validation_GSE144735/KUL_gene_coverage.csv")
dm  <- rd("04_processed_data/GSE178341/donor_metadata_tumour.csv")

S1a <- data.table(
  Dataset = c("GSE178341", "GSE132465", "GSE166555", "GSE200997", "GSE144735",
              "Valdeolivas et al. (Zenodo 7739700)", "GSE280318", "GSE226997"),
  `Source publication` = c("Pelka et al. (2021)", "Lee et al. (2020)", "Uhlitz et al. (2021)", "Khaliq et al. (2022)",
                           "Lee et al. (2020)", "Valdeolivas et al. (2024)", "Oliveira et al. (2025)", "Park et al. (2024)"),
  Country = c("USA", "South Korea", "Germany", "USA", "Belgium", "Germany (tissue source per publication)", "USA", "South Korea"),
  Technology = c("scRNA-seq, 10x Chromium 3′ v2 and v3", "scRNA-seq, 10x Chromium 3′ v2",
                 "scRNA-seq, 10x Chromium 3′ v3", "scRNA-seq, 10x Chromium 5′",
                 "scRNA-seq, 10x Chromium 3′ v2", "Visium (fresh-frozen)", "Visium HD (FFPE), 8-µm bins aggregated to 56-µm pseudo-spots",
                 "Visium for FFPE"),
  `Cell selection` = c("Unsorted, CD45⁺ MACS and live-cell MACS fractions (sorting fraction recorded per channel)",
                       "Unsorted", "Unsorted", "Viable singlets (FACS, DAPI exclusion); no lineage enrichment",
                       "Unsorted", "Not applicable (intact sections)", "Not applicable (intact sections)", "Not applicable (intact sections)"),
  `Specimens in dataset` = c(sprintf("%d patients; %d tumour and 36 adjacent-normal specimens", nrow(dm), sum(dm$n_tumour_specimens)),
                             "23 patients; 23 tumour and 10 adjacent-mucosa specimens",
                             "12 patients; tumour and adjacent normal", "16 patients; 16 tumour and 8 adjacent normal",
                             "6 patients; tumour core, border and adjacent normal", "7 patients; 14 sections",
                             "5 patients; Visium HD sections for 3 tumours (P1, P2, P5)", "4 patients; 4 sections"),
  `Specimens analysed` = c("Tumour only; unsorted and live-cell MACS fractions (CD45+ and mixed fractions excluded); two patients with two tumour specimens pooled", "Tumour only", "Tumour only", "Tumour only",
                           "Tumour only", "Tumour sections; patient A798015 (non-neoplastic tissue) excluded",
                           "Tumour sections", "Tumour sections"),
  `Treatment status` = c(rep("Untreated at resection (source publication)", 7), "Not reported in GEO; retained as pre-specified (results shown with and without this cohort)"),
  `Count data` = c("Raw UMI counts (10x HDF5)", "Raw UMI matrix (normalised TPM file not used)", "Raw UMI counts", "Raw UMI counts",
                   "Raw UMI matrix", "Space Ranger filtered counts", "Space Ranger filtered counts (binned outputs)", "Space Ranger filtered counts"),
  `Cell annotation` = c("Author labels", "Author labels", "Author labels (no macrophage label; myeloid = monocytes)",
                        "None; pre-specified marker rule on Seurat clusters (Supplementary Table S13)",
                        "Author labels", "Pathologist spot annotations", "None used; compartment-proportion quartiles", "None used; compartment-proportion quartiles"),
  `Study role` = c("Discovery", "Independent validation", "Independent validation", "Independent validation",
                   "Descriptive validation (shares a publication with GSE132465)", "Spatial assessment (primary)",
                   "Spatial assessment (added)", "Spatial assessment (added)"),
  `Patients analysed` = c(nrow(dm), nd["GSE132465"], nd["GSE166555"], nd["GSE200997"], max(kul$kul_donors),
                          spn["Valdeolivas"], spn["GSE280318"], spn["GSE226997"]),
  `Lymphoid (CD8) analysis` = c("Included", "Included", "Included", "Excluded: T and NK cells annotated jointly",
                                "Attempted; 2 patients (not evaluable)", "Included", "Included", "Included"))
S1a[, `Patients analysed` := as.integer(`Patients analysed`)]
S1a_note <- paste("Patients analysed: single-cell cohorts, distinct patients contributing an eligible profile (≥20 cells in a tumour specimen) to at least one program;",
                  "spatial cohorts, patients with at least one analysed section. Considered but not audited: Joanito et al. (2022; controlled access, likely overlap with validation cohorts)",
                  "and GSE146771 (CD45-sorted, mixed platforms). No shared donors were identified between cohorts.")

## ---- S1b: eligible donors per cohort x compartment -------------------------------------------------
cmp_of <- c(EPI_HYPOXIA = "Epithelial", EPI_APOPTOSIS = "Epithelial", STR_TGFB = "Stromal", MYE_IFNG = "Myeloid",
            MYE_INFLAM = "Myeloid", END_ANGIO = "Endothelial", END_TNFA_NFKB = "Endothelial", LYM_IL2_STAT5 = "CD8 T cells")
el <- rd("06_results/GSE178341/pseudobulk_eligibility.csv")[min_cells == 20]
ell <- rd("06_results/GSE178341/pseudobulk_eligibility_lymphoid.csv")[min_cells == 20]
cov <- rbind(
  data.table(cohort = "GSE178341", compartment = c("Epithelial", "Stromal", "Myeloid", "Endothelial", "CD8 T cells"),
             n = c(el$Epithelial, el$Stromal, el$Myeloid, el$Endothelial, ell$Lymphoid_CD8)),
  rd("06_results/validation_GSE132465/validation_gene_coverage.csv")[, .(cohort = "GSE132465", compartment = cmp_of[program], n = n_donors)],
  rd("06_results/validation_GSE132465/validation_gene_coverage_lymphoid.csv")[, .(cohort = "GSE132465", compartment = "CD8 T cells", n = n_donors)],
  rd("06_results/validation_GSE166555/validation_gene_coverage.csv")[, .(cohort = "GSE166555", compartment = cmp_of[program], n = n_donors)],
  rd("06_results/validation_GSE200997/validation_gene_coverage.csv")[, .(cohort = "GSE200997", compartment = cmp_of[program], n = n_donors)],
  kul[, .(cohort = "GSE144735", compartment = cmp_of[program], n = kul_donors)],
  rd("06_results/validation_lymphoid/validation_coverage_extra_cohorts.csv")[, .(cohort = sub("_KUL", "", cohort), compartment = "CD8 T cells", n = n_donors)])
cov <- unique(cov)
S1b <- dcast(cov, cohort ~ compartment, value.var = "n", fun.aggregate = max, fill = NA)
S1b[cohort == "GSE200997", `CD8 T cells` := NA]
setcolorder(S1b, c("cohort", intersect(c("Epithelial", "Stromal", "Myeloid", "Endothelial", "CD8 T cells"), names(S1b))))
S1b <- S1b[match(c("GSE178341", "GSE132465", "GSE166555", "GSE200997", "GSE144735"), cohort)]
S1b[, Role := c("Discovery", "Validation", "Validation", "Validation", "Descriptive")]
setnames(S1b, "cohort", "Cohort")
S1b_note <- "Eligible donors: at least 20 quality-controlled cells of the compartment in tumour specimens. Profiles from cohorts with fewer than 10 eligible donors were not evaluable for validation. Blank: not definable (GSE200997 annotates T and NK cells jointly)."

## ---- S1c: validation patient roster -----------------------------------------------------------------
S1c <- roster[, .(`Patients (n)` = uniqueN(donor), `Patient identifiers` = paste(sort(unique(donor)), collapse = ", ")),
              by = .(Cohort = cohort, program)]
S1c <- add_label(S1c)
tot <- roster[!program %in% "LYM_IL2_STAT5", .(n = uniqueN(donor)), by = .(cohort, program)]
per_prog_max <- tot[, .(n = sum(n)), by = program][, max(n)]
S1c_note <- sprintf(paste("Frozen-scaling analysis; GSE144735 (descriptive) not included. Distinct independent validation patients: %d",
                          "(GSE132465 %d, GSE166555 %d, GSE200997 %d). The largest total for a single program is %d",
                          "(epithelial programs: 21 + 12 + 12); myeloid 44, lymphoid 32, stromal 23 and endothelial 11 (GSE132465 only)."),
                    sum(nd), nd["GSE132465"], nd["GSE166555"], nd["GSE200997"], per_prog_max)

## ---- S2: provenance and discovery-donor characteristics ---------------------------------------------
pv <- rd("02_metadata/cohort_provenance.csv")
site <- c(GSE178341 = "Massachusetts General Hospital and Dana-Farber Cancer Institute, Boston", GSE132465 = "Samsung Medical Center, Seoul",
          GSE144735 = "KU Leuven, Leuven", GSE166555 = "Charité – Universitätsmedizin Berlin", GSE200997 = "Indiana University, Indianapolis",
          GSE280318 = "10x Genomics (commercially sourced tissue)", GSE226997 = "Ajou University School of Medicine, Suwon")
S2a <- pv[, .(Dataset = gse, `Recruitment site (source publication)` = site[gse],
              `GEO submitting institution` = paste0(institute, ", ", tools::toTitleCase(tolower(city))), Country = country, PMID = pubmed)]
ctry <- c(GSE178341 = "USA", GSE132465 = "South Korea", GSE144735 = "Belgium", GSE166555 = "Germany",
          GSE200997 = "USA", GSE280318 = "USA", GSE226997 = "South Korea")
S2a[, `GEO submitter country` := Country][, Country := ctry[Dataset]]
setcolorder(S2a, c("Dataset", "Recruitment site (source publication)", "Country", "GEO submitting institution", "GEO submitter country", "PMID"))
S2a <- rbind(S2a, data.table(Dataset = "Valdeolivas et al. (Zenodo 7739700)",
                             `Recruitment site (source publication)` = "Human Tissue and Cell Research (HTCR) Foundation, Munich",
                             `GEO submitting institution` = "Not applicable (Zenodo)", `GEO submitter country` = NA_character_,
                             Country = "Germany", PMID = NA_integer_), fill = TRUE)
S2a_note <- "Cohort independence was judged from the recruitment site in each source publication, not the depositing institution (GSE144735 was deposited by Samsung Medical Center but recruited at KU Leuven)."
S2b <- dm[, .(Patient = PID, Hospital = hospital, `MMR status` = MMRStatus, `Tumour site` = site, `Metastasis at surgery` = metastasis,
              Sex = sex, Age = age, Chemistry = chemistry, `Tumour specimens` = n_tumour_specimens,
              `Epithelial cells` = n_Epithelial, `Stromal cells` = n_Stromal, `Myeloid cells` = n_Myeloid, `Endothelial cells` = n_Endothelial)]
cd8 <- rd("04_processed_data/GSE178341/donor_CD8_cellcounts.csv")
cd8_col <- intersect(c("n_CD8", "n_cells", "n_primary", "CD8"), names(cd8))[1]
if (!is.na(cd8_col)) S2b <- merge(S2b, cd8[, .(Patient = get(names(cd8)[1]), `CD8 T cells` = get(cd8_col))], by = "Patient", all.x = TRUE)
S2b_note <- "Discovery cohort (GSE178341), quality-controlled cells in tumour specimens per compartment (primary sorting fractions). Clinical fields as deposited."

## ---- S3: software ------------------------------------------------------------------------------------
si <- list.files("06_results", pattern = "^sessionInfo_.*\\.txt$", recursive = TRUE, full.names = TRUE)
KEY <- c("Matrix", "data.table", "Seurat", "SeuratObject", "edgeR", "limma", "glmnet", "msigdbr", "spacexr",
         "SpatialExperiment", "SummarizedExperiment", "SingleCellExperiment", "hdf5r", "rhdf5", "arrow", "GEOquery",
         "digest", "ggplot2", "patchwork", "openxlsx")
sw <- rbindlist(lapply(si, function(f) { x <- readLines(f, warn = FALSE)
  m <- regmatches(x, gregexpr("[A-Za-z][A-Za-z0-9.]*_[0-9][0-9.-]*", x))
  v <- unlist(m); if (!length(v)) return(NULL)
  data.table(package = sub("_.*$", "", v), version = sub("^[^_]*_", "", v), stage = sub("^sessionInfo_(.*)\\.txt$", "\\1", basename(f))) }))
sw <- unique(sw[package %in% KEY])[, .(Version = paste(sort(unique(version)), collapse = "; "),
                                       `Pipeline stages` = paste(sort(unique(stage)), collapse = ", ")), by = .(Package = package)]
rv <- unique(unlist(lapply(si, function(f) sub("^R version ([0-9.]+).*$", "\\1", grep("^R version", readLines(f, 1), value = TRUE)))))
gv <- rd("04_processed_data/programs/geneset_version.csv")
S3 <- rbind(data.table(Package = "R", Version = paste(rv, collapse = "; "), `Pipeline stages` = "all"),
            sw[order(match(Package, KEY))],
            data.table(Package = "MSigDB (via msigdbr)", Version = gv$db_version[1], `Pipeline stages` = paste("retrieved", gv$retrieved[1])))
S3_note <- paste("Versions parsed from the sessionInfo() record written beside each stage's output; stage numbers refer to scripts in 05_analysis_scripts.","Seurat and edgeR were updated during the project (Seurat 5.5.1 to 5.6.0; edgeR 4.10.5 to 4.10.6). All primary analyses, including the GSE200997 clustering (stage 14), used Seurat 5.5.1 and edgeR 4.10.5;","later versions appear only for re-runs of stages that use edgeR::cpm() (unchanged between these releases), the random-reference sensitivity analysis and figure rebuilding.","Versions for the download stages (01, 05) were recorded after the downloads.")

## ---- S4-S12: analysis results -------------------------------------------------------------------------
S4 <- rd2("04_processed_data/programs/program_definitions.csv", "04_processed_data/programs/program_definitions_lymphoid.csv"); lc <- rd("06_results/programs/lymphoid_coverage.csv")[retained == FALSE]; if (nrow(lc)) S4 <- rbind(S4, lc[, .(program, compartment = "Lymphoid", gene_set, description = "T-cell-mediated cytotoxicity transcription", set_size, expressed, coverage, n_evaluation = NA_integer_, n_predictor = NA_integer_, retained)], fill = TRUE)
S4 <- ren(add_label(tidy(S4)), c(compartment = "Compartment", gene_set = "Gene set (MSigDB)", description = "Description",
          set_size = "Genes in set", expressed = "Expressed in discovery", coverage = "Coverage", n_evaluation = "Evaluation genes",
          n_predictor = "Predictor genes", retained = "Retained"))

S5 <- rd2("04_processed_data/programs/gene_partition.tsv", "04_processed_data/programs/gene_partition_lymphoid.tsv")
S5 <- ren(add_label(tidy(S5)), c(compartment = "Compartment", ens_clean = "Ensembl ID", gene_symbol = "Gene symbol",
          row_id = "Matrix row ID", expressed = "Expressed", role = "Role"))

frz <- readRDS("04_processed_data/frozen_panel/frozen_panel_v2.rds")
frz$programs <- c(frz$programs, readRDS("04_processed_data/frozen_panel/frozen_panel_LYM_v1.rds")$programs)
S6 <- rbindlist(lapply(names(frz$programs), function(pg) { P <- frz$programs[[pg]]
  data.table(program = pg, Compartment = P$compartment, Gene = P$panel$gene_symbol, `Matrix row ID` = P$panel$row_id,
             Coefficient = P$panel$coefficient, `Scaling mean` = P$panel$pred_mean, `Scaling SD` = P$panel$pred_sd,
             Intercept = P$intercept, `Evaluation genes` = nrow(P$evaluation_genes)) }))
S6 <- add_label(tidy(S6, 4)); S6[abs(Intercept) < 1e-8, Intercept := 0]

MLAB <- c(panel_refit = "Compact panel (relaxed refit; frozen)", panel_shrunk = "Compact panel (penalised coefficients)",
          mean = "Discovery mean (baseline)", equal_weight = "Equal-weight top-4 score", conventional = "Conventional markers",
          larger_model = "Unrestricted elastic net")
S7 <- rd2("06_results/panel/cv_performance_discovery_v2.csv", "06_results/panel/cv_performance_discovery_lymphoid.csv")
S7[, method := ifelse(method %in% names(MLAB), MLAB[method], method)]
S7 <- ren(add_label(tidy(S7)), c(n_donors = "Donors", method = "Model", MAE = "MAE", MAE_rep_min = "MAE, lowest repeat",
          MAE_rep_max = "MAE, highest repeat", RMSE = "RMSE", spearman = "Spearman ρ", cal_slope = "Calibration slope",
          cal_intercept = "Calibration intercept", MAE_vs_mean = "MAE ratio vs discovery mean"))
S7b <- rd2("06_results/panel/selection_stability_v2.csv", "06_results/panel/selection_stability_lymphoid.csv")
S7b <- ren(add_label(tidy(S7b)), c(row_id = "Matrix row ID", gene_symbol = "Gene symbol", selection_freq = "Selection frequency"))
in_panel <- S6[, paste(`Program code`, Gene)]
S7b[, `In frozen panel` := ifelse(paste(`Program code`, `Gene symbol`) %in% in_panel, "Yes", "No")]
setorder(S7b, `Program code`, -`Selection frequency`)
S7b[, .o := match(`Program code`, ORD)]; setorder(S7b, .o, -`Selection frequency`); S7b[, .o := NULL]

S8 <- rd2("06_results/cross_cohort_meta/per_cohort_results_long.csv", "06_results/validation_lymphoid/per_cohort_results_long_lymphoid.csv")
S8[, cohort := sub("_KUL$", "", cohort)]
S8[, scaling := fifelse(scaling == "frozen", "Frozen discovery scaling (primary)", "Within-cohort standardisation (secondary)")]
S8 <- ren(add_label(tidy(S8)), c(cohort = "Cohort", scaling = "Scaling", n_donors = "Donors", MAE_ratio = "MAE ratio",
          ratio_lo = "95% CI, lower", ratio_hi = "95% CI, upper", spearman = "Spearman ρ", status = "Status"))
S8b <- rd2("06_results/cross_cohort_meta/meta_analysis_results.csv", "06_results/validation_lymphoid/meta_analysis_results_lymphoid.csv")
S8b[is.na(analysis), analysis := "main"]
S8b[, analysis := fifelse(analysis == "main", "Main (independent cohorts)", "Including GSE144735 (sensitivity)")]
S8b[, scaling := fifelse(scaling == "frozen", "Frozen", "Within-cohort")]
S8b[, cohorts := gsub(";", ", ", cohorts)]
S8b <- ren(add_label(tidy(S8b)), c(analysis = "Analysis", scaling = "Scaling", k_cohorts = "Cohorts (k)", cohorts = "Cohorts",
          total_donors = "Donors", n_supported = "Cohorts supported", pooled_MAE_ratio = "Pooled MAE ratio",
          ratio_lo = "95% CI, lower", ratio_hi = "95% CI, upper", tau2_ratio = "τ² (log ratio)", I2_ratio = "I² (%), ratio",
          pooled_rho = "Pooled ρ", rho_lo = "ρ 95% CI, lower", rho_hi = "ρ 95% CI, upper", I2_rho = "I² (%), ρ",
          meta_status = "Status"))

S9 <- rd("06_results/sensitivity_discovery/discovery_sensitivity_results.csv")
S9 <- ren(add_label(tidy(S9)), c(scenario = "Scenario", n_donors = "Donors", MAE_ratio = "MAE ratio vs mean",
          ratio_range = "Range over repeats", spearman = "Spearman ρ", frozen_genes_reselected = "Frozen genes reselected",
          selected = "Genes selected", note = "Note", primary = "Primary-analysis MAE ratio"))

S10a <- add_label(tidy(spA)); S10a <- ren(S10a, c(cohort = "Cohort", patient = "Patient", effect = "Criterion A coefficient", se = "Bootstrap SE"))
S10b <- rd2("06_results/spatial_support_combined/criterionB_patients_combined.csv", "06_results/spatial_support_lymphoid_combined/criterionB_patients_combined.csv")
S10b <- ren(add_label(tidy(S10b)), c(cohort = "Cohort", patient = "Patient", diff = "Criterion B difference"))
S10c <- rd2("06_results/spatial_support_combined/spatial_support_summary_combined.csv", "06_results/spatial_support_lymphoid_combined/spatial_support_summary_combined.csv")
S10c <- ren(add_label(tidy(S10c)), c(analysis = "Analysis", n_patients = "Patients", cohorts = "Cohorts", pooled_effect = "Pooled criterion A",
          lo = "95% CI, lower", hi = "95% CI, upper", I2 = "I² (%)", A_positive = "Patients positive (A)", A_met = "Criterion A met",
          B_positive = "Patients positive (B)", B_met = "Criterion B met", spatial_label = "Spatial label"))
cv1 <- rd2("06_results/spatial_support/spatial_gene_coverage.csv", "06_results/spatial_support_lymphoid/spatial_gene_coverage.csv")[, cohort := "Valdeolivas"]
cv2 <- rd("06_results/spatial_support_extra/spatial_gene_coverage_extra.csv")
cv2[, cohort := fifelse(grepl("^GSM7089", section), "GSE226997", "GSE280318")]
cv3 <- rd("06_results/spatial_support_lymphoid/spatial_gene_coverage_extra.csv")
S10d <- rbind(cv1, cv2, cv3, fill = TRUE)
setcolorder(S10d, c("cohort", "section", "program"))
S10d <- ren(add_label(tidy(S10d)), c(cohort = "Cohort", section = "Section", panel_present = "Panel genes present",
          panel_genes = "Panel genes", missing_panel = "Missing panel genes", eval_coverage = "Evaluation-gene coverage"))

S11 <- rd("08_tables/Table_program_evidence_8.csv")
S11 <- ren(add_label(tidy(S11)), c(description = "Description", compartment = "Compartment", panel = "Panel genes",
          sc_source = "Single-cell evidence source", frozen_MAE_ratio_CI = "Frozen MAE ratio (95% CI)",
          sc_independent = "Independent single-cell support", within_scaling_supported = "Supported under within-cohort scaling",
          spatial_label = "Spatial label", spatial_effect_CI = "Pooled criterion A (95% CI)", evidence_category = "Evidence category",
          lopo_min_effect = "Leave-one-patient-out minimum effect", lopo_min_lo = "LOPO minimum lower CI",
          A_robust_all_LOPO = "Criterion A robust to every LOPO"))
S12 <- rd("08_tables/Table_chip_assay_proposal_8.csv")
S12 <- ren(tidy(S12), c(gene = "Gene", coefficient_sign = "Coefficient sign", compartment = "Compartment", chip_cell_model = "Chip cell model",
          program = "Program", single_cell_evidence = "Single-cell evidence", spatial_evidence = "Spatial evidence",
          evidence_category = "Evidence category", proposed_future_assay = "Proposed future assay", interpretation_limits = "Interpretation limits"))

S13 <- rd("06_results/GSE200997_annotation/cluster_marker_assignment.csv")
S13[, assigned := sub("^ctx_", "Context: ", assigned)]
setnames(S13, names(S13), sub("^ctx_", "", names(S13)))
S13 <- ren(tidy(S13, 2), c(cluster = "Seurat cluster", n_cells = "Cells", assigned = "Assigned compartment", margin = "Margin (best minus second)"))

## ---- table catalogue ------------------------------------------------------------------------------------
CAT <- list(
 list("S1a", "S1a_cohort_eligibility", "Dataset audit and eligibility", "One row per dataset audited before analysis. Roles were fixed before expression data were analysed.", S1a, S1a_note),
 list("S1b", "S1b_compartment_coverage", "Eligible donors per cohort and compartment", "Single-cell cohorts; source data for Fig. 1B.", S1b, S1b_note),
 list("S1c", "S1c_validation_patients", "Independent validation patients by program", "Patients contributing to the frozen-panel validation of each program.", S1c, S1c_note),
 list("S2a", "S2a_cohort_provenance", "Cohort provenance and independence", "Recruitment site, depositing institution and publication for every cohort.", S2a, S2a_note),
 list("S2b", "S2b_discovery_donors", "Discovery donors: identifiers, clinical characteristics and cell counts", "GSE178341, one row per patient.", S2b, S2b_note),
 list("S3",  "S3_software", "Software versions", "R and package versions used by the pipeline.", S3, S3_note),
 list("S4",  "S4_program_definitions", "Program definitions", "All ten candidate programs, gene-set source and coverage in the discovery cohort; retention required at least 80% coverage, 10 evaluation and 5 predictor genes. STR_ECM (stromal extracellular-matrix organisation) and LYM_CYTOTOX (CD8 T-cell-mediated cytotoxicity) failed and were not modelled.", S4, ""),
 list("S5",  "S5_gene_partition", "Gene partition", "Assignment of each expressed program gene to the evaluation (target) or predictor (panel-eligible) pool; the pools are disjoint across all programs.", S5, ""),
 list("S6",  "S6_frozen_panels", "Frozen panels", "Genes, relaxed-refit coefficients and discovery scaling values; panel score = intercept + Σ coefficient × (log-CPM − mean)/SD. Both panel files are SHA-256 verified before use.", S6, ""),
 list("S7",  "S7_discovery_CV", "Discovery cross-validation", "Donor-grouped five-fold cross-validation repeated ten times. The repeat range is the source of the bars in Supplementary Fig. S2A.", S7, ""),
 list("S7b", "S7b_selection_stability", "Gene selection stability", "Frequency with which each predictor gene was selected across 50 cross-validation fits (source data for Supplementary Fig. S2B).", S7b, ""),
 list("S8",  "S8_validation_by_cohort", "Frozen-panel validation by cohort", "MAE ratio against the discovery-mean baseline with 95% donor-bootstrap CI (2,000 resamples). Supported: ratio and upper CI below 1.", S8, ""),
 list("S8b", "S8b_meta_analysis", "Random-effects meta-analysis", "DerSimonian–Laird pooling over cohorts with at least 10 donors; main analysis excludes GSE144735.", S8b, ""),
 list("S9",  "S9_sensitivity", "Discovery sensitivity analyses", "Cross-validation repeated under alternative cell thresholds, partitions, normalisation, chemistry and sorting rules; the frozen panels were not modified.", S9, ""),
 list("S10a","S10a_spatial_criterionA", "Spatial criterion A by patient", "Composition-adjusted standardised panel coefficient; SE from a spatial block bootstrap (8 × 8 array blocks, 1,000 resamples).", S10a, ""),
 list("S10b","S10b_spatial_criterionB", "Spatial criterion B by patient", "Mean target difference between expected-high and comparison regions.", S10b, ""),
 list("S10c","S10c_spatial_summary", "Spatial support summary", "Pooled criterion A and patient consistency for the main analysis and by cohort group (source data for Supplementary Fig. S4).", S10c, ""),
 list("S10d","S10d_spatial_coverage", "Spatial gene coverage", "Panel-gene presence and evaluation-gene coverage per section and program.", S10d, ""),
 list("S11", "S11_evidence_classification", "Evidence classification", "Evidence category per program under the pre-specified rule (source data for Fig. 4).", S11, ""),
 list("S12", "S12_chip_assay_table", "Proposed chip measurements", "All 30 panel transcripts with evidence, proposed assay and interpretation limits (extends Table 1).", S12, ""),
 list("S13", "S13_GSE200997_annotation", "GSE200997 cluster annotation", "Pre-specified marker rule applied to Seurat clusters: mean marker-set score per cluster, assigned compartment and margin between the best and second-best set.", S13, ""))

REFS <- c(
 "Cable DM et al. (2022) Robust decomposition of cell type mixtures in spatial transcriptomics. Nat. Biotechnol., 40, 517–526. https://doi.org/10.1038/s41587-021-00830-w",
 "Chen Y, Chen L, Lun ATL, Baldoni PL and Smyth GK (2025) edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets. Nucleic Acids Res., 53, gkaf018. https://doi.org/10.1093/nar/gkaf018",
 "DerSimonian R and Laird N (1986) Meta-analysis in clinical trials. Control. Clin. Trials, 7, 177–188. https://doi.org/10.1016/0197-2456(86)90046-2",
 "Friedman J, Hastie T and Tibshirani R (2010) Regularization paths for generalized linear models via coordinate descent. J. Stat. Softw., 33, 1–22. https://doi.org/10.18637/jss.v033.i01",
 "Hao Y et al. (2024) Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nat. Biotechnol. https://doi.org/10.1038/s41587-023-01767-y",
 "Higgins JPT and Thompson SG (2002) Quantifying heterogeneity in a meta-analysis. Stat. Med., 21, 1539–1558. https://doi.org/10.1002/sim.1186",
 "Joanito I et al. (2022) Single-cell and bulk transcriptome sequencing identifies two epithelial tumor cell states and refines the consensus molecular classification of colorectal cancer. Nat. Genet., 54, 963–975. https://doi.org/10.1038/s41588-022-01100-4",
 "Khaliq AM et al. (2022) Refining colorectal cancer classification and clinical stratification through a single-cell atlas. Genome Biol., 23. https://doi.org/10.1186/s13059-022-02677-z",
 "Künsch HR (1989) The jackknife and the bootstrap for general stationary observations. Ann. Stat., 17, 1217–1241. https://doi.org/10.1214/aos/1176347265",
 "Lee HO et al. (2020) Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer. Nat. Genet., 52, 594–603. https://doi.org/10.1038/s41588-020-0636-z",
 "Liberzon A et al. (2015) The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst., 1, 417–425. https://doi.org/10.1016/j.cels.2015.12.004",
 "Meinshausen N (2007) Relaxed lasso. Comput. Stat. Data Anal., 52, 374–393. https://doi.org/10.1016/j.csda.2006.12.019",
 "Oliveira MF et al. (2025) High-definition spatial transcriptomic profiling of immune cell populations in colorectal cancer. Nat. Genet. https://doi.org/10.1038/s41588-025-02193-3",
 "Park SS, Lee Y-K, Choi YW, et al. (2024) Cellular senescence is associated with the spatial evolution toward a higher metastatic phenotype in colorectal cancer. Cell Rep. PMID: 38446659",
 "Pelka K et al. (2021) Spatially organized multicellular immune hubs in human colorectal cancer. Cell, 184, 4734–4752.e20. https://doi.org/10.1016/j.cell.2021.08.003",
 "Uhlitz F et al. (2021) Mitogen-activated protein kinase activity drives cell trajectories in colorectal cancer. EMBO Mol. Med., 13, e14123. https://doi.org/10.15252/emmm.202114123",
 "Valdeolivas A et al. (2024) Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics. NPJ Precis. Oncol., 8, 10. https://doi.org/10.1038/s41698-023-00488-4",
 "Zhang L et al. (2020) Single-cell analyses inform mechanisms of myeloid-targeted therapies in colon cancer. Cell, 181, 442–459.e29. https://doi.org/10.1016/j.cell.2020.03.048",
 "Zou H and Hastie T (2005) Regularization and variable selection via the elastic net. J. R. Stat. Soc. Series B Stat. Methodol., 67, 301–320. https://doi.org/10.1111/j.1467-9868.2005.00503.x")

ABBR <- c("CI: confidence interval", "CPM: counts per million", "FFPE: formalin-fixed paraffin-embedded",
          "LOPO: leave one patient out", "MAE: mean absolute error", "MMR: mismatch repair (MMRd deficient, MMRp proficient)",
          "SD/SE: standard deviation / standard error",
          "Criterion A: composition-adjusted association between panel score and program target in spots",
          "Criterion B: higher program target in expected-high than comparison regions",
          "MAE ratio: panel MAE divided by the MAE of predicting the discovery mean; below 1 favours the panel")

## ---- write workbook --------------------------------------------------------------------------------------
wb <- createWorkbook(creator = AUTH, title = TITLE)
st_t <- createStyle(fontSize = 12, textDecoration = "bold")
st_l <- createStyle(wrapText = TRUE, valign = "top", fontColour = "#333333")
st_h <- createStyle(textDecoration = "bold", fgFill = "#DDE6F0", border = "bottom", wrapText = TRUE, valign = "top")
st_n <- createStyle(fontColour = "#555555", textDecoration = "italic", wrapText = TRUE, valign = "top")
modifyBaseFont(wb, fontSize = 10, fontName = "Arial")

addWorksheet(wb, "Contents")
writeData(wb, "Contents", "Supplementary Tables", startRow = 1); addStyle(wb, "Contents", st_t, rows = 1, cols = 1)
writeData(wb, "Contents", TITLE, startRow = 2); writeData(wb, "Contents", AUTH, startRow = 3)
contents <- rbindlist(lapply(CAT, function(z) data.table(Table = paste("Supplementary Table", z[[1]]), Sheet = z[[2]],
                                                         Title = z[[3]], Rows = nrow(z[[5]]))))
contents <- rbind(contents, data.table(Table = "References", Sheet = "References", Title = "Works cited in these tables", Rows = length(REFS)))
writeData(wb, "Contents", contents, startRow = 5, headerStyle = st_h)
r0 <- 5 + nrow(contents) + 2
writeData(wb, "Contents", "Abbreviations and definitions", startRow = r0); addStyle(wb, "Contents", createStyle(textDecoration = "bold"), rows = r0, cols = 1)
writeData(wb, "Contents", ABBR, startRow = r0 + 1)
setColWidths(wb, "Contents", cols = 1:4, widths = c(28, 30, 70, 8))
for (i in seq_len(nrow(contents)))
  writeFormula(wb, "Contents", startRow = 5 + i, startCol = 2,
               x = sub("^=", "", makeHyperlinkString(sheet = contents$Sheet[i], row = 1, col = 1, text = contents$Sheet[i])))

for (z in CAT) {
  sh <- z[[2]]; d <- z[[5]]; addWorksheet(wb, sh)
  writeData(wb, sh, paste0("Supplementary Table ", z[[1]], ". ", z[[3]]), startRow = 1); addStyle(wb, sh, st_t, rows = 1, cols = 1)
  writeData(wb, sh, z[[4]], startRow = 2); addStyle(wb, sh, st_l, rows = 2, cols = 1)
  if (nzchar(z[[6]])) { writeData(wb, sh, z[[6]], startRow = 3); addStyle(wb, sh, st_n, rows = 3, cols = 1) }
  writeData(wb, sh, d, startRow = 5, headerStyle = st_h, keepNA = FALSE)
  freezePane(wb, sh, firstActiveRow = 6)
  w <- pmin(60, pmax(10, vapply(names(d), function(n) max(nchar(n) * 0.9, quantile(nchar(as.character(d[[n]])), 0.9, na.rm = TRUE), na.rm = TRUE), numeric(1)) + 2))
  setColWidths(wb, sh, cols = seq_along(d), widths = w)
  mergeCells(wb, sh, cols = 1:max(6, min(ncol(d), 10)), rows = 2); mergeCells(wb, sh, cols = 1:max(6, min(ncol(d), 10)), rows = 3)
  setRowHeights(wb, sh, rows = 2:3, heights = c(30, if (nzchar(z[[6]])) 45 else 15))
}
addWorksheet(wb, "References")
writeData(wb, "References", "References", startRow = 1); addStyle(wb, "References", st_t, rows = 1, cols = 1)
writeData(wb, "References", REFS, startRow = 3); setColWidths(wb, "References", cols = 1, widths = 180)
saveWorkbook(wb, OUT, overwrite = TRUE)

cat("Written:", OUT, "\n"); print(contents)
cat("\nValidation patients:", S1c_note, "\n")
