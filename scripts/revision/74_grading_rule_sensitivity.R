###############################################################################
# 74_grading_rule_sensitivity.R - Protocol amendment v2.3, grading-rule sensitivity
#
# The evidence-category rule (protocol v1.3) was written after the single-cell and
# spatial results were known. This script re-derives every program's category under
# the published rule (checked against Supplementary Table S11) and under the three
# alternative rules fixed in protocol v2.3:
#   R1  single-cell support needs pooled frozen support AND >= 2 cohorts individually
#       Supported (single-cohort programs cannot meet it)
#   R2  within-cohort scaling (secondary) replaces frozen scaling throughout
#   R3a spatial support = criterion A alone
#   R3b spatial support = criteria A and B, each positive in >= 3/4 of patients
# Uses existing published results only (Supplementary Tables S8, S8b, S10c, S11).
# No new data are analysed; nothing is overwritten.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision()
suppressPackageStartupMessages({ library(readxl); library(data.table) })

WB <- file.path(PROJECT_ROOT, "10_submission", "Supplementary_Tables_V14_aem.xlsx")
stopifnot(file.exists(WB))
rd <- function(sheet) {                       # header row = row 5 in these sheets
  raw <- read_excel(WB, sheet = sheet, col_names = FALSE, .name_repair = "minimal")
  h <- which(as.character(raw[[1]]) == "Program code")[1]
  d <- as.data.table(read_excel(WB, sheet = sheet, skip = h - 1, .name_repair = "unique"))
  d[!is.na(`Program code`) & grepl("^[A-Z]{3}_", `Program code`)]
}
s8  <- rd("S8_validation_by_cohort")
s8b <- rd("S8b_meta_analysis")
s10 <- rd("S10c_spatial_summary")
s11 <- rd("S11_evidence_classification")

CAT <- c(top = "Independently supported with spatial evidence",
         sc  = "Independently supported in single-cell data",
         exp = "Exploratory", uns = "Unsupported")
frac <- function(x) { v <- as.numeric(strsplit(x, "/")[[1]]); v[1] / v[2] }

single_cell <- function(p, scal) {             # scal: "Frozen" or "Within"
  pat <- if (scal == "Frozen") "^Frozen" else "^Within"
  co <- s8[`Program code` == p & grepl(pat, Scaling) & Cohort != "GSE144735" &
           Status %in% c("Supported", "Partially supported", "Not supported")]
  m  <- s8b[`Program code` == p & grepl("^Main", Analysis) & grepl(pat, Scaling)]
  k  <- nrow(co)
  if (k >= 2) list(k = k, n_sup = sum(co$Status == "Supported"),
                   support = m$Status == "Supported", below1 = as.numeric(m$`Pooled MAE ratio`) < 1)
  else list(k = k, n_sup = sum(co$Status == "Supported"),
            support = k == 1 && co$Status == "Supported", below1 = k == 1 && as.numeric(co$`MAE ratio`) < 1)
}

category <- function(sc_support, explore, spatial_support) {
  if (sc_support && spatial_support) CAT[["top"]] else if (sc_support) CAT[["sc"]]
  else if (explore) CAT[["exp"]] else CAT[["uns"]]
}

out <- list()
for (p in s11$`Program code`) {
  fz <- single_cell(p, "Frozen"); wi <- single_cell(p, "Within")
  sp <- s10[`Program code` == p & grepl("main", Analysis)]
  spat_pub <- sp$`Spatial label` == "Spatial support"
  A_met <- sp$`Criterion A met` == "Yes"
  A_frac <- if (is.na(sp$`Patients positive (A)`)) NA else frac(sp$`Patients positive (A)`)
  B_frac <- if (is.na(sp$`Patients positive (B)`)) NA else frac(sp$`Patients positive (B)`)

  pub <- category(fz$support, fz$below1 || wi$support, spat_pub)
  r1_sc <- fz$k >= 2 && fz$support && fz$n_sup >= 2
  r1  <- category(r1_sc, fz$below1 || wi$support, spat_pub)
  r2  <- category(wi$support, wi$below1 || fz$support, spat_pub)
  r3a <- category(fz$support, fz$below1 || wi$support, isTRUE(A_met))
  r3b <- category(fz$support, fz$below1 || wi$support, isTRUE(A_frac >= 0.75 & B_frac >= 0.75))

  out[[p]] <- data.table(program = p, external_cohorts = fz$k, cohorts_supported_frozen = fz$n_sup,
    spatial_A = sp$`Patients positive (A)`, spatial_B = sp$`Patients positive (B)`,
    published_S11 = s11[`Program code` == p, `Evidence category`], rederived = pub,
    R1_two_cohorts = r1, R2_within_scaling = r2, R3a_criterionA_only = r3a, R3b_three_quarters = r3b)
}
res <- rbindlist(out)

## Gate: the published rule must reproduce Supplementary Table S11 exactly
bad <- res[published_S11 != rederived]
if (nrow(bad)) { print(bad[, .(program, published_S11, rederived)]); stop("Re-derived categories do not match S11") }
cat("Published categories reproduced for all", nrow(res), "programs\n")

alts <- c("R1_two_cohorts", "R2_within_scaling", "R3a_criterionA_only", "R3b_three_quarters")
chg <- sapply(alts, function(a) res[[a]] != res$published_S11)          # programs x rules
res[, changed_under := apply(chg, 1, function(v) if (any(v)) paste(alts[v], collapse = "; ") else "none")]
res[, rederived := NULL]
write_out(res, "74_grading_rule_sensitivity")

short <- function(x) c("Independently supported with spatial evidence" = "Top (sc + spatial)",
                       "Independently supported in single-cell data" = "Single-cell",
                       "Exploratory" = "Exploratory", "Unsupported" = "Unsupported")[x]
print(res[, .(program, published = short(published_S11), R1 = short(R1_two_cohorts), R2 = short(R2_within_scaling),
              R3a = short(R3a_criterionA_only), R3b = short(R3b_three_quarters))])
cat("\nPrograms whose category changes under at least one alternative rule:",
    sum(res$changed_under != "none"), "of", nrow(res), "\n")
print(res[changed_under != "none", .(program, changed_under)])
