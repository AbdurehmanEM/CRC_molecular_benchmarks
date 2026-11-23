###############################################################################
# 58_spatial_table_sensitivity.R - Reviewer points 11 and 14 (and code-review B2)
#
# Works from the published per-patient spatial estimates (S10a criterion A
# coefficients with bootstrap SE; S10b criterion B differences), so no spatial
# data need to be reloaded. It first reproduces S10c (gate), then reports:
#   (a) the two-thirds rule varied: 1/2, 2/3 (published), 3/4, all patients;
#   (b) a stricter patient-level rule for criterion A (lower 95% limit > 0);
#   (c) exclusion of GSE226997 (treatment status unknown)    -> point 14;
#   (d) leave-one-spatial-cohort-out;
#   (e) criterion B split by region rule (pathologist vs composition quartile) -> B2.
# For each scenario, the spatial label and the resulting evidence category
# (rule v1.3; single-cell support taken from S11) are recomputed.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision()

A <- read_supp("S10a_spatial_criterionA")
setnames(A, c("Program code", "Cohort", "Patient", "Criterion A coefficient", "Bootstrap SE"),
         c("program", "cohort", "patient", "coefA", "seA"))
B <- read_supp("S10b_spatial_criterionB")
setnames(B, c("Program code", "Cohort", "Patient", "Criterion B difference"), c("program", "cohort", "patient", "diffB"))
A[, `:=`(coefA = as.numeric(coefA), seA = as.numeric(seA))]; B[, diffB := as.numeric(diffB)]
A <- A[is.finite(coefA)]; B <- B[is.finite(diffB)]
S11 <- read_supp("S11_evidence_classification")
sc_support <- setNames(S11$`Independent single-cell support` == "Yes", S11$`Program code`)
S10c <- read_supp("S10c_spatial_summary")

rule_B <- function(program, cohort) ifelse(cohort == "Valdeolivas" & !grepl("^(END_|LYM_)", program),
                                           "pathologist annotation", "composition quartile")
B[, rule := rule_B(program, cohort)]

summarise_spatial <- function(a, b, thr = 2/3, strictA = FALSE) {
  progs <- union(a$program, b$program)
  rbindlist(lapply(progs, function(p) {
    aa <- a[program == p]; bb <- b[program == p]
    nA <- nrow(aa); nB <- nrow(bb)
    posA <- if (strictA) sum(aa$coefA - 1.959964 * aa$seA > 0) else sum(aa$coefA > 0)
    posB <- sum(bb$diffB > 0)
    pooled <- if (nA >= 1) dl_pool(aa$coefA, aa$seA) else NULL
    metA <- nA >= 3 && posA >= ceiling(thr * nA - 1e-9)
    metB <- nB >= 3 && posB >= ceiling(thr * nB - 1e-9)
    lab <- if (max(nA, nB) < 3) "Not evaluable" else if (metA && metB) "Spatial support" else if (metA || metB) "Partial spatial support" else "No spatial support"
    cat_ <- if (isTRUE(sc_support[[p]])) {
      if (lab == "Spatial support") "Independently supported with spatial evidence" else "Independently supported in single-cell data"
    } else S11[`Program code` == p, `Evidence category`]
    data.table(program = p, nA = nA, posA = posA, pooledA = if (!is.null(pooled)) pooled$est else NA_real_,
               loA = if (!is.null(pooled)) pooled$lo else NA_real_, hiA = if (!is.null(pooled)) pooled$hi else NA_real_,
               I2 = if (!is.null(pooled)) pooled$I2 else NA_real_, metA = metA, nB = nB, posB = posB, metB = metB,
               spatial_label = lab, evidence_category = cat_)
  }))
}

## ---- Gate: reproduce S10c main rows
main <- summarise_spatial(A, B)
pub <- S10c[Analysis == "All spatial cohorts (main)", .(program = `Program code`, pooled_pub = as.numeric(`Pooled criterion A`),
                                                        label_pub = `Spatial label`, A_pub = `Patients positive (A)`, B_pub = `Patients positive (B)`)]
chk <- merge(main[, .(program, pooledA = round(pooledA, 3), A = paste0(posA, "/", nA), B = paste0(posB, "/", nB), spatial_label)],
             pub, by = "program")
chk[, ok := abs(pooledA - pooled_pub) <= 0.002 & A == A_pub & B == B_pub & spatial_label == label_pub]
print(chk)
if (!all(chk$ok)) warning("S10c not reproduced exactly from S10a/S10b; inspect 58_spatial_gate.csv before using results.")
write_out(chk, "58_spatial_gate")

## ---- Scenarios
sc <- list()
for (t in list(c("1/2", 1/2), c("2/3 (published)", 2/3), c("3/4", 3/4), c("all patients", 1)))
  sc[[length(sc) + 1]] <- cbind(scenario = paste("threshold", t[1]), summarise_spatial(A, B, as.numeric(t[2])))
sc[[length(sc) + 1]] <- cbind(scenario = "criterion A: patient lower 95% limit > 0", summarise_spatial(A, B, strictA = TRUE))
sc[[length(sc) + 1]] <- cbind(scenario = "excluding GSE226997", summarise_spatial(A[cohort != "GSE226997"], B[cohort != "GSE226997"]))
for (co in unique(A$cohort))
  sc[[length(sc) + 1]] <- cbind(scenario = paste("without", co), summarise_spatial(A[cohort != co], B[cohort != co]))
sc[[length(sc) + 1]] <- cbind(scenario = "criterion B: pathologist-annotation patients only",
                              summarise_spatial(A, B[rule == "pathologist annotation"]))
sc[[length(sc) + 1]] <- cbind(scenario = "criterion B: composition-quartile patients only",
                              summarise_spatial(A, B[rule == "composition quartile"]))
sc <- rbindlist(sc)
num <- c("pooledA", "loA", "hiA", "I2"); sc[, (num) := lapply(.SD, round, 3), .SDcols = num]
pubcat <- setNames(S11$`Evidence category`, S11$`Program code`)
publab <- setNames(pub$label_pub, pub$program)
sc[, `:=`(label_changed = spatial_label != publab[program], category_changed = evidence_category != pubcat[program])]
write_out(sc, "58_spatial_table_sensitivity")
print(sc[label_changed | category_changed, .(scenario, program, A = paste0(posA, "/", nA), B = paste0(posB, "/", nB),
                                             spatial_label, evidence_category, category_changed)])
cat("\nOnly scenarios that change a label or category are printed above; the full table is in the CSV.\n")
