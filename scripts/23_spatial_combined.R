
# 23_spatial_combined.R  (protocol v1.5) - pooled spatial support over all eligible patients from all spatial cohorts
suppressPackageStartupMessages(library(data.table))
od <- "06_results/spatial_support_combined"; dir.create(od, recursive = TRUE, showWarnings = FALSE)
va <- fread("06_results/spatial_support/criterionA_by_patient.csv")[, .(cohort = "Valdeolivas", patient, program, effect, se)]
vb <- fread("06_results/spatial_support/criterionB_by_patient.csv")[is.finite(diff), .(cohort = "Valdeolivas", patient, program, diff)]
ea <- fread("06_results/spatial_support_extra/criterionA_by_section_extra.csv")[role == "main", .(effect = mean(effect), se = sqrt(mean(se^2))), by = .(cohort, patient, program)]
eb <- fread("06_results/spatial_support_extra/criterionB_by_section_extra.csv")[role == "main", .(diff = mean(diff)), by = .(cohort, patient, program)]
A <- rbind(va, ea, use.names = TRUE); B <- rbind(vb, eb, use.names = TRUE)
dl <- function(est, se) { w <- 1 / se^2; fe <- sum(w * est) / sum(w); k <- length(est); Q <- sum(w * (est - fe)^2)
  t2 <- if (k > 1) max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))) else 0; wr <- 1 / (se^2 + t2); m <- sum(wr * est) / sum(wr); s <- sqrt(1 / sum(wr))
  c(est = m, lo = m - 1.96 * s, hi = m + 1.96 * s, I2 = if (k > 1 && Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0) }
summ <- function(a, b, label) rbindlist(fill = TRUE, lapply(sort(unique(a$program)), function(pg) {
  x <- a[program == pg & is.finite(se) & se > 0]; y <- b[program == pg]; np <- nrow(x)
  if (np < 3) return(data.table(analysis = label, program = pg, n_patients = np, spatial_label = "Not evaluable"))
  m <- dl(x$effect, x$se); A_met <- m[["lo"]] > 0 && sum(x$effect > 0) >= ceiling(2 / 3 * np)
  B_met <- nrow(y) >= 3 && sum(y$diff > 0) >= ceiling(2 / 3 * nrow(y))
  data.table(analysis = label, program = pg, n_patients = np, cohorts = paste(unique(x$cohort), collapse = "+"),
    pooled_effect = round(m[["est"]], 3), lo = round(m[["lo"]], 3), hi = round(m[["hi"]], 3), I2 = round(m[["I2"]], 1),
    A_positive = paste0(sum(x$effect > 0), "/", np), A_met = A_met, B_positive = paste0(sum(y$diff > 0), "/", nrow(y)), B_met = B_met,
    spatial_label = if (A_met && B_met) "Spatial support" else if (A_met || B_met) "Partial spatial support" else "No spatial support") }))
S <- rbind(summ(A, B, "All spatial cohorts (main)"), summ(A[cohort == "Valdeolivas"], B[cohort == "Valdeolivas"], "Valdeolivas only"),
           summ(A[cohort != "Valdeolivas"], B[cohort != "Valdeolivas"], "Additional cohorts only"), fill = TRUE)
fwrite(S, file.path(od, "spatial_support_summary_combined.csv")); fwrite(A, file.path(od, "criterionA_patients_combined.csv")); fwrite(B, file.path(od, "criterionB_patients_combined.csv"))
# 19 Oct 2026: scripts 20, 21 and 40 read spatial_support_summary_main.csv, which no
# script wrote. It is now derived here from the main analysis rows, so the evidence
# categories are always regenerated from the current spatial results.
fwrite(S[analysis == "All spatial cohorts (main)"], file.path(od, "spatial_support_summary_main.csv"))
cat("Patients per cohort:\n"); print(A[, .(patients = uniqueN(patient)), by = cohort])
for (an in unique(S$analysis)) { cat("\n===", an, "===\n"); print(S[analysis == an, .(program, n_patients, pooled_effect, lo, hi, I2, A_positive, B_positive, spatial_label)]) }
 

dir.create("06_results/sessionInfo", recursive = TRUE, showWarnings = FALSE); writeLines(capture.output(sessionInfo()), "06_results/sessionInfo/sessionInfo_23.txt")
