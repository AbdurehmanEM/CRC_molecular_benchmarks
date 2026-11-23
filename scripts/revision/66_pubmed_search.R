# 66_pubmed_search.R - documented PubMed search for the prior-art comparison (reviewer point 1)
# Runs the broad and focused queries through NCBI E-utilities, saves every focused record
# (PMID, year, journal, title, abstract, DOI) and flags records that may compete with this study.
# Output: 02_metadata/literature/pubmed_search_<date>/
#   pubmed_focused_records.csv   all focused records, sorted by flag score (highest first)
#   pubmed_search_log.txt        queries, date, counts (for the Methods / S17 note)
# Run from the project root:  source("05_analysis_scripts/revision/66_pubmed_search.R")
for (p in c("rentrez", "xml2", "data.table")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(rentrez); library(xml2); library(data.table) })

Q1 <- paste0('("colorectal neoplasms"[MeSH Terms] OR "colorectal cancer"[tiab] OR "colon cancer"[tiab] OR "rectal cancer"[tiab]) AND ',
             '("lab-on-a-chip devices"[MeSH Terms] OR "microphysiological systems"[MeSH Terms] OR "organ-on-a-chip"[tiab] OR "organ-on-chip"[tiab] OR ',
             '"tumor-on-a-chip"[tiab] OR "tumour-on-a-chip"[tiab] OR "cancer-on-a-chip"[tiab] OR "microphysiological system*"[tiab] OR microfluidic*[tiab])')
Q2 <- paste0(Q1, ' AND ("single-cell"[tiab] OR "single cell"[tiab] OR "spatial transcriptom*"[tiab] OR "RNA-seq"[tiab] OR ',
             '"RNA sequencing"[tiab] OR transcriptom*[tiab] OR biomarker*[tiab])')

today <- format(Sys.Date())
OUT <- file.path("02_metadata", "literature", paste0("pubmed_search_", today))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

s1 <- entrez_search(db = "pubmed", term = Q1, retmax = 0)
s2 <- entrez_search(db = "pubmed", term = Q2, retmax = 1000)
cat("Broad query:  ", s1$count, "records\nFocused query:", s2$count, "records\n")
ids <- s2$ids

## ---- fetch records in batches --------------------------------------------------------
txt <- function(node, xp) { v <- xml_text(xml_find_all(node, xp)); if (length(v)) paste(v, collapse = " ") else NA_character_ }
recs <- list()
for (b in split(ids, ceiling(seq_along(ids) / 100))) {
  x <- read_xml(entrez_fetch(db = "pubmed", id = b, rettype = "xml"))
  for (a in xml_find_all(x, ".//PubmedArticle")) {
    recs[[length(recs) + 1]] <- data.table(
      pmid     = txt(a, ".//MedlineCitation/PMID"),
      year     = substr(txt(a, ".//Article/Journal/JournalIssue/PubDate/Year | .//Article/Journal/JournalIssue/PubDate/MedlineDate"), 1, 4),
      journal  = txt(a, ".//Article/Journal/ISOAbbreviation"),
      first_author = txt(xml_find_first(a, ".//AuthorList/Author"), "./LastName"),
      title    = txt(a, ".//Article/ArticleTitle"),
      pub_type = txt(a, ".//PublicationTypeList/PublicationType"),
      doi      = txt(xml_find_first(a, "./PubmedData/ArticleIdList"), "./ArticleId[@IdType='doi']"),
      abstract = txt(a, ".//Abstract/AbstractText"))
  }
  Sys.sleep(0.4)   # stay within NCBI's request limit
}
d <- rbindlist(recs)

## ---- flags (to prioritise reading; the decision is made by reading) ------------------
has <- function(pat) grepl(pat, paste(d$title, d$abstract), ignore.case = TRUE, perl = TRUE)
d[, f_patient_sc_spatial := has("single[- ]cell|scRNA|snRNA|spatial transcriptom|visium|GSE[0-9]{4,}|TCGA|patient[- ]derived (scRNA|single)")]
d[, f_marker_panel       := has("panel|signature|gene set|marker genes|biomarker panel|benchmark|reference")]
d[, f_multicompartment   := has("fibroblast|stromal|endothelial|HUVEC|macrophage|monocyte|THP-1|immune")]
d[, f_compare_to_patient := has("compared (with|to) (patient|primary|human) (tissue|tumou?r)|recapitulat|fidelity|resembl|similar(ity)? to (patient|primary)")]
d[, is_review := grepl("Review", pub_type)]
d[, flag_score := f_patient_sc_spatial * 2 + f_marker_panel + f_multicompartment + f_compare_to_patient - is_review * 2]
setorder(d, -flag_score, -year)
d[, `:=`(decision = "", reason = "")]   # to be filled by the screener

fwrite(d, file.path(OUT, "pubmed_focused_records.csv"))
writeLines(c(paste("Search date:", today), "Database: PubMed (NCBI E-utilities via rentrez)",
             paste("Broad query:", Q1), paste("Broad count:", s1$count),
             paste("Focused query:", Q2), paste("Focused count:", s2$count),
             paste("Records saved:", nrow(d)), paste("Reviews among them:", sum(d$is_review))),
           file.path(OUT, "pubmed_search_log.txt"))

cat("\nSaved", nrow(d), "records to", file.path(OUT, "pubmed_focused_records.csv"), "\n")
cat("Reviews:", sum(d$is_review), " | flag score >= 3:", sum(d$flag_score >= 3), "\n\nTop 15 by flag score:\n")
print(d[1:min(15, .N), .(pmid, year, first_author, flag_score, title = substr(title, 1, 90))])
