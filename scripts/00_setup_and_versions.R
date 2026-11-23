# 00_setup_and_versions.R
# Install packages once and record the exact software environment.
# Updated 19 Oct 2026: the list now covers every package the pipeline calls
# (hdf5r, rhdf5, arrow, jsonlite, SpatialExperiment, SummarizedExperiment,
# S4Vectors, Biobase and spacexr were used by scripts 05-41 but were neither
# installed here nor recorded in software_versions.csv).
#
# spacexr MUST be the Bioconductor release: scripts 18 and 22 call createRctd() /
# runRctd(), which do not exist in the GitHub (dmcable/spacexr) version
# (that version provides create.RCTD() / run.RCTD()).

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c("Matrix", "data.table", "ggplot2", "patchwork", "Seurat",
               "glmnet", "digest", "openxlsx", "R.utils",
               "hdf5r", "jsonlite", "arrow")
bioc_pkgs <- c("edgeR", "GEOquery", "Biobase", "scDblFinder", "SingleCellExperiment",
               "msigdbr", "rhdf5", "S4Vectors", "SummarizedExperiment",
               "SpatialExperiment", "spacexr")

for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in bioc_pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE)

if (!exists("createRctd", envir = asNamespace("spacexr")))
  stop("Installed spacexr lacks createRctd(): remove it and install the Bioconductor release ",
       "with BiocManager::install('spacexr').")

dir.create("06_results", showWarnings = FALSE)
pkgs <- c(cran_pkgs, bioc_pkgs)
src  <- function(p) { d <- packageDescription(p)
  if (!is.null(d$RemoteSha)) paste0("GitHub ", d$RemoteUsername, "/", d$RemoteRepo, "@", substr(d$RemoteSha, 1, 10))
  else if (!is.null(d$biocViews)) "Bioconductor" else "CRAN" }
versions <- data.frame(
  package = pkgs,
  version = vapply(pkgs, function(p) as.character(packageVersion(p)), character(1)),
  source  = vapply(pkgs, src, character(1)))
versions <- rbind(data.frame(package = "R",
                             version = paste(R.version$major, R.version$minor, sep = "."),
                             source = "R"),
                  data.frame(package = "Bioconductor",
                             version = as.character(BiocManager::version()), source = "BiocManager"),
                  versions)
write.csv(versions, "06_results/software_versions.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "06_results/sessionInfo_setup.txt")
print(versions, row.names = FALSE)
