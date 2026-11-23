# install_packages_new_pc.R
# One-time installation of every R package used by the CRC biomarker pipeline
# (scripts 00-41) on a fresh Windows PC with R 4.6.x.
#
# Run in RStudio:  source("install_packages_new_pc.R")
# Takes ~30-45 min on an empty library. Safe to re-run: installed packages are skipped.
# If it stops, re-run it - it continues where it left off.

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 3600)
if (.Platform$OS.type == "windows") options(pkgType = "binary")   # pre-built Windows packages, no compiler needed
cat("R version:", R.version.string, "\nInstalling into:", .libPaths()[1], "\n\n")

# ---- 1. CRAN packages --------------------------------------------------------
cran <- c(
  "data.table", "Matrix", "digest", "R.utils", "jsonlite",   # data handling, checksums, downloads
  "ggplot2", "patchwork",                                    # figures
  "glmnet",                                                  # elastic net (scripts 10, 24, 32)
  "openxlsx",                                                # supplementary workbook (27, 41)
  "hdf5r",                                                   # reading .h5 files (06, 08, 17, 29; Seurat::Read10X_h5)
  "arrow",                                                   # Visium HD parquet positions (22)
  "Seurat",                                                  # GSE200997 annotation (14), Read10X_h5 (18, 19, 22, 36)
  "BiocManager", "remotes")
todo <- setdiff(cran, rownames(installed.packages()))
if (length(todo)) { cat("CRAN:", paste(todo, collapse = ", "), "\n"); install.packages(todo, dependencies = TRUE) }

# ---- 2. Bioconductor packages ------------------------------------------------
# BiocManager picks the Bioconductor release that matches this R version.
bioc <- c(
  "edgeR",                                                   # CPM / pseudobulk normalisation
  "GEOquery", "Biobase",                                     # GEO downloads and sample metadata (01, 04, 05)
  "msigdbr",                                                 # Hallmark / Reactome / GO gene sets (09, 28, 30)
  "rhdf5",                                                   # chunked h5 reading (17)
  "S4Vectors", "SummarizedExperiment", "SpatialExperiment",  # RCTD input objects (18, 22)
  "spacexr",                                                 # RCTD deconvolution - Bioconductor version (createRctd/runRctd)
  "scDblFinder", "SingleCellExperiment")                     # doublet detection (listed in the protocol)
todo <- setdiff(bioc, rownames(installed.packages()))
if (length(todo)) { cat("\nBioconductor", as.character(BiocManager::version()), ":", paste(todo, collapse = ", "), "\n")
  BiocManager::install(todo, update = FALSE, ask = FALSE) }

# ---- 3. Verify every package loads -------------------------------------------
all_pk <- c(cran, bioc)
chk <- data.frame(package = all_pk,
                  loads   = vapply(all_pk, function(p) suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE))), logical(1)))
chk$version <- vapply(seq_along(all_pk), function(i) if (chk$loads[i]) as.character(packageVersion(all_pk[i])) else "-", character(1))
print(chk, row.names = FALSE)

if (chk$loads[chk$package == "spacexr"] && !exists("createRctd", envir = asNamespace("spacexr")))
  warning("spacexr is installed but lacks createRctd(): it is the GitHub version. Run remove.packages('spacexr'); BiocManager::install('spacexr')")
if (chk$loads[chk$package == "hdf5r"]) {
  ok_h5 <- tryCatch({ f <- tempfile(fileext = ".h5"); h <- hdf5r::H5File$new(f, "w"); h$close_all(); TRUE }, error = function(e) FALSE)
  cat("\nhdf5r can create/read .h5 files:", ok_h5, "\n")
}

missing <- chk$package[!chk$loads]
if (length(missing)) {
  cat("\n*** NOT INSTALLED:", paste(missing, collapse = ", "), "***\n",
      "Re-run this script once. If they still fail, copy the error messages above and send them.\n")
} else cat("\nAll", length(all_pk), "packages installed and loading. Next: run 05_analysis_scripts/00_setup_and_versions.R from the project folder to record versions.\n")
