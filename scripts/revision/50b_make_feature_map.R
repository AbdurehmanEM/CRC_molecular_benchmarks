###############################################################################
# 50b_make_feature_map.R - writes the discovery feature map (row_id -> symbol)
# needed by scripts 53 and 61, from the GSE178341 h5 file.
# Discovery pseudobulk rownames look like "ENSG00000243485.5_4". Two matching
# rules are tried and checked against the row_id/symbol pairs in S5:
#   (1) position: h5 feature i = pseudobulk row i (both have 43,113 rows);
#   (2) Ensembl ID: strip the "_n" suffix (and version if needed) and match.
# The rule that agrees best with S5 (>= 95%) is used; S5 symbols then take precedence.
# Needs rhdf5 (Bioconductor) or hdf5r.
###############################################################################
source("50_revision_config.R"); source("50_revision_helpers.R"); init_revision()

h5 <- list.files(GSE178341_H5_DIR, pattern = "\\.h5$", full.names = TRUE, recursive = TRUE)
if (!length(h5)) stop("No .h5 file under ", GSE178341_H5_DIR)
h5 <- h5[1]; msg("h5:", h5)

read_h5 <- function(path) {
  if (requireNamespace("rhdf5", quietly = TRUE)) return(as.character(rhdf5::h5read(h5, path)))
  if (requireNamespace("hdf5r", quietly = TRUE)) { f <- hdf5r::H5File$new(h5, mode = "r"); on.exit(f$close_all())
    return(as.character(f[[path]]$read())) }
  stop("Install rhdf5 (BiocManager::install('rhdf5')) or hdf5r")
}
ls_paths <- if (requireNamespace("rhdf5", quietly = TRUE)) rhdf5::h5ls(h5) else NULL
grp <- if (!is.null(ls_paths)) unique(ls_paths$group[ls_paths$name %in% c("id", "name") & grepl("features", ls_paths$group)])[1] else "/matrix/features"
if (is.na(grp)) grp <- "/matrix/features"
ids <- read_h5(paste0(grp, "/id")); nms <- read_h5(paste0(grp, "/name"))
msg("h5 features:", length(ids), "| e.g.", ids[1], nms[1])

disc <- readRDS(file.path(PB_DIR, PB_FILES[[DISCOVERY]]$Epithelial))
rid <- rownames(disc$counts %||% disc)
msg("pseudobulk rows:", length(rid), "| e.g.", rid[1])

partition <- load_partition()
truth <- unique(partition[, .(row_id, symbol)])
check <- function(map) { m <- map[truth$row_id]; mean(!is.na(m) & m == truth$symbol) }

maps <- list()
maps$exact_id <- setNames(nms[match(rid, ids)], rid)              # IDs may already carry the "_n" suffix
if (length(ids) == length(rid)) maps$position <- setNames(nms, rid)
rid_base <- sub("_\\d+$", "", rid)
m2 <- setNames(nms[match(rid_base, ids)], rid)
if (mean(is.na(m2)) > 0.5) m2 <- setNames(nms[match(sub("\\.\\d+$", "", rid_base), sub("\\.\\d+$", "", ids))], rid)
maps$ensembl_id <- m2
agree <- sapply(maps, check); print(round(agree, 4))
best <- names(which.max(agree))
# Renamed genes (GENCODE name in h5 vs current HGNC symbol in S5) make ~1-2% disagree; 95% is the bar.
if (max(agree) < 0.95) stop("No rule agrees with S5 for >= 95% of genes (best: ", best, " ", round(max(agree), 3),
                            "). Send the printed examples so the rule can be adapted.")
maps[[best]][truth$row_id] <- truth$symbol                       # S5 symbols take precedence
out <- data.table(row_id = rid, symbol = unname(maps[[best]]))
fwrite(out, DISCOVERY_FEATURES)
msg("Feature map written with rule '", best, "' (agreement with S5:", round(max(agree), 4), ") ->", DISCOVERY_FEATURES)
