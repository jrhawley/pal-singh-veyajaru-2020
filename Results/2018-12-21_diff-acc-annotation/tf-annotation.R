# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
# read CisBP information
cisbp = fread("../../Data/External/CisBP/TF_Information_all_motifs_plus.txt")
# make colnames unique (this has 'DBID' as the name for 2 different columns)
colnames(cisbp)[14] = "DBID_Named"

# read differentially accessible loci from previous work
dars = fread("../2018-11-29_differential-accessibility/DiffBind/Filtered/EdgeR/diff-bound-sites.significant-threshold.annotated.tsv")

# ==============================================================================
# Analysis
# ==============================================================================
# clean up Ensembl IDs in `dars` since they contain versions, whereas CisBP
#   does not (not doing so makes it hard to merge IDs)
# remove '.#' in Ensembl IDs
dars[, CleanedIDs := gsub("\\.\\d+", "", NearestGeneEnsemblID)]
# for loci with multiple Ensembl IDs, split and copy them
#   (the most annotations that any locus has is 3 Ensembl IDs)
dars_duplicate = dars[grep(";", CleanedIDs)]
dars_triplicate = dars[grep("^ENSG.*; ENSG.*; ENSG.*", CleanedIDs)]
#   take first ID in original `dars` object
dars[, CleanedIDs := gsub("; .*", "", CleanedIDs)]
#   take second ID in original `dars` object
dars_duplicate[, CleanedIDs := gsub("^ENSG.*; ", "", CleanedIDs)]
dars_triplicate[, CleanedIDs := gsub("^ENSG\\d+; ENSG\\d+; ", "", CleanedIDs)]

dars_all = rbindlist(list(dars, dars_duplicate, dars_triplicate))

# find TFs from annotated list
tfs = merge(
    dars_all,
    cisbp[, .(DBID, TF_Name, Family_Name)],
    by.x = "CleanedIDs",
    by.y = "DBID"
)

# remove duplicates and keep unqiue loci
unique_idx = c(1, 3, 11, 21, 31, 32)

# save unique TFs
fwrite(
    tfs[unique_idx],
    "diff-bound-sites.TFs.tsv",
    col.names = TRUE,
    sep = "\t"
)
