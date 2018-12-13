# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("GenomicRanges"))

# ==============================================================================
# Functions
# ==============================================================================
#' Parse out attributes from the "attribute" column in the GFF3 files
#'
#' @param text Text to be parsed
#' @param att Attribute to return
#' @return Character corresponding to the input attribute
getAtt <- function(text, att) {
    # separate attributes
    separated = strsplit(text, ";")[[1]]
    # get attribute matching `att`
    pos = grep(att, separated)
    # find position of "="
    idx = regexpr("=", separated[pos])
    return(substring(separated[pos], idx + 1))
}

# ==============================================================================
# Data
# ==============================================================================
# read pre-parsed GENCODE file for hg19
gencode_dt = fread(
    "../../Data/External/Gencode/gencode.v19.annotation.genes.bed",
    header = FALSE,
    sep = "\t",
    select = c(1:3, 7),
    col.names = c("chr", "start", "end", "attributes")
)
# extract names out of `attributes` column
gencode_dt$name = sapply(
    gencode_dt$attributes,
    getAtt,
    att = "gene_name"
)
# extract Ensembl IDs out of `attributes` column
gencode_dt$id = sapply(
    gencode_dt$attributes,
    getAtt,
    att = "gene_id"
)

# read in significantly differentially accessible sites
sigsites_dt = fread(
    "DiffBind/Filtered/EdgeR/diff-bound-sites.significant.tsv",
    header = TRUE,
    sep = "\t"
)

# convert data.tables to GRanges objects for easier comparison
gencode = GRanges(
    seqnames = gencode_dt$chr,
    ranges = IRanges(
        start = gencode_dt$start + 1,
        end = gencode_dt$end
    ),
    name = gencode_dt$name,
    id = gencode_dt$id,
)
sigsites = GRanges(
    seqnames = sigsites_dt$seqnames,
    ranges = IRanges(
        start = sigsites_dt$start + 1,
        end = sigsites_dt$end
    )
)

# ==============================================================================
# Analysis
# ==============================================================================
# find the nearest sites, and their distances
hits = nearest(sigsites, gencode, select = "all")
dist = distanceToNearest(sigsites, gencode)

# record these findings in the data.table
sigsites_dt$NearestGeneName = sapply(
    1:sigsites_dt[, .N],
    function(i) {
        # get hits for the current query site, return all the names from GENCODE
        # and concatenate them as a ;-separated string
        paste(gencode_dt[subjectHits(hits[queryHits(hits) == i]), name], collapse = "; ")
    }
)
sigsites_dt$NearestGeneEnsemblID = sapply(
    1:sigsites_dt[, .N],
    function(i) {
        # get hits for the current query site, return all the names from GENCODE
        # and concatenate them as a ;-separated string
        paste(gencode_dt[subjectHits(hits[queryHits(hits) == i]), id], collapse = "; ")
    }
)
sigsites_dt$DistanceToNearest = mcols(dist)$distance

# save data
fwrite(
    sigsites_dt,
    "DiffBind/Filtered/EdgeR/diff-bound-sites.significant.annotated.tsv",
    col.names = TRUE,
    sep = "\t"
)
fwrite(
    sigsites_dt[abs(Fold) >= 1.5],
    "DiffBind/Filtered/EdgeR/diff-bound-sites.significant-threshold.annotated.tsv",
    col.names = TRUE,
    sep = "\t"
)
