# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("GenomicRanges"))

# ==============================================================================
# Functions
# ==============================================================================
#' Convert narrowPeak files to GRanges
#'
#' @param filename
#' @return GRanges
bed2gr <- function(filename, zipped = FALSE) {
    if (zipped) {
        filename = paste("zcat", filename)
    }
    dt = fread(
        filename,
        header = FALSE,
        sep = "\t",
        select = 1:3,
        col.names = c("chr", "start", "end")
    )
    return(GRanges(dt[, paste0(chr, ":", start, "-", end)]))
}

#' Get binary matrix of whether a list of query GRanges overlaps a subject GRange
#'
#' @param queries vector of file names to be compared against `subject`
#' @param subject file name of basis for comparison
#' @return data.frame of the consensus with the same number of columns as length
#' of the queries list
binarymat <- function(queries, subject) {
    # load subject
    subject_gr = bed2gr(subject)
    subject_df = data.frame(subject)
    # remove unnecesary columns
    subject_df$width = NULL
    subject_df$strand = NULL
    for (i in 1:length(queries)) {
        # load query
        query = bed2gr(queries[i])
        # find overlaps between query and subject
        hits = findOverlaps(query, subject_gr)
        # convert to data.frame
        hits_df = data.frame(hits)
        # in the data.frame, set as 1 if overlap, 0 if not
        subject_df[hits_df$subjectHits, names(queries)[i]] = 1
        subject_df[-hits_df$subjectHits, names(queries)[i]] = 0
    }
    return(subject_df)
}

# ==============================================================================
# Data
# ==============================================================================
consensus_sets = c(
    "all" = "Consensus/consensus.all.bed",
    "Ctrl-MB6" = "Consensus/consensus.Ctrl-MB6.bed",
    "Ctrl-MB6Pen" = "Consensus/consensus.Ctrl-MB6Pen.bed",
    "MB6-MB6Pen" = "Consensus/consensus.MB6-MB6Pen.bed"
)
comparisons = list(
    "all" = c(
        "Ctrl" = "Consensus/Ctrl.filtered.bed",
        "MB6" = "Consensus/MB6.filtered.bed",
        "MB6Pen" = "Consensus/MB6Pen.filtered.bed"
    ),
    "Ctrl-MB6" = c(
        "Ctrl" = "Consensus/Ctrl.filtered.bed",
        "MB6" = "Consensus/MB6.filtered.bed"
    ),
    "Ctrl-MB6Pen" = c(
        "Ctrl" = "Consensus/Ctrl.filtered.bed",
        "MB6Pen" = "Consensus/MB6Pen.filtered.bed"
    ),
    "MB6-MB6Pen" = c(
        "MB6" = "Consensus/MB6.filtered.bed",
        "MB6Pen" = "Consensus/MB6Pen.filtered.bed"
    )
)

# ==============================================================================
# Analysis
# ==============================================================================
for (i in 1:length(consensus_sets)){
    # get comparison
    comp = names(consensus_sets)[i]
    print(comp)
    # generate binary matrix for this comparison
    df = binarymat(comparisons[[comp]], consensus_sets[comp])
    # save to output file
    fwrite(
        df,
        paste("Consensus/consensus-mtx", comp,"tsv", sep = "."),
        col.names = TRUE,
        sep = "\t"
    )
}

