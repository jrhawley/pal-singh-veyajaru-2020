    # ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("DiffBind"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Description"
    )
    PARSER$add_argument(
        "samplesheet",
        type = "character",
        help = "Path to input sample sheet to feed to DiffBind"
    )
    PARSER$add_argument(
        "-o", "--outdir",
        type = "character",
        help = "Path to output file directory. Default '.'",
        default = getwd()
    )
    PARSER$add_argument(
        "-m", "--method",
        type = "character",
        help = "DESeq2 or edgeR. Default 'DESeq2'",
        default = "DESeq2"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS = list(
        samplesheet = "diffbind-config-filtered.csv",
        outdir = "DiffBind/Filtered",
        method = "DESeq2"
    )
}

if (ARGS$method == "DESeq2") {
    ARGS$method = DBA_DESEQ2
} else if (ARGS$method == "edgeR") {
    ARGS$method = DBA_EDGER_GLM
}


# ==============================================================================
# Data
# ==============================================================================
cat("Loading data\n")
db = dba(sampleSheet=ARGS$samplesheet)
db$config$AnalysisMethod = ARGS$method

cat("Counting data\n")
db_counts = dba.count(db)

# ==============================================================================
# Analysis
# ==============================================================================
db_contrast = dba.contrast(db_counts, categories = DBA_CONDITION)
cat("Analyzing data\n")
db_analysis = dba.analyze(db_contrast)

# save significantly differential sites
sigsites_gr = dba.report(db_analysis)
sigsites = as.data.table(sigsites_gr)

fwrite(
    sigsites,
    paste0(ARGS$outdir, "/diff-bound-sites.significant.tsv"),
    sep = "\t"
)

# save all sites
sites_gr = dba.report(db_analysis, th = 1)
sites = as.data.table(sites_gr)

fwrite(
    sites,
    paste0(ARGS$outdir, "/diff-bound-sites.all.tsv"),
    sep = "\t"
)

# ==============================================================================
# Plots
# ==============================================================================
cat("Plotting data\n")
png(
    paste0(ARGS$outdir, "/correlation_heatmap_peakcallerscore.png"),
    width = 20,
    height = 20,
    units = "cm",
    res = 300
)
plot(db)
dev.off()

png(
    paste0(ARGS$outdir, "/correlation_heatmap_readcount.png"),
    width = 20,
    height = 20,
    units = "cm",
    res = 300
)
plot(db_counts)
dev.off()

png(
    paste0(ARGS$outdir, "/correlation_heatmap_diffbind.png"),
    width = 20,
    height = 20,
    units = "cm",
    res = 300
)
plot(db_analysis, contrast = 1)
dev.off()

png(
    paste0(ARGS$outdir, "/pca.png"),
    width = 20,
    height = 20,
    units = "cm",
    res = 300
)
dba.plotPCA(db_analysis, contrast = 1, label = DBA_CONDITION)
dev.off()

png(
    paste0(ARGS$outdir, "/ma.png"),
    width = 20,
    height = 20,
    units = "cm",
    res = 300
)
dba.plotMA(db_analysis, contrast = 1)
dev.off()

png(
    paste0(ARGS$outdir, "/volcano.png"),
    width = 20,
    height = 20,
    units = "cm",
    res = 300
)
dba.plotVolcano(db_analysis, contrast = 1)
dev.off()

png(
    paste0(ARGS$outdir, "/boxplot.png"),
    width = 20,
    height = 20,
    units = "cm",
    res = 300
)
dba.plotBox(db_analysis, contrast = 1, label = DBA_CONDITION)
dev.off()

png(
    paste0(ARGS$outdir, "/bindingaffinity_heatmap.png"),
    width = 20,
    height = 20,
    units = "cm",
    res = 300
)
dba.plotHeatmap(db_analysis, contrast = 1, correlations = FALSE)
dev.off()

gg <- (
    ggplot(data = sites)
    + geom_histogram(aes(x = p.value), binwidth = 0.01)
    + labs(x = "p-value", y = "Frequency")
)
ggsave(
    paste0(ARGS$outdir, "/p-value-histogram.png"),
    height = 12,
    width = 20,
    units = "cm"
)

cat("Done\n")
