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
mask_mb6 = dba.mask(db_counts, DBA_CONDITION, "MB6")
mask_ctrl = dba.mask(db_counts, DBA_CONDITION, "Ctrl")
# ensure statistics will be of the form MB6 - Ctrl
db_contrast = dba.contrast(db_counts, mask_mb6, mask_ctrl, "MB6", "Ctrl")
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
fwrite(
    sigsites[abs(Fold) >= 1.5, .SD],
    paste0(ARGS$outdir, "/diff-bound-sites.significant-threshold.tsv"),
    sep = "\t"
)
fwrite(
    sigsites[, .SD, .SDcols = c("seqnames", "start", "end")],
    paste0(ARGS$outdir, "/diff-bound-sites.significant.bed"),
    sep = "\t",
    col.names = FALSE
)
fwrite(
    sigsites[abs(Fold) >= 1.5, .SD, .SDcols = c("seqnames", "start", "end")],
    paste0(ARGS$outdir, "/diff-bound-sites.significant-threshold.bed"),
    sep = "\t",
    col.names = FALSE
)

# save all sites
sites_gr = dba.report(db_analysis, th = 1)
sites = as.data.table(sites_gr)

fwrite(
    sites,
    paste0(ARGS$outdir, "/diff-bound-sites.all.tsv"),
    sep = "\t"
)

# calculate the number of significant sites at a certain fold-change threshold
deltas = seq(0, 3, 0.05)
sizes = data.table(
    Threshold = deltas,
    Increased = sapply(
        deltas,
        function(d) {sigsites[Fold >= d, .N]}
    ),
    Decreased = sapply(
        deltas,
        function(d) {sigsites[-Fold >= d, .N]}
    )
)
sizes[, "% Increased" := Increased / (Increased + Decreased)]
sizes[, "% Decreased" := Decreased / (Increased + Decreased)]
sizes_counts = melt(
    sizes,
    id.vars = "Threshold",
    measure.vars = c("Increased", "Decreased"),
    variable.name = "Direction",
    value.name = "Count"
)
sizes_prop = melt(
    sizes,
    id.vars = "Threshold",
    measure.vars = c("% Increased", "% Decreased"),
    variable.name = "Direction",
    value.name = "Percentage"
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

gg <- (
    ggplot(data = sizes_counts)
    + geom_col(aes(x = Threshold, y = Count, fill = Direction))
    + labs(x = "log(FC) Threshold", y = "Number of significant sites")
)
ggsave(
    paste0(ARGS$outdir, "/sigsites-fc-threshold.count.png"),
    height = 12,
    width = 20,
    units = "cm"
)
gg <- (
    ggplot(data = sizes_prop)
    + geom_col(aes(x = Threshold, y = Percentage, fill = Direction))
    + labs(x = "log(FC) Threshold", y = "Fraction of significant sites")
    # + scale_x_discrete(breaks = seq(0, 3, 0.5), labels = seq(0, 3, 0.5))
    # + scale_fill_continuous(labels = c("Increased", "Decreased"))
)
ggsave(
    paste0(ARGS$outdir, "/sigsites-fc-threshold.fraction.png"),
    height = 12,
    width = 20,
    units = "cm"
)

cat("Done\n")
