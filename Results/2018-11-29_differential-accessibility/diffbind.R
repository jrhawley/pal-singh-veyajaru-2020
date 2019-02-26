# ==============================================================================
# Command Line
# ==============================================================================
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
        help = "DESeq2 or edgeR. Default 'edgeR'",
        default = "edgeR"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS = list(
        samplesheet = "config.csv",
        outdir = "DiffBind",
        method = "edgeR"
    )
}

# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("gtable"))
suppressMessages(library("grid"))
suppressMessages(library("DiffBind"))
suppressMessages(library("edgeR"))

if (ARGS$method == "DESeq2") {
    ARGS$method = DBA_DESEQ2
} else if (ARGS$method == "edgeR") {
    ARGS$method = DBA_EDGER_GLM
}

CHRS = paste0('chr', c(1:22, 'X', 'Y'))
FOLD_THRESH = 1.5
FDR_THRESH = 0.05

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
# mask_mb6 = dba.mask(db_counts, DBA_CONDITION, "MB6")
# mask_ctrl = dba.mask(db_counts, DBA_CONDITION, "Ctrl")
# mask_mb6pen = dba.mask(db_counts, DBA_CONDITION, "MB6Pen")
# ensure statistics will be of the form MB6 - Ctrl
db_contrast = dba.contrast(
    db_counts,
    categories = DBA_CONDITION,
    minMembers = 2
)

cat("Analyzing data\n")
db_analysis = dba.analyze(db_contrast, bReduceObjects = FALSE)

# save R objects for future loading, if needed
saveRDS(db_analysis, paste0(ARGS$outdir, "/diffbind-analysis.rds"))

# Extract DARs
# -------------------------------------
# get DARs from analyses
#   getting both normalized (for analysis)
#   and unnormalized (for plotting on shared axis)
sites_mb6_ctrl_norm = dba.report(
    db_analysis,
    contrast = 1,
    th = 1,
    bCounts = TRUE,
    bNormalized = TRUE
)
sites_mb6pen_ctrl_norm = dba.report(
    db_analysis,
    contrast = 2,
    th = 1,
    bCounts = TRUE,
    bNormalized = TRUE
)
sites_mb6_mb6pen_norm = dba.report(
    db_analysis,
    contrast = 3,
    th = 1,
    bCounts = TRUE,
    bNormalized = TRUE
)
sites_mb6_ctrl_unnorm = dba.report(
    db_analysis,
    contrast = 1,
    th = 1,
    bCounts = TRUE,
    bNormalized = FALSE
)
sites_mb6pen_ctrl_unnorm = dba.report(
    db_analysis,
    contrast = 2,
    th = 1,
    bCounts = TRUE,
    bNormalized = FALSE
)
sites_mb6_mb6pen_unnorm = dba.report(
    db_analysis,
    contrast = 3,
    th = 1,
    bCounts = TRUE,
    bNormalized = FALSE
)

# set order for GRanges objects
seqlevels(sites_mb6_ctrl_norm) = CHRS
seqlevels(sites_mb6pen_ctrl_norm) = CHRS
seqlevels(sites_mb6_mb6pen_norm) = CHRS
seqlevels(sites_mb6_ctrl_unnorm) = CHRS
seqlevels(sites_mb6pen_ctrl_unnorm) = CHRS
seqlevels(sites_mb6_mb6pen_unnorm) = CHRS
sites_mb6_ctrl_norm = sort(sites_mb6_ctrl_norm)
sites_mb6pen_ctrl_norm = sort(sites_mb6pen_ctrl_norm)
sites_mb6_mb6pen_norm = sort(sites_mb6_mb6pen_norm)
sites_mb6_ctrl_unnorm = sort(sites_mb6_ctrl_unnorm)
sites_mb6pen_ctrl_unnorm = sort(sites_mb6pen_ctrl_unnorm)
sites_mb6_mb6pen_unnorm = sort(sites_mb6_mb6pen_unnorm)

# convert to data.tables
mb6_ctrl = as.data.table(sites_mb6_ctrl_norm)
mb6pen_ctrl = as.data.table(sites_mb6pen_ctrl_norm)
mb6_mb6pen = as.data.table(sites_mb6_mb6pen_norm)

# convert start coordinates back to 0-indexed for writing
mb6_ctrl[, start := start - 1]
mb6pen_ctrl[, start := start - 1]
mb6_mb6pen[, start := start - 1]

# switch comparison direction so:
#   Fold = log2(MB6) - log2(Ctrl)
#   Fold = log2(MB6Pen) - log2(Ctrl)
#   Fold = log2(MB6) - log2(MB6Pen)
mb6_ctrl[, Fold := -Fold]
mb6pen_ctrl[, Fold := -Fold]

# Save tested regions
# -------------------------------------
agg_sites = list(mb6_ctrl, mb6pen_ctrl, mb6_mb6pen)
comp_names = c('mb6-ctrl', 'mb6pen-ctrl', 'mb6-mb6pen')
for (i in 1:3) {
    # save full dataset
    fwrite(
        agg_sites[[i]],
        paste0(ARGS$outdir, "/", comp_names[i], ".all.tsv"),
        sep = "\t"
    )
    # save significant sites that pass the fold-change threshold
    fwrite(
        agg_sites[[i]][FDR < FDR_THRESH & abs(Fold) >= FOLD_THRESH, .SD],
        paste0(ARGS$outdir, "/", comp_names[i], ".significant-threshold.tsv"),
        sep = "\t"
    )
    # save significant sites that pass the fold-change threshold as BED file
    fwrite(
        agg_sites[[i]][
            FDR < FDR_THRESH & abs(Fold) >= FOLD_THRESH,
            .SD,
            .SDcols = c("seqnames", "start", "end", "Fold", "FDR")
        ],
        paste0(ARGS$outdir, "/", comp_names[i], ".significant-threshold.bedGraph"),
        sep = "\t"
    )
}

# Calculate the number of significant sites at a certain fold-change threshold
# -------------------------------------
deltas = seq(0, 3, 0.05)
sizes = data.table(
    Threshold = deltas,
    Increased = sapply(
        deltas,
        function(d) {mb6_ctrl[Fold >= d & FDR < FDR_THRESH, .N]}
    ),
    Decreased = sapply(
        deltas,
        function(d) {mb6_ctrl[-Fold >= d & FDR < FDR_THRESH, .N]}
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

# save results
fwrite(
    sizes_counts,
    paste0(ARGS$outdir, "/mb6-ctrl.fdr-counts.tsv"),
    col.names = TRUE,
    sep = "\t"
)
fwrite(
    sizes_prop,
    paste0(ARGS$outdir, "/mb6-ctrl.fdr-proportions.tsv"),
    col.names = TRUE,
    sep = "\t"
)

# Get read counts for all samples (DARs and all sites)
# -------------------------------------
# get indices for DARs that are significantly more accessible in MB6 than Ctrl
sig_up_idx = mb6_ctrl[, which(FDR < FDR_THRESH & Fold >= FOLD_THRESH)]
# get indices for DARs that are significantly less accessible in MB6 than Ctrl
sig_dn_idx = mb6_ctrl[, which(FDR < FDR_THRESH & Fold <= -FOLD_THRESH)]

# get read counts for all sites in all samples
#   extract counts from reports
all_counts_ctrl = mcols(sites_mb6_ctrl_unnorm)[, 7:9]
all_counts_mb6 = mcols(sites_mb6_ctrl_unnorm)[, 10:12]
all_counts_mb6pen = mcols(sites_mb6pen_ctrl_unnorm)[, 10:11]
#   calculate log2(mean(counts)) for each peak
all_counts_ctrl = as.data.table(apply(
    all_counts_ctrl,
    1,
    function(x) log2(mean(x))
))
all_counts_mb6 = as.data.table(apply(
    all_counts_mb6,
    1,
    function(x) log2(mean(x))
))
all_counts_mb6pen = as.data.table(apply(
    all_counts_mb6pen,
    1,
    function(x) log2(mean(x))
))
#   add which sample and site they come from
loci = paste0(mb6_ctrl[, paste0(seqnames, ":", start, "-", end)])
all_counts_ctrl[, Condition := "Ctrl"]
all_counts_mb6[, Condition := "MB6"]
all_counts_mb6pen[, Condition := "MB6Pen"]
all_counts_ctrl[, Locus := loci]
all_counts_mb6[, Locus := loci]
all_counts_mb6pen[, Locus := loci]
colnames(all_counts_ctrl)[1] = "Log2MeanCount"
colnames(all_counts_mb6)[1] = "Log2MeanCount"
colnames(all_counts_mb6pen)[1] = "Log2MeanCount"

#   classify which type of DAR the sites belong to
all_counts_ctrl[sig_up_idx, DAR := "More Accessible in MB6"]
all_counts_ctrl[sig_dn_idx, DAR := "More Accessible in Ctrl"]
all_counts_mb6[sig_up_idx, DAR := "More Accessible in MB6"]
all_counts_mb6[sig_dn_idx, DAR := "More Accessible in Ctrl"]
all_counts_mb6pen[sig_up_idx, DAR := "More Accessible in MB6"]
all_counts_mb6pen[sig_dn_idx, DAR := "More Accessible in Ctrl"]

#   create long form data.table for plotting
all_counts = rbind(
    all_counts_ctrl,
    all_counts_mb6,
    all_counts_mb6pen
)


# ==============================================================================
# Plots
# ==============================================================================
cat("Plotting data\n")
# png(
#     paste0(ARGS$outdir, "/correlation_heatmap_peakcallerscore.png"),
#     width = 20,
#     height = 20,
#     units = "cm",
#     res = 300
# )
# plot(db)
# dev.off()

# png(
#     paste0(ARGS$outdir, "/correlation_heatmap_readcount.png"),
#     width = 20,
#     height = 20,
#     units = "cm",
#     res = 300
# )
# plot(db_counts)
# dev.off()

# png(
#     paste0(ARGS$outdir, "/correlation_heatmap_diffbind.png"),
#     width = 20,
#     height = 20,
#     units = "cm",
#     res = 300
# )
# plot(db_analysis, contrast = 1)
# dev.off()

# png(
#     paste0(ARGS$outdir, "/pca.png"),
#     width = 20,
#     height = 20,
#     units = "cm",
#     res = 300
# )
# dba.plotPCA(db_analysis, contrast = 1, label = DBA_CONDITION)
# dev.off()

# png(
#     paste0(ARGS$outdir, "/ma.png"),
#     width = 20,
#     height = 20,
#     units = "cm",
#     res = 300
# )
# dba.plotMA(db_analysis, contrast = 1)
# dev.off()

# Volcano plot showing DARs and fold change
cat("  Volcano\n")
mb6_ctrl[, Colouring := ifelse(abs(Fold) > FOLD_THRESH & FDR < FDR_THRESH, "firebrick", "gray")]
volcano = (
    ggplot()
    + geom_point(
        data = mb6_ctrl,
        aes(x = Fold, y = -log10(FDR), colour = Colouring)
    )
    + geom_vline(aes(xintercept = c(-1.5, 1.5)), lty = "dashed")
    + geom_hline(aes(yintercept = -log10(FDR_THRESH)), lty = "dashed")
    + scale_colour_identity()
    + labs(
        x = "log2(MB6) - log2(Ctrl)",
        y = "-log10(FDR)",
        title = "MB6 vs Ctrl Accessibility"
    )
    + theme_minimal()
)
ggsave(
    paste0(ARGS$outdir, "/mb6-ctrl.volcano.pdf"),
    height = 20,
    width = 12,
    units = "cm"
)

# top_density = ggplotGrob(
#     ggplot(data = sites)
#     + geom_density(aes(x = Fold))
#     + theme_void()
# )
# side_density = ggplotGrob(
#     ggplot(data = sites)
#     + geom_density(aes(x = -log10(FDR)))
#     + coord_flip()
#     + theme_void()
# )
# png(
#     paste0(ARGS$outdir, "/volcano.png"),
#     width = 12,
#     height = 20,
#     units = "cm",
#     res = 300
# )
# gg = rbind(
#     c(top_density, NA),
#     c(volcano, side_density)
# )
# gg$widths = unit.pmax(top_density$widths, volcano$widths)
# grid.newpage()
# grid.draw(gg)
# dev.off()

# png(
#     paste0(ARGS$outdir, "/boxplot.png"),
#     width = 20,
#     height = 20,
#     units = "cm",
#     res = 300
# )
# dba.plotBox(db_analysis, contrast = 1, label = DBA_CONDITION)
# dev.off()
cat("  Boxplots\n")
# read count boxplot in all peaks
gg_boxplot_all = (
    ggplot(data = all_counts)
    + geom_violin(
        aes(x = Condition, y = Log2MeanCount, fill = Condition)
    )
    + geom_boxplot(
        aes(x = Condition, y = Log2MeanCount),
        width = 0.1
    )
    + labs(
        x = NULL,
        y = "log2(mean read count)",
        title = "Accessibility in All Peaks"
    )
    + guides(fill = FALSE)
    + theme_minimal()
)
ggsave(
    paste0(ARGS$outdir, "/boxplot-all.pdf"),
    height = 12,
    width = 12,
    units = "cm"
)

# read count boxplot in DARs
gg_boxplot_dar = (
    ggplot(data = all_counts[!is.na(DAR)])
    + geom_violin(
        aes(x = Condition, y = Log2MeanCount, fill = Condition)
    )
    + geom_boxplot(
        aes(x = Condition, y = Log2MeanCount),
        width = 0.1
    )
    + labs(
        x = NULL,
        y = "log2(mean read count)",
        title = "Accessibility in Differentially Accessible Peaks"
    )
    + guides(fill = FALSE)
    + facet_wrap(~ DAR)
    + theme_minimal()
)
ggsave(
    paste0(ARGS$outdir, "/boxplot-dar.pdf"),
    height = 12,
    width = 20,
    units = "cm"
)

cat("  Heatmaps")
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
    ggplot(data = mb6_ctrl)
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
