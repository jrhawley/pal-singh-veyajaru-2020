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
db_counts = dba.count(db, score = DBA_SCORE_TMM_READS_FULL_CPM)
db_peaks = dba.peakset(db_counts, bRetrieve = TRUE)

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
# db_contrast_mb6_ctrl = dba.contrast(
#     db_counts,
#     mask_mb6,
#     mask_ctrl,
#     "MB6",
#     "Ctrl"
# )
# # ensure statistics will be of the form MB6Pen - Ctrl
# db_contrast_mb6pen_ctrl = dba.contrast(
#     db_counts,
#     mask_mb6pen,
#     mask_ctrl,
#     "MB6Pen",
#     "Ctrl",
#     minMembers = 2
# )
# # ensure statistics will be of the form MB6 - MB6Pen
# db_contrast_mb6_mb6pen = dba.contrast(
#     db_counts,
#     mask_mb6,
#     mask_mb6pen,
#     "MB6",
#     "MB6Pen",
#     minMembers = 2
# )

cat("Analyzing data\n")
db_analysis = dba.analyze(db_contrast, bReduceObjects = FALSE)

# Extract DARs
# -------------------------------------
# get DARs from analyses
sites_gr_mb6_ctrl = dba.report(db_analysis, contrast = 1, th = 1)
sites_gr_mb6pen_ctrl = dba.report(db_analysis, contrast = 2, th = 1)
sites_gr_mb6_mb6pen = dba.report(db_analysis, contrast = 3, th = 1)

# set order for GRanges objects
seqlevels(sites_gr_mb6_ctrl) = CHRS
seqlevels(sites_gr_mb6pen_ctrl) = CHRS
seqlevels(sites_gr_mb6_mb6pen) = CHRS
sites_gr_mb6_ctrl = sort(sites_gr_mb6_ctrl)
sites_gr_mb6pen_ctrl = sort(sites_gr_mb6pen_ctrl)
sites_gr_mb6_mb6pen = sort(sites_gr_mb6_mb6pen)

# convert to data.tables
sites_mb6_ctrl = as.data.table(sites_gr_mb6_ctrl)
sites_mb6pen_ctrl = as.data.table(sites_gr_mb6pen_ctrl)
sites_mb6_mb6pen = as.data.table(sites_gr_mb6_mb6pen)

# convert start coordinates back to 0-indexed for writing
sites_mb6_ctrl[, start := start - 1]
sites_mb6pen_ctrl[, start := start - 1]
sites_mb6_mb6pen[, start := start - 1]

# switch comparison direction so:
#   Fold = log2(MB6) - log2(Ctrl)
#   Fold = log2(MB6Pen) - log2(Ctrl)
#   Fold = log2(MB6) - log2(MB6Pen)
sites_mb6_ctrl[, Fold := -Fold]
sites_mb6pen_ctrl[, Fold := -Fold]

# save R objects for future loading, if needed
saveRDS(db_analysis, paste0(ARGS$outdir, "/diffbind-analysis.rds"))
saveRDS(db_peaks, paste0(ARGS$outdir, "/diffbind-peaks-cpm.rds"))

# Get CPM calues for DARs
# -------------------------------------
# convert to data.table
sig_cpm = as.data.table(db_peaks[sig_idx])

# average CPM across conditions
sig_cpm[, Ctrl := mapply(mean, Ctrl1, Ctrl2, Ctrl3)]
sig_cpm[, c("Ctrl1", "Ctrl2", "Ctrl3") := NULL]
sig_cpm[, MB6 := mapply(mean, MB61, MB62, MB63)]
sig_cpm[, c("MB61", "MB62", "MB63") := NULL]
sig_cpm[, MB6Pen := mapply(mean, MB6Pen1, MB6Pen2)]
sig_cpm[, c("MB6Pen1", "MB6Pen2") := NULL]

# convert to long form
sig_cpm_ctrl = sig_cpm[, .SD, .SDcols = c("seqnames", "start", "end", "Ctrl")]
sig_cpm_ctrl[, Condition := "Ctrl"]
colnames(sig_cpm_ctrl) = c("chr", "start", "end", "CPM", "Condition")

sig_cpm_mb6 = sig_cpm[, .SD, .SDcols = c("seqnames", "start", "end", "MB6")]
sig_cpm_mb6[, Condition := "MB6"]
colnames(sig_cpm_mb6) = c("chr", "start", "end", "CPM", "Condition")

sig_cpm_mb6pen = sig_cpm[, .SD, .SDcols = c("seqnames", "start", "end", "MB6Pen")]
sig_cpm_mb6pen[, Condition := "MB6Pen"]
colnames(sig_cpm_mb6pen) = c("chr", "start", "end", "CPM", "Condition")

sig_cpm_long = rbind(
    sig_cpm_ctrl,
    sig_cpm_mb6,
    sig_cpm_mb6pen
)

# Save tested regions
# -------------------------------------
agg_sites = list(sites_mb6_ctrl, sites_mb6pen_ctrl, sites_mb6_mb6pen)
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

# get indices for DARs that pass fold threshold
sig_idx = sites_mb6_ctrl[, which(FDR < FDR_THRESH & abs(Fold) >= FOLD_THRESH)]

# Calculate the number of significant sites at a certain fold-change threshold
# -------------------------------------
deltas = seq(0, 3, 0.05)
sizes = data.table(
    Threshold = deltas,
    Increased = sapply(
        deltas,
        function(d) {sites_mb6_ctrl[Fold >= d & FDR < FDR_THRESH, .N]}
    ),
    Decreased = sapply(
        deltas,
        function(d) {sites_mb6_ctrl[-Fold >= d & FDR < FDR_THRESH, .N]}
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
sites_mb6_ctrl[, Colouring := ifelse(abs(Fold) > FOLD_THRESH & FDR < FDR_THRESH, "firebrick", "gray")]
volcano = (
    ggplot()
    + geom_point(
        data = sites_mb6_ctrl,
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
cat("  Boxplot\n")
gg_boxplot = (
    ggplot(data = sig_cpm_long)
    + geom_boxplot(aes(x = Condition, y = CPM, fill = Condition))
    + labs(
        x = "Condition",
        y = "CPM",
        title = "Chromatin Accessibility"
    )
    + theme_minimal()
)
ggsave(
    paste0(ARGS$outdir, "/boxplot.pdf"),
    height = 12,
    width = 20,
    units = "cm"
)


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
