# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("DESeq2"))

# ==============================================================================
# Data
# ==============================================================================
cases = c(
    "Ctrl_Rep1", "Ctrl_Rep2", "Ctrl_Rep3",
    "MB6_Rep1", "MB6_Rep2", "MB6_Rep3"
)

# read in loci
consensus = fread(
    "../2018-11-28_global-accessbility/Consensus/consensus.Ctrl-MB6.bed",
    header = FALSE,
    sep = "\t",
    col.names = c("Chr", "Start", "End")
)
for (case in cases) {
    # read in counts for each case
    dt = fread(
        paste0("Counts/", case, ".sorted.bed"),
        header = FALSE,
        sep = "\t",
        col.names = c("Chr", "Start", "End", case)
    )
    # append columns of counts for each case to the consensus
    consensus = merge(
        consensus,
        dt,
        by = c("Chr", "Start", "End")
    )
}

# consensus_melt = melt(consensus, id.vars=c("Chr", "Start", "End"), variable.name = "Case", value.name = "Count")

# gg <- (
#     ggplot(data = consensus_melt[Count < 1000])
#     + geom_violin(aes(x = Case, y = Count))
#     + labs(x = "Case", y = "Count")
# )
# ggsave(
#     "filename.png",
#     height = 12,
#     width = 20,
#     units = "cm"
# )

# just the count matrix
count_mtx = as.matrix(consensus[, -(1:3)])
loci = consensus[, paste0(Chr, ":", Start, "-", End)]
rownames(count_mtx) = loci

# column data information with rownames in the same order as the counts' colnames
annotation = data.frame(
    Condition = factor(
        rep(c("Ctrl", "MB6"), c(3, 3)),
        levels = c("Ctrl", "MB6")
    ),
    Replicate = factor(c(1, 2, 3, 1, 2, 3)),
    row.names = cases
)

# create DESeqDataSet object
dds = DESeqDataSetFromMatrix(
    countData = count_mtx,
    colData = annotation,
    design = ~ Replicate + Condition
)

# ==============================================================================
# Analysis
# ==============================================================================
# pre-filter out peaks that contain < 60 reads across all samples (avg of 10/sample)
low_counts = which(rowSums(count_mtx) < 60)
loci = loci[-low_counts]
consensus_loci = consensus[-low_counts, .SD, .SDcols = c("Chr", "Start", "End")]
dds = dds[-low_counts, ]

# perform differential calculations
dds = DESeq(dds)

# save results with genomic coordinates
res = results(dds, contrast = c("Condition", "Ctrl", "MB6"))
res_dt = as.data.table(cbind(
    consensus_loci,
    as.data.frame(res)
))

fwrite(
    res_dt,
    "DEseq/Ctrl-vs-MB6.results.tsv",
    sep = "\t",
    col.names = TRUE
)

res_up = res_dt[padj < 0.05 & log2FoldChange > log2(1.5), ]
res_dn = res_dt[padj < 0.05 & log2FoldChange < -log2(1.5), ]

# ==============================================================================
# Plots
# ==============================================================================
# estimates dispersion matrix
png("DEseq/dispersion.png", width = 12, height = 12, units = "cm", res = 300)
plotDispEsts(dds)
dev.off()

# histogram of p-values
gg <- (
    ggplot(data = res_dt)
    + geom_histogram(aes(x = pvalue))
    + labs(
        title = "Ctrl vs MB6 Differential Accessibility",
        subtitle = "Histogram of p-values",
        x = "p-value",
        y = "Frequency"
    )
)
ggsave(
    "DEseq/Ctrl-vs-MB6.pvalues.png",
    height = 12,
    width = 20,
    units = "cm"
)

png("DEseq/Ctrl-vs-MB6.MA.png", width = 12, height = 12, units = "cm", res = 300)
plotMA(res, ylim = c(-5, 5))
dev.off()

# volcano plot
gg <- (
    ggplot(data = res_dt)
    + geom_point(aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.1)
    + geom_hline(aes(yintercept = -log10(0.05)))
    + geom_vline(aes(xintercept = -log2(1.5)))
    + geom_vline(aes(xintercept = log2(1.5)))
    + labs(
        x = "log2(FoldChange)",
        y = "-log10(adjusted p-value)",
        title = "Ctrl vs MB6 Differential Accessibility",
        subtitle = "Volcano plot of peaks"
    )
)
ggsave(
    "Ctrl-vs-MB6.volcano.png",
    height = 12,
    width = 20,
    units = "cm"
)


# RLD Plot
png("RLD.png", width = 12, height = 12, units = "cm", res = 300)
rld = rlog(dds, blind = FALSE)
vsn::meanSdPlot(assay(rld))
dev.off()