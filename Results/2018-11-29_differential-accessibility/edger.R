# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("edgeR"))

# ==============================================================================
# Data
# ==============================================================================
cases = data.table(
    Sample = c(
        "Ctrl_Rep1", "Ctrl_Rep2", "Ctrl_Rep3",
        "MB6_Rep1", "MB6_Rep2", "MB6_Rep3"
    ),
    Condition = factor(rep(c("Ctrl", "MB6"), each = 3)),
    Replicate = factor(rep(1:3, 2))
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

# ==============================================================================
# Analysis
# ==============================================================================
# set up edgeR structure
y = DGEList(counts = consensus[, .SD, .SDcols = 4:9], genes=consensus[, .SD, .SDcols = 1:3])
design = model.matrix(~ cases$Replicate + cases$Condition)

# estimate dispersions
y = estimateDisp(y, design, robust=TRUE)
# general linear model fitting
fit = glmFit(y, design)
# likelihood ratio test
lrt = glmLRT(fit)

lrt$table$QValue = p.adjust(lrt$table$PValue, method = "fdr")

# summary of differential test results
ord = order(lrt$table$QValue)
# get indexes of significantly differentially accessible loci
sig = as.integer(rownames(lrt$table[lrt$table$QValue <= 0.05, ]))
sig_loci = consensus[sig, .(Chr, Start, End)]
decideTestsDGE(lrt)

# ==============================================================================
# Plots
# ==============================================================================
png("EdgeR/md.png", width = 12, height = 12, units = "cm", res = 300)
plotMD(lrt)
abline(h=c(-1, 1), col="blue")
dev.off()

# biological coefficient of variation
png("EdgeR/bcv.png", width = 12, height = 12, units = "cm", res = 300)
plotBCV(y)
dev.off()

# histogram of p-values
gg <- (
    ggplot(data = lrt$table)
    + geom_histogram(aes(x = PValue), binwidth = 0.01)
    + labs(
        title = "Ctrl vs MB6 Differential Accessibility",
        subtitle = "Histogram of p-values",
        x = "p-value",
        y = "Frequency"
    )
)
ggsave(
    "EdgeR/p-values.png",
    height = 12,
    width = 20,
    units = "cm"
)

# histogram of q-values
gg <- (
    ggplot(data = lrt$table)
    + geom_histogram(aes(x = QValue), binwidth = 0.01)
    + labs(
        title = "Ctrl vs MB6 Differential Accessibility",
        subtitle = "Histogram of q-values",
        x = "q-value",
        y = "Frequency"
    )
)
ggsave(
    "EdgeR/q-values.png",
    height = 12,
    width = 20,
    units = "cm"
)
