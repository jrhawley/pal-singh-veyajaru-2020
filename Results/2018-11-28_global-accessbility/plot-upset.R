# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("UpSetR"))

# ==============================================================================
# Data
# ==============================================================================
consensus = fread(
    "Consensus/consensus-mtx.all.tsv",
    header = TRUE,
    sep = "\t"
)

# ==============================================================================
# Plots
# ==============================================================================
png("Consensus/upset.png", width = 20, height = 12, units = "cm", res = 300)
upset(
    consensus,
    nsets = 8,
    queries = list(
        list(
            query = intersects,
            params = list("Ctrl"),
            active = T
        ),
        list(
            query = intersects,
            params = list("MB6"),
            active = T
        ),
        list(
            query = intersects,
            params = list("MB6Pen"),
            active = T
        ),
        list(
            query = intersects,
            params = list("Ctrl", "MB6Pen"),
            active = T
        ),
        list(
            query = intersects,
            params = list("Ctrl", "MB6", "MB6Pen"),
            active = T
        )
    )
)
dev.off()
