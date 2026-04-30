#!/usr/bin/env Rscript

# Part 4: PCA visualization of NOOB-normalized beta values
#
# Usage: Rscript Part4_pca.R <combined_dir> <outdir> [col1 col2 ...]
#
# One png is produced per metadata column. Columns must exist in pData of
# gRatioSet_combined.rds. Numeric columns get a continuous color scale;
# character/factor columns get a discrete palette.
# Defaults to "Batch" if no columns are supplied.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: Rscript Part4_pca.R <combined_dir> <outdir> [col1 col2 ...]")
}

combined_dir <- args[1]
outdir       <- args[2]
meta_cols    <- if (length(args) > 2) args[3:length(args)] else "Batch"

cat("=== Part 4: PCA Visualization ===\n")
cat("Combined directory:", combined_dir, "\n")
cat("Output directory  :", outdir, "\n")
cat("Metadata columns  :", paste(meta_cols, collapse = ", "), "\n\n")

# Verify Part 3 completed cleanly before proceeding
flag_p3    <- file.path(combined_dir, ".completed")
gratio_rds <- file.path(combined_dir, "gRatioSet_combined.rds")

if (!file.exists(flag_p3)) {
    stop("Part 3 completion flag not found — rerun Part 3 before Part 4: ", flag_p3)
}
if (!file.exists(gratio_rds)) {
    stop("gRatioSet not found: ", gratio_rds)
}

# Load required libraries
suppressPackageStartupMessages({
    if (!require("minfi", quietly = TRUE))
        stop("Required package 'minfi' not installed. Install with: BiocManager::install('minfi')")
    if (!require("matrixStats", quietly = TRUE))
        stop("Required package 'matrixStats' not installed. Install with: install.packages('matrixStats')")
    if (!require("ggplot2", quietly = TRUE))
        stop("Required package 'ggplot2' not installed. Install with: install.packages('ggplot2')")
})

# Create output directory
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Skip if already completed (rerun with --force-part4 to change columns)
flag_done <- file.path(outdir, ".completed_part4")
if (file.exists(flag_done)) {
    cat("Part 4 already completed — skipping.\n")
    cat("Use --force-part4 in run_pipeline.R to rerun with different columns.\n")
    quit(status = 0)
}

tryCatch({

    cat("Loading gRatioSet...\n")
    gRatioSet <- readRDS(gratio_rds)

    beta <- getBeta(gRatioSet)
    pd   <- as.data.frame(pData(gRatioSet))

    cat("Beta matrix    :", nrow(beta), "probes x", ncol(beta), "samples\n")

    # Validate requested columns
    missing_cols <- setdiff(meta_cols, colnames(pd))
    if (length(missing_cols) > 0) {
        stop("Column(s) not found in pData: ", paste(missing_cols, collapse = ", "),
             "\nAvailable columns: ", paste(colnames(pd), collapse = ", "))
    }

    # Subset to top 40k most variable probes
    n_probes <- min(40000L, nrow(beta))
    cat(sprintf("Selecting top %d most variable probes...\n", n_probes))
    probe_var  <- matrixStats::rowVars(beta, na.rm = TRUE)
    top_probes <- order(probe_var, decreasing = TRUE)[seq_len(n_probes)]
    beta_sub   <- beta[top_probes, ]

    # prcomp expects samples as rows
    cat("Running PCA...\n")
    cat("Setting seed to 42 for reproducibility...\n")
    set.seed(42)
    pca     <- prcomp(t(beta_sub), center = TRUE, scale. = FALSE)
    pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

    x_lab <- sprintf("PC1 (%.1f%% variance)", pct_var[1])
    y_lab <- sprintf("PC2 (%.1f%% variance)", pct_var[2])

    # Build base data frame with PCs and all requested columns
    pca_df <- data.frame(
        PC1    = pca$x[, 1],
        PC2    = pca$x[, 2],
        Sample = colnames(beta),
        stringsAsFactors = FALSE
    )
    for (col in meta_cols) pca_df[[col]] <- pd[[col]]

    # Save coordinates (includes all metadata columns)
    pca_csv <- file.path(outdir, "pca_coordinates.csv")
    write.csv(pca_df, file = pca_csv, row.names = FALSE)
    cat("Coordinates saved to:", pca_csv, "\n")

    # One plot per metadata column
    for (col in meta_cols) {
        cat(sprintf("Generating PCA plot colored by '%s'...\n", col))

        color_val <- pd[[col]]

        # Discrete palette for character/factor; continuous for numeric
        if (is.numeric(color_val)) {
            color_scale <- scale_color_viridis_c(name = col)
        } else {
            color_val   <- factor(color_val)
            color_scale <- scale_color_brewer(name = col, palette = "Set1",
                                              drop = FALSE)
        }

        p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = color_val)) +
            geom_point(size = 2.5, alpha = 0.8) +
            color_scale +
            labs(
                title = sprintf("PCA colored by %s", col),
                x     = x_lab,
                y     = y_lab
            ) +
            theme_bw() +
            theme(
                plot.title      = element_text(hjust = 0.5, face = 'bold'),
                legend.position = "bottom"
            )

        out_png <- file.path(outdir, sprintf("pca_by_%s.png", col))
        ggsave(out_png, plot = p, width = 7, height = 7, dpi = 600)
        cat("  Saved:", out_png, "\n")
    }

    writeLines(paste(meta_cols, collapse = ","), flag_done)
    cat("Part 4 completed successfully!\n")

}, error = function(e) {
    cat("ERROR in Part 4:", conditionMessage(e), "\n")
    writeLines(paste("FAILED:", conditionMessage(e)), file.path(outdir, ".failed_part4"))
    quit(status = 1)
})
