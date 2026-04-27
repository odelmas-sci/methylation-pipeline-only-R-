#!/usr/bin/env Rscript

# Part 2: Quality control analysis

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: Rscript Part2_qc.R <batch> <batch_dir>")
}

batch <- as.integer(args[1])
batch_dir <- args[2]

cat("=== Part 2: Quality Control ===\n")
cat("Batch:", batch, "\n")
cat("Batch directory:", batch_dir, "\n\n")

# Skip if already completed
flag_done  <- file.path(batch_dir, ".completed_part2")
qc_csv_chk <- file.path(batch_dir, paste0("qc_summary_batch", batch, ".csv"))
detp_chk   <- file.path(batch_dir, paste0("detP_batch", batch, ".rds"))
if (file.exists(flag_done) && file.exists(qc_csv_chk) && file.exists(detp_chk)) {
    cat("Part 2 already completed for batch", batch, "— skipping.\n")
    quit(status = 0)
}

# Load required libraries
suppressPackageStartupMessages({
    if (!require("minfi", quietly = TRUE)) {
        stop("Required package 'minfi' not installed")
    }
})

# Load RGChannelSet from Part 1
rgset_file <- file.path(batch_dir, paste0("rgSet_batch", batch, ".rds"))
if (!file.exists(rgset_file)) {
    stop("RGChannelSet file not found: ", rgset_file)
}

rgSet <- readRDS(rgset_file)
cat("Loaded rgSet with", ncol(rgSet), "samples\n")

# Perform QC checks
tryCatch({
    # Calculate detection p-values
    detP <- detectionP(rgSet)
    
    # Calculate mean detection p-values per sample
    mean_detP <- colMeans(detP)
    
    # Identify failed samples (mean detection p-value > 0.05)
    failed <- mean_detP > 0.05
    cat("Number of failed samples:", sum(failed), "\n")
    
    # Create QC summary
    qc_summary <- data.frame(
        Sample = colnames(rgSet),
        Mean_DetP = mean_detP,
        Failed = failed,
        stringsAsFactors = FALSE
    )
    
    # Save QC summary
    write.csv(qc_summary, 
              file = file.path(batch_dir, paste0("qc_summary_batch", batch, ".csv")),
              row.names = FALSE)
    
    # Save detection p-values
    saveRDS(detP, file = file.path(batch_dir, paste0("detP_batch", batch, ".rds")))
    
    # Generate QC plots
    pdf(file.path(batch_dir, paste0("qc_plots_batch", batch, ".pdf")), width = 10, height = 8)
    
    # Plot 1: Mean detection p-value barplot
    barplot(mean_detP, 
            main = paste("Mean Detection P-values - Batch", batch),
            las = 2,
            ylab = "Mean Detection P-value",
            col = ifelse(failed, "red", "green"))
    abline(h = 0.05, col = "blue", lty = 2)
    
    # Plot 2: Sample-wise detection p-value distribution
    boxplot(detP, 
            main = paste("Detection P-value Distribution - Batch", batch),
            las = 2,
            ylab = "Detection P-value",
            outline = FALSE)
    abline(h = 0.05, col = "blue", lty = 2)
    
    dev.off()
    
    # create success/failure flag for clear data prrovenance
    writeLines("SUCCESS", file.path(batch_dir, ".completed_part2"))
    cat("Part 2 QC completed successfully!\n")
    cat("QC summary saved to:", file.path(batch_dir, paste0("qc_summary_batch", batch, ".csv")), "\n")

}, error = function(e) {
    cat("Error in Part 2 QC:", conditionMessage(e), "\n")
    writeLines(paste("FAILED:", conditionMessage(e)), file.path(batch_dir, ".failed_part2"))
    quit(status = 1)
})
