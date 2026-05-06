#!/usr/bin/env Rscript

# Part 1: Preprocessing for methylation QC pipeline: Get IDAT files

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
    stop("Usage: Rscript Part1_preprocessing.R <batch> <sample_sheet> <datadir> <outdir>")
}

batch <- as.integer(args[1])
sample_sheet <- args[2]
datadir <- args[3]
outdir <- args[4]

cat("=== Part 1: Preprocessing ===\n")
cat("Batch:", batch, "\n")
cat("Sample sheet:", sample_sheet, "\n")
cat("Data directory:", datadir, "\n")
cat("Output directory:", outdir, "\n\n")

# Load required libraries
suppressPackageStartupMessages({
    if (!require("minfi", quietly = TRUE)) {
        stop("Required package 'minfi' not installed. Install with: BiocManager::install('minfi')")
    }
})

# Create output directory if needed
if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
}

# Read sample sheet
sheet   <- read.csv(sample_sheet, skip = 8, header = TRUE,  stringsAsFactors = FALSE)
sheet$Basename <- paste0(sheet$Sentrix_ID, "_", sheet$Sentrix_Position)

# Filter for target samples (those starting with "M")
targets <- subset(sheet, substr(Sample_Plate, 1, 1) == "M")
targets$Batch <- as.integer(factor(targets$Sample_Plate))


# Filter for current batch
batch_targets <- subset(targets, Batch == batch)

cat("Processing", nrow(batch_targets), "samples in batch", batch, "\n")

# Read IDAT files
tryCatch({
    rgSet <- read.metharray.exp(targets = batch_targets, base = datadir)
    
    # Save raw RGChannelSet
    saveRDS(rgSet, file = file.path(outdir, paste0("rgSet_batch", batch, ".rds")))
    cat("Saved rgSet to:", file.path(outdir, paste0("rgSet_batch", batch, ".rds")), "\n")
    
    # Generate QC report
    qcReport(rgSet, 
             sampNames = batch_targets$Sample_Name,
             pdf = file.path(outdir, paste0("qcReport_batch", batch, ".pdf")))
    
    # Save batch targets
    write.csv(batch_targets, 
              file = file.path(outdir, paste0("targets_batch", batch, ".csv")),
              row.names = FALSE)
    
    # Write success flag
    writeLines("SUCCESS", file.path(outdir, ".completed"))
    cat("Batch", batch, "completed successfully\n")
    
}, error = function(e) {
    cat("ERROR in batch", batch, ":", conditionMessage(e), "\n")
    writeLines(paste("FAILED:", conditionMessage(e)), 
               file.path(outdir, ".failed"))
    quit(status = 1)
})
