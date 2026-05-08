#!/usr/bin/env Rscript

# Part 2b: Per-batch NOOB normalization

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: Rscript Part2b_noob.R <batch> <batch_dir> <outdir>")
}

batch     <- as.integer(args[1])
batch_dir <- args[2]
outdir    <- args[3]

cat("=== Part 2b: NOOB normalization ===\n")
cat("Batch:", batch, "\n")
cat("Batch dir:", batch_dir, "\n")
cat("Output dir:", outdir, "\n\n")

suppressPackageStartupMessages({
    library(minfi)
})

if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
}

flag_done <- file.path(outdir, ".completed")
mset_file <- file.path(outdir,
                       paste0("mSet_noob_batch", batch, ".rds"))

if (file.exists(flag_done) && file.exists(mset_file)) {
    cat("Part2b already completed — skipping.\n")
    quit(status = 0)
}

tryCatch({

    rg_file <- file.path(batch_dir,
                         paste0("rgSet_batch", batch, ".rds"))

    if (!file.exists(rg_file)) {
        stop("RGSet not found: ", rg_file)
    }

    cat("Loading RGSet...\n")
    rgSet <- readRDS(rg_file)
    cat("Samples:", ncol(rgSet), "\n")

    cat("Performing NOOB normalization...\n")
    mSet <- preprocessNoob(rgSet)

    saveRDS(mSet, mset_file)
    cat("Saved:", mset_file, "\n")

    writeLines("SUCCESS", flag_done)
    cat("Batch", batch, "completed successfully\n")

}, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    writeLines(
        paste("FAILED:", conditionMessage(e)),
        file.path(outdir, ".failed")
    )
    quit(status = 1)
})
