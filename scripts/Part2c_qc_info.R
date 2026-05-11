#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Part 2c: Combine NOOB-normalized batches and run ENmix QCinfo
# -------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: Rscript Part2c_qc_info.R <batch_list> <batch_dirs...> <output_dir>")
}

# Parse arguments
batch_list <- as.integer(strsplit(args[1], ",")[[1]])
batch_dirs <- args[2:(length(args) - 1)]
output_dir <- args[length(args)]

cat("=== Part 2c: ENmix QCinfo (Combined MethylSet) ===\n")
cat("Batches:", paste(batch_list, collapse = ", "), "\n")
cat("Output directory:", output_dir, "\n\n")

suppressPackageStartupMessages({
    library(minfi)
    library(ENmix)
})

# Output checks
qcinfo_rds <- file.path(output_dir, "QCinfo_combined.rds")
flag_done  <- file.path(output_dir, ".completed")

if (file.exists(qcinfo_rds) && file.exists(flag_done)) {
    cat("Part 2c already completed — skipping.\n")
    quit(status = 0)
}

tryCatch({

    # -------------------------------------------------------------
    # Load per-batch NOOB-normalized MethylSets
    # -------------------------------------------------------------
    mSet_list <- vector("list", length(batch_list))

    for (i in seq_along(batch_list)) {

        batch     <- batch_list[i]
        batch_dir <- batch_dirs[i]
        mset_file <- file.path(batch_dir,
                               paste0("mSet_noob_batch", batch, ".rds"))

        if (!file.exists(mset_file)) {
            stop("Missing NOOB-normalized MethylSet for batch ", batch)
        }

        cat("Loading batch", batch, "...\n")
        mSet_list[[i]] <- readRDS(mset_file)
        cat("  Samples:", ncol(mSet_list[[i]]), "\n")
    }

    # -------------------------------------------------------------
    # Combine all batches
    # -------------------------------------------------------------
    cat("\nCombining NOOB-normalized batches...\n")
    mSet <- do.call(cbind, mSet_list)
    cat("Combined", ncol(mSet), "samples\n")

    rm(mSet_list)
    gc()

    # -------------------------------------------------------------
    # Run ENmix QCinfo
    # -------------------------------------------------------------
    cat("Running ENmix QCinfo...\n")
    QCinfo <- ENmix::QCinfo(mSet)

    saveRDS(QCinfo, qcinfo_rds)
    cat("Saved QCinfo to:", qcinfo_rds, "\n")

    # Success flag
    writeLines("SUCCESS", flag_done)

    cat("QCinfo completed — intentionally stopping here.\n")
    quit(status = 0)

}, error = function(e) {
    cat("ERROR in Part 2c:", conditionMessage(e), "\n")
    writeLines(
        paste("FAILED:", conditionMessage(e)),
        file.path(output_dir, ".failed")
    )
    quit(status = 1)
})
