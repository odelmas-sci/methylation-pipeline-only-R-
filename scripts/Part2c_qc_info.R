#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Part 2c: Combine per-batch RGChannelSets and MethylSets,
#           run ENmix QCinfo on RGChannelSet
# -------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: Rscript Part2c_qc_info.R <batch_list> <batch_dirs...> <output_dir>")
}

batch_list <- as.integer(strsplit(args[1], ",")[[1]])
batch_dirs <- args[2:(length(args) - 1)]
output_dir <- args[length(args)]

cat("=== Part 2c: ENmix QCinfo + Combined MethylSet ===\n")
cat("Batches:", paste(batch_list, collapse = ", "), "\n")
cat("Output directory:", output_dir, "\n\n")

suppressPackageStartupMessages({
    library(minfi)
    library(ENmix)
})

qcinfo_rds <- file.path(output_dir, "QCinfo_combined.rds")
mset_rds   <- file.path(output_dir, "mSet_noob_combined.rds")
flag_done  <- file.path(output_dir, ".completed")

if (file.exists(qcinfo_rds) && file.exists(mset_rds) && file.exists(flag_done)) {
    cat("Part 2c already completed — skipping.\n")
    quit(status = 0)
}

tryCatch({

    # -------------------------------------------------------------
    # Load and combine RGChannelSets (RAW DATA for QCinfo)
    # -------------------------------------------------------------
    rg_list <- vector("list", length(batch_list))

    for (i in seq_along(batch_list)) {

        batch     <- batch_list[i]
        batch_dir <- batch_dirs[i]
        rg_file   <- file.path(batch_dir,
                               paste0("rgSet_batch", batch, ".rds"))

        if (!file.exists(rg_file)) {
            stop("Missing RGChannelSet for batch ", batch)
        }

        cat("Loading RGSet for batch", batch, "...\n")
        rg_list[[i]] <- readRDS(rg_file)
        cat("  Samples:", ncol(rg_list[[i]]), "\n")
    }

    cat("\nCombining RGChannelSets...\n")
    RGSet_combined <- do.call(cbind, rg_list)
    cat("Combined", ncol(RGSet_combined), "samples (RGSet)\n")

    rm(rg_list)
    gc()

    # -------------------------------------------------------------
    # Run ENmix QCinfo on RGChannelSet
    # -------------------------------------------------------------
    cat("Running ENmix QCinfo on combined RGChannelSet...\n")
    QCinfo <- ENmix::QCinfo(RGSet_combined)

    saveRDS(QCinfo, qcinfo_rds)
    cat("Saved QCinfo to:", qcinfo_rds, "\n")

    # -------------------------------------------------------------
    # Load and combine NOOB-normalized MethylSets
    # -------------------------------------------------------------
    mset_list <- vector("list", length(batch_list))

    for (i in seq_along(batch_list)) {

        batch     <- batch_list[i]
        batch_dir <- batch_dirs[i]
        mset_file <- file.path(batch_dir,
                               paste0("mSet_noob_batch", batch, ".rds"))

        if (!file.exists(mset_file)) {
            stop("Missing NOOB-normalized MethylSet for batch ", batch)
        }

        cat("Loading MethylSet for batch", batch, "...\n")
        mset_list[[i]] <- readRDS(mset_file)
        cat("  Samples:", ncol(mset_list[[i]]), "\n")
    }

    cat("\nCombining NOOB-normalized MethylSets...\n")
    mSet_noob_combined <- do.call(cbind, mset_list)
    cat("Combined", ncol(mSet_noob_combined), "samples (MethylSet)\n")

    saveRDS(mSet_noob_combined, mset_rds)
    cat("Saved combined MethylSet to:", mset_rds, "\n")

    rm(mset_list)
    gc()

    writeLines("SUCCESS", flag_done)
    cat("Part 2c completed successfully.\n")
    quit(status = 0)

}, error = function(e) {
    cat("ERROR in Part 2c:", conditionMessage(e), "\n")
    writeLines(
        paste("FAILED:", conditionMessage(e)),
        file.path(output_dir, ".failed")
    )
    quit(status = 1)
})
