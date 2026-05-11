#!/usr/bin/env Rscript

# Part 3: Calculate Beta values (post-QCinfo, combined data)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: Rscript Part3_betas.R <batch_list> <batch_dirs...> <output_dir>")
}

# Parse arguments
batch_list <- as.integer(strsplit(args[1], ",")[[1]])
batch_dirs <- args[2:(length(args) - 1)]
output_dir <- args[length(args)]

cat("=== Part 3: Calculating Betas (Combined) ===\n")
cat("Batches:", paste(batch_list, collapse = ", "), "\n")
cat("Output directory:", output_dir, "\n\n")

suppressPackageStartupMessages({
    library(minfi)
})

# Define required inputs and outputs
qcinfo_file <- file.path(output_dir, "QCinfo_combined.rds")
mset_file   <- file.path(output_dir, "mSet_noob_combined.rds")

beta_csv_final <- file.path(output_dir, "beta_values_combined_filtered.csv")
gratio_rds     <- file.path(output_dir, "gRatioSet_combined.rds")
flag_done      <- file.path(output_dir, ".completed")

# Skip logic: Part 3 is complete if FINAL outputs exist
if (file.exists(flag_done) && file.exists(beta_csv_final) && file.exists(gratio_rds)) {
    cat("Part 3 already completed — skipping.\n")
    quit(status = 0)
}

# Hard guards: upstream stages must exist
if (!file.exists(qcinfo_file)) {
    stop("QCinfo not found. Run Part2c_qc_info.R first.")
}
if (!file.exists(mset_file)) {
    stop("Combined MethylSet not found. Run Part2b + Part2c first.")
}

tryCatch({

    # Load inputs
    cat("Loading combined MethylSet...\n")
    mSet <- readRDS(mset_file)
    cat("  Samples:", ncol(mSet), "\n")

    cat("Loading QCinfo (diagnostic)...\n")
    QCinfo <- readRDS(qcinfo_file)
    cat("  QCinfo samples:", nrow(QCinfo$sample), "\n")

    
    # Convert to RatioSet and extract beta values
    cat("Converting to RatioSet...\n")
    gRatioSet <- ratioConvert(mSet, what = "both", keepCN = TRUE)

    cat("Extracting beta values...\n")
    beta <- getBeta(gRatioSet)

    # Annotation-based probe filtering
    cat("\nLoading annotation dataset...\n")
    annot <- as.data.frame(
        read.table(
            "scripts/EPIC.hg38.manifest.tsv",
            sep = "\t", header = TRUE, stringsAsFactors = FALSE
        )
    )
    rownames(annot) <- annot$probeID

    cat("Applying annotation filters...\n")
    annot_filt <- annot[
        annot$CpG_chrm != "chrX" &
        annot$CpG_chrm != "chrY" &
        annot$CpG_chrm != "chrM" &
        annot$probeType == "cg" &
        annot$MASK_general == FALSE,
    ]

    # Record removed probes
    removed_probes_annot <- setdiff(rownames(annot), rownames(annot_filt))
    write.table(
        removed_probes_annot,
        file = file.path(output_dir, "probes_removed_annotation_filters.txt"),
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )
    cat("Removed", length(removed_probes_annot), "probes via annotation filtering\n")

    # Subset beta matrix
    keep_probes <- intersect(rownames(beta), rownames(annot_filt))
    beta <- beta[keep_probes, , drop = FALSE]

    cat("Final beta matrix:",
        nrow(beta), "probes x", ncol(beta), "samples\n")

    # Save outputs (Part 3 products only)
    cat("Saving outputs...\n")
    saveRDS(gRatioSet, gratio_rds)
    write.csv(beta, beta_csv_final, row.names = TRUE)

    writeLines("SUCCESS", flag_done)
    cat("Part 3 completed successfully.\n")

}, error = function(e) {
    cat("ERROR in Part 3:", conditionMessage(e), "\n")
    writeLines(
        paste("FAILED:", conditionMessage(e)),
        file.path(output_dir, ".failed")
    )
    quit(status = 1)
})
