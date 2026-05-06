#!/usr/bin/env Rscript

# Part 3: NOOB normalization (Memory optimized)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: Rscript noob_normalize.R <batch_list> <batch_dirs...> <output_dir>")
}

# Parse arguments
batch_list <- as.integer(strsplit(args[1], ",")[[1]])
output_dir <- args[length(args)]
batch_dirs <- args[2:(length(args)-1)]

cat("=== Part 3: NOOB Normalization (Combined) ===\n")
cat("Batches:", paste(batch_list, collapse=", "), "\n")
cat("Output directory:", output_dir, "\n\n")

# Skip if already completed
flag_done    <- file.path(output_dir, ".completed")
beta_csv_chk <- file.path(output_dir, "beta_values_combined.csv")
mset_rds_chk <- file.path(output_dir, "mSet_noob_combined.rds")
if (file.exists(flag_done) && file.exists(beta_csv_chk) && file.exists(mset_rds_chk)) {
    cat("Part 3 already completed — skipping.\n")
    quit(status = 0)
}

tryCatch({

# Load all RGChannelSets and combine
rgSet_list <- list()
for (i in seq_along(batch_list)) {
    batch <- batch_list[i]
    batch_dir <- batch_dirs[i]

    rgset_file <- file.path(batch_dir, paste0("rgSet_batch", batch, ".rds"))
    if (file.exists(rgset_file)) {
        cat("Loading batch", batch, "...\n")
        rgSet_list[[i]] <- readRDS(rgset_file)
        cat("Loaded batch", batch, "with", ncol(rgSet_list[[i]]), "samples\n")
    }
}

# Combine all batches
cat("\nCombining batches...\n")
rgSet_combined <- do.call(cbind, rgSet_list)
cat("Combined", ncol(rgSet_combined), "samples across", length(batch_list), "batches\n")

# Free memory from list
rm(rgSet_list)
gc()

# Perform NOOB normalization on combined data
cat("Performing NOOB normalization on combined dataset...\n")
mSet <- preprocessNoob(rgSet_combined)

# Free memory from raw data
rm(rgSet_combined)
gc()

cat("Converting to ratio set...\n")
gRatioSet <- ratioConvert(mSet, what = "both", keepCN = TRUE)

cat("Extracting beta values...\n")
beta <- getBeta(gRatioSet)

cat("Loading annotation...")
annot <- as.data.frame(read.table"scripts/EPIC.hg38.manifest.tsv",
                                  sep = "\t", header = T)
rownames(annot) = annot$probeID

annot<-annot[rownames(annot) %in% rownames(beta),]

#filter crossreactive, polymorphic probes at the single base extension, X, Y, M and ch
betas2 <- beta[rownames(beta) %in% rownames(annot[annot$CpG_chrm!="chrY" &
                                            annot$CpG_chrm!="chrX" &
                                            annot$CpG_chrm!="chrM" &
                                            annot$probeType=="cg" &
                                            annot$MASK_general==FALSE,]),]
print(dim(betas2))


# Save results
cat("Saving results...\n")
saveRDS(mSet, file = file.path(output_dir, "mSet_noob_combined.rds"))
saveRDS(gRatioSet, file = file.path(output_dir, "gRatioSet_combined.rds"))
write.csv(beta, file = file.path(output_dir, "beta_values_combined.csv"), row.names = TRUE)
write.csv(betas2, file = file.path(output_dir,"beta_values_combined_filtered.csv"), row.names = TRUE)

# create success/failure flag for clear data provenance
writeLines("SUCCESS", file.path(output_dir, ".completed"))
cat("NOOB normalization completed!\n")

}, error = function(e) {
    cat("ERROR in Part 3:", conditionMessage(e), "\n")
    writeLines(paste("FAILED:", conditionMessage(e)), file.path(output_dir, ".failed"))
    quit(status = 1)
})
