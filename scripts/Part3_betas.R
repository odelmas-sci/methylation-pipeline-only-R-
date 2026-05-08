#!/usr/bin/env Rscript

# Part 3: NOOB normalization (Memory optimized)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: Rscript betas.R <batch_list> <batch_dirs...> <output_dir>")
}

# Parse arguments
batch_list <- as.integer(strsplit(args[1], ",")[[1]])
output_dir <- args[length(args)]
batch_dirs <- args[2:(length(args)-1)]

cat("=== Part 3: Calculating Betas (Combined) ===\n")
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
mSet_list <- list()

for (i in seq_along(batch_list)) {
    batch     <- batch_list[i]
    batch_dir <- batch_dirs[i]

    mset_file <- file.path(
        batch_dir,
        paste0("mSet_noob_batch", batch, ".rds")
    )

    if (!file.exists(mset_file)) {
        stop("Missing mSet file for batch ", batch)
    }

    cat("Loading NOOB-normalized batch", batch, "...\n")
    mSet_list[[i]] <- readRDS(mset_file)
    cat("  Samples:", ncol(mSet_list[[i]]), "\n")
}

cat("\nCombining NOOB-normalized batches...\n")
mSet <- do.call(cbind, mSet_list)
cat("Combined", ncol(mSet), "samples\n")

rm(mSet_list)
gc()

cat("Running ENmix QCinfo...\n")
QCinfo <- ENmix::QCinfo(mSet)

saveRDS(QCinfo, file.path(output_dir, "QCinfo_combined.rds"))

cat("QCinfo completed — stopping here temporarily.\n")
quit(status = 0) #--------------------------------------------------------> TEMPORARY STOP IN THE PIPELINE, comment/uncomment when needed

cat("Converting to ratio set...\n")
gRatioSet <- ratioConvert(mSet, what = "both", keepCN = TRUE)

cat("Extracting beta values...\n")
beta <- getBeta(gRatioSet)

cat("\nApplying annotation-based probe filtering...\n")

# Remove certain annotations 
cat("\nLoading annotation dataset...\n")
annot <- as.data.frame(read.table("scripts/EPIC.hg38.manifest.tsv",
                                  sep = "\t", header = T))
rownames(annot) = annot$probeID

cat("\nFiltering annotation...\n")
    
# Filter crossreactive, polymorphic probes at the single base extension, X, Y, M and ch
annot_filt <- annot[
  annot$CpG_chrm != "chrX" &
  annot$CpG_chrm != "chrY" &
  annot$CpG_chrm != "chrM" &
  annot$probeType == "cg" &
  annot$MASK_general == FALSE,
]

cat("Saving annotation-filtered probe lists...\n")

# All probes in original annotation
all_probes <- rownames(annot)
# Probes that pass annotation filtering
kept_probes_annot <- rownames(annot_filt)
# Probes removed by annotation rules
removed_probes_annot <- setdiff(all_probes, kept_probes_annot)

write.table(
  removed_probes_annot,
  file = file.path(output_dir, "probes_removed_annotation_filters.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("Removed", length(removed_probes_annot), "probes via annotation filtering\n")

cat("Subsetting Beta to filtered probes...\n")
    
keep_probes <- intersect(rownames(beta), rownames(annot_filt))

beta <- beta[keep_probes, ]

cat("Beta matrix now contains", nrow(beta), "probes and", ncol(beta), "samples\n")

# Save results
cat("Saving results...\n")
saveRDS(mSet, file = file.path(output_dir, "mSet_noob_combined.rds"))
saveRDS(gRatioSet, file = file.path(output_dir, "gRatioSet_combined.rds"))
write.csv(beta, file = file.path(output_dir, "beta_values_combined.csv"), row.names = TRUE)

# create success/failure flag for clear data provenance
writeLines("SUCCESS", file.path(output_dir, ".completed"))
cat("NOOB normalization completed!\n")

}, error = function(e) {
    cat("ERROR in Part 3:", conditionMessage(e), "\n")
    writeLines(paste("FAILED:", conditionMessage(e)), file.path(output_dir, ".failed"))
    quit(status = 1)
})
