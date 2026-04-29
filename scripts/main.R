#!/usr/bin/env Rscript
# main.R
#
# Orchestrates the 4-part methylation preprocessing pipeline:
#   Stage 1 — Part1_read_idat.R   : run per batch, in parallel
#   Stage 2 — Part2_qc.R          : run per batch, in parallel (after Stage 1)
#   Stage 3 — Part3_betas.R       : single combined run (after Stage 2)
#   Stage 4 — Part4_pca.R         : PCA visualization colored by batch (after Stage 3)
#
# Checkpointing: each stage checks for its output files before running. (Caching to avoid unnecessary reruns.))
# Existing outputs are skipped; rerun specific stages with --force-part{N}.
#
# Usage:
#   Rscript main.R <sample_sheet> <datadir> <outdir> [options]
#
# Options:
#   --batches 1,2,3   Only process these batch numbers (default: all from sample sheet)
#   --cores N         Max parallel workers for Stage 1 and 2 (default: nCores - 1) 
#   --force-part1     Force rerun Stage 1 even if cached
#   --force-part2     Force rerun Stage 2 even if cached
#   --force-part3     Force rerun Stage 3 even if cached
#   --force-part4     Force rerun Stage 4 even if cached
#   --meta-cols A,B   Comma-separated pData column names to color PCA by (default: Batch)
#
# Output layout:
#   <outdir>/batch<N>/     Per-batch outputs from Part 1 and Part 2
#   <outdir>/combined/     Combined outputs from Part 3
#   <outdir>/pca/          PCA plot and coordinates from Part 4
#   <outdir>/pipeline.log  Orchestrator log
#
# Required packages: minfi, parallel, matrixStats, ggplot2

suppressPackageStartupMessages(library(parallel))

# -- Helpers ------------------------------------------------------------------

get_script_dir <- function() {
    f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else normalizePath(".")
}

log_msg <- function(..., file = "") {
    line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste(...))
    cat(line, "\n", sep = "")
    if (nzchar(file)) cat(line, "\n", sep = "", file = file, append = TRUE)
    invisible(NULL)
}

# -- Argument parsing ----------------------------------------------------------

argv <- commandArgs(trailingOnly = TRUE)

if (length(argv) < 3) {
    cat(paste0(
        "Usage: Rscript main.R <sample_sheet> <datadir> <outdir> [options]\n\n",
        "Options:\n",
        "  --batches 1,2,3   Run only specific batches (default: all)\n",
        "  --cores N         Parallel workers (default: physical cores - 1)\n",
        "  --force-part1     Force rerun Stage 1 (Part 1) even if outputs exist\n",
        "  --force-part2     Force rerun Stage 2 (Part 2) even if outputs exist\n",
        "  --force-part3     Force rerun Stage 3 (Part 3) even if outputs exist\n",
        "  --force-part4     Force rerun Stage 4 (Part 4) even if outputs exist\n",
        "  --meta-cols A,B   pData columns to color PCA by, comma-separated (default: Batch)\n"
    ))
    quit(status = 1)
}

sample_sheet <- argv[1]
datadir      <- argv[2]
outdir       <- argv[3]

n_cores     <- max(1L, detectCores(logical = FALSE) - 1L)
batches_arg <- NULL
meta_cols   <- "Batch"
force1 <- force2 <- force3 <- force4 <- FALSE

i <- 4L
while (i <= length(argv)) {
    if (argv[i] == "--batches" && i < length(argv)) {
        batches_arg <- as.integer(strsplit(argv[i + 1L], ",")[[1]])
        i <- i + 2L
    } else if (argv[i] == "--cores" && i < length(argv)) {
        n_cores <- as.integer(argv[i + 1L])
        i <- i + 2L
    } else if (argv[i] == "--meta-cols" && i < length(argv)) {
        meta_cols <- argv[i + 1L]
        i <- i + 2L
    } else if (argv[i] == "--force-part1") { force1 <- TRUE; i <- i + 1L
    } else if (argv[i] == "--force-part2") { force2 <- TRUE; i <- i + 1L
    } else if (argv[i] == "--force-part3") { force3 <- TRUE; i <- i + 1L
    } else if (argv[i] == "--force-part4") { force4 <- TRUE; i <- i + 1L
    } else { i <- i + 1L }
}

# -- Validate inputs -----------------------------------------------------------

if (!file.exists(sample_sheet)) stop("Sample sheet not found: ", sample_sheet)
if (!dir.exists(datadir))       stop("Data directory not found: ", datadir)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pipeline_log <- file.path(outdir, "pipeline.log")
# Reset log file for this run
writeLines(character(0), pipeline_log)

script_dir   <- get_script_dir()
part1_script <- normalizePath(file.path(script_dir, "Part1_read_idat.R"), mustWork = TRUE)
part2_script <- normalizePath(file.path(script_dir, "Part2_qc.R"),        mustWork = TRUE)
part3_script <- normalizePath(file.path(script_dir, "Part3_betas.R"),     mustWork = TRUE)
part4_script <- normalizePath(file.path(script_dir, "Part4_pca.R"),       mustWork = TRUE)

# -- Determine batches from sample sheet ---------------------------------------

sheet   <- read.csv(sample_sheet, stringsAsFactors = FALSE)
targets <- subset(sheet, substr(Sample_Plate, 1L, 1L) == "M")
targets$Batch <- as.integer(factor(targets$Sample_Plate))
all_batches   <- sort(unique(targets$Batch))

batches <- if (!is.null(batches_arg)) batches_arg else all_batches

if (length(invalid <- setdiff(batches, all_batches)) > 0)
    stop("Batch numbers not in sample sheet: ", paste(invalid, collapse = ", "))

combined_dir <- file.path(outdir, "combined")
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("=== Methylation Pipeline Starting ===", file = pipeline_log)
log_msg("Sample sheet : ", sample_sheet,          file = pipeline_log)
log_msg("Data dir     : ", datadir,               file = pipeline_log)
log_msg("Output dir   : ", outdir,                file = pipeline_log)
log_msg("Batches      : ", paste(batches, collapse = ", "), file = pipeline_log)
log_msg("Max workers  : ", n_cores,               file = pipeline_log)
log_msg("Meta columns : ", meta_cols,                file = pipeline_log)
log_msg("Force flags  : part1=", force1, " part2=", force2, " part3=", force3, " part4=", force4, file = pipeline_log)

pipeline_start <- proc.time()[["elapsed"]]

# --------------------------------------------------------------------------------
# STAGE 1 — Part 1: Read IDAT files (parallel per batch)
#
# Cached when:  <outdir>/batch<N>/rgSet_batch<N>.rds  AND  .completed  both exist
# Log file:     <outdir>/batch<N>/part1.log
# --------------------------------------------------------------------------------

log_msg("", file = pipeline_log)
log_msg("--- Stage 1: Reading IDAT files ---", file = pipeline_log)

# Worker runs inside cluster worker (no shared env) — receives all state via list p
worker_part1 <- function(p) {
    batch_dir <- file.path(p$outdir, paste0("batch", p$batch))
    dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)

    rds_out  <- file.path(batch_dir, paste0("rgSet_batch",  p$batch, ".rds"))
    flag_out <- file.path(batch_dir, ".completed")

    if (!p$force && file.exists(rds_out) && file.exists(flag_out)) {
        return(list(batch = p$batch, status = "cached", dir = batch_dir, elapsed = 0L))
    }

    log_file <- file.path(batch_dir, "part1.log")
    t0 <- proc.time()[["elapsed"]]

    rc <- system2("Rscript",
        args = c(p$script, p$batch, p$sample_sheet, p$datadir, batch_dir),
        stdout = log_file, stderr = log_file, wait = TRUE
    )

    elapsed <- as.integer(round(proc.time()[["elapsed"]] - t0))
    ok      <- rc == 0L && file.exists(rds_out) && file.exists(flag_out)
    list(batch = p$batch, status = if (ok) "success" else "failed",
         dir = batch_dir, elapsed = elapsed, log = log_file)
}

params1 <- lapply(batches, function(b) list(
    batch        = b,
    outdir       = outdir,
    sample_sheet = sample_sheet,
    datadir      = datadir,
    script       = part1_script,
    force        = force1
))

n_w1 <- min(length(batches), n_cores)
if (n_w1 > 1L) {
    log_msg(sprintf("  Launching %d batches across %d workers", length(batches), n_w1),
            file = pipeline_log)
    cl <- makeCluster(n_w1, type = "PSOCK")
    part1_res <- tryCatch(
        parLapply(cl, params1, worker_part1),
        finally = stopCluster(cl)
    )
} else {
    part1_res <- lapply(params1, worker_part1)
}

for (r in part1_res) {
    msg <- if (r$status == "cached") {
        sprintf("  Batch %d: cached (skipped)", r$batch)
    } else if (r$status == "success") {
        sprintf("  Batch %d: SUCCESS (%ds)  log: %s", r$batch, r$elapsed, r$log)
    } else {
        sprintf("  Batch %d: FAILED  (%ds)  log: %s", r$batch, r$elapsed, r$log)
    }
    log_msg(msg, file = pipeline_log)
}

part1_ok  <- sapply(part1_res, function(r) r$status %in% c("success", "cached"))
batches_p2 <- batches[part1_ok]
log_msg(sprintf("Stage 1 complete: %d/%d batches OK", sum(part1_ok), length(batches)),
        file = pipeline_log)

if (sum(!part1_ok) > 0) {
    failed_b <- sapply(part1_res[!part1_ok], `[[`, "batch")
    log_msg("  Skipping failed batches in subsequent stages: ",
            paste(failed_b, collapse = ", "), file = pipeline_log)
}

# --------------------------------------------------------------------------------
# STAGE 2 — Part 2: QC analysis (parallel per batch, after Stage 1)
#
# Cached when:  qc_summary_batch<N>.csv  AND  detP_batch<N>.rds  both exist
# Log file:     <outdir>/batch<N>/part2.log
# --------------------------------------------------------------------------------

log_msg("", file = pipeline_log)
log_msg("--- Stage 2: QC Analysis ---", file = pipeline_log)

worker_part2 <- function(p) {
    batch_dir <- file.path(p$outdir, paste0("batch", p$batch))

    qc_csv   <- file.path(batch_dir, paste0("qc_summary_batch", p$batch, ".csv"))
    detp_rds <- file.path(batch_dir, paste0("detP_batch",       p$batch, ".rds"))
    flag_p2  <- file.path(batch_dir, ".completed_part2")

    if (!p$force && file.exists(qc_csv) && file.exists(detp_rds) && file.exists(flag_p2)) {
        return(list(batch = p$batch, status = "cached", dir = batch_dir, elapsed = 0L))
    }

    log_file <- file.path(batch_dir, "part2.log")
    t0 <- proc.time()[["elapsed"]]

    rc <- system2("Rscript",
        args   = c(p$script, p$batch, batch_dir),
        stdout = log_file, stderr = log_file, wait = TRUE
    )

    elapsed <- as.integer(round(proc.time()[["elapsed"]] - t0))
    ok      <- rc == 0L && file.exists(qc_csv) && file.exists(flag_p2)
    list(batch = p$batch, status = if (ok) "success" else "failed",
         dir = batch_dir, elapsed = elapsed, log = log_file)
}

params2 <- lapply(batches_p2, function(b) list(
    batch  = b,
    outdir = outdir,
    script = part2_script,
    force  = force2
))

n_w2 <- min(length(batches_p2), n_cores)
if (n_w2 > 1L) {
    log_msg(sprintf("  Launching %d batches across %d workers", length(batches_p2), n_w2),
            file = pipeline_log)
    cl <- makeCluster(n_w2, type = "PSOCK")
    part2_res <- tryCatch(
        parLapply(cl, params2, worker_part2),
        finally = stopCluster(cl)
    )
} else {
    part2_res <- lapply(params2, worker_part2)
}

for (r in part2_res) {
    msg <- if (r$status == "cached") {
        sprintf("  Batch %d: cached (skipped)", r$batch)
    } else if (r$status == "success") {
        sprintf("  Batch %d: SUCCESS (%ds)  log: %s", r$batch, r$elapsed, r$log)
    } else {
        sprintf("  Batch %d: FAILED  (%ds)  log: %s", r$batch, r$elapsed, r$log)
    }
    log_msg(msg, file = pipeline_log)
}

part2_ok   <- sapply(part2_res, function(r) r$status %in% c("success", "cached"))
batches_p3 <- batches_p2[part2_ok]
log_msg(sprintf("Stage 2 complete: %d/%d batches OK", sum(part2_ok), length(batches_p2)),
        file = pipeline_log)

# --------------------------------------------------------------------------------
# STAGE 3 — Part 3: NOOB normalization (single combined run)
#
# Cached when:  mSet_noob_combined.rds  AND  beta_values_combined.csv  both exist
# Log file:     <outdir>/combined/part3.log
# --------------------------------------------------------------------------------

log_msg("", file = pipeline_log)
log_msg("--- Stage 3: NOOB Normalization ---", file = pipeline_log)

if (length(batches_p3) == 0L) {
    log_msg("ERROR: No batches passed Stage 1 + 2. Cannot run Stage 3.", file = pipeline_log)
    quit(status = 1L)
}

beta_csv <- file.path(combined_dir, "beta_values_combined.csv")
mset_rds <- file.path(combined_dir, "mSet_noob_combined.rds")
flag_p3  <- file.path(combined_dir, ".completed")

if (!force3 && file.exists(beta_csv) && file.exists(mset_rds) && file.exists(flag_p3)) {
    log_msg("  Outputs already exist — skipping (use --force-part3 to rerun)",
            file = pipeline_log)
} else {
    batch_dirs_p3 <- file.path(outdir, paste0("batch", batches_p3))
    log_file3     <- file.path(combined_dir, "part3.log")

    log_msg(sprintf("  Normalizing %d batches: %s",
                    length(batches_p3), paste(batches_p3, collapse = ", ")),
            file = pipeline_log)

    t0 <- proc.time()[["elapsed"]]

    rc <- system2("Rscript",
        args   = c(part3_script,
                   paste(batches_p3, collapse = ","),
                   batch_dirs_p3,
                   combined_dir),
        stdout = log_file3, stderr = log_file3, wait = TRUE
    )

    elapsed <- as.integer(round(proc.time()[["elapsed"]] - t0))

    if (rc != 0L || !file.exists(beta_csv)) {
        log_msg(sprintf("  FAILED (%ds)  log: %s", elapsed, log_file3), file = pipeline_log)
        quit(status = 1L)
    }
    log_msg(sprintf("  SUCCESS (%ds)  log: %s", elapsed, log_file3), file = pipeline_log)
}

# --------------------------------------------------------------------------------
# STAGE 4 — Part 4: PCA visualization (single run, after Stage 3)
#
# Cached when:  pca/pca_by_batch.pdf  AND  pca/.completed_part4  both exist
# Log file:     <outdir>/pca/part4.log
# --------------------------------------------------------------------------------

log_msg("", file = pipeline_log)
log_msg("--- Stage 4: PCA Visualization ---", file = pipeline_log)

pca_dir <- file.path(outdir, "pca")
flag_p4 <- file.path(pca_dir, ".completed_part4")

dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)

# Split comma-separated column string into individual args for Part 4
meta_col_args <- strsplit(meta_cols, ",")[[1]]

log_msg(sprintf("  Metadata columns: %s", paste(meta_col_args, collapse = ", ")),
        file = pipeline_log)

if (!force4 && file.exists(flag_p4)) {
    log_msg("  Outputs already exist — skipping (use --force-part4 to rerun)",
            file = pipeline_log)
} else {
    log_file4 <- file.path(pca_dir, "part4.log")
    t0 <- proc.time()[["elapsed"]]

    rc <- system2("Rscript",
        args   = c(part4_script, combined_dir, pca_dir, meta_col_args),
        stdout = log_file4, stderr = log_file4, wait = TRUE
    )

    elapsed <- as.integer(round(proc.time()[["elapsed"]] - t0))

    if (rc != 0L || !file.exists(flag_p4)) {
        log_msg(sprintf("  FAILED (%ds)  log: %s", elapsed, log_file4), file = pipeline_log)
        quit(status = 1L)
    }
    log_msg(sprintf("  SUCCESS (%ds)  log: %s", elapsed, log_file4), file = pipeline_log)
}

# -- Final summary -------------------------------------------------------------

total_elapsed <- as.integer(round(proc.time()[["elapsed"]] - pipeline_start))

log_msg("", file = pipeline_log)
log_msg("=== Pipeline Complete ===", file = pipeline_log)
log_msg(sprintf("Total wall time : %dm %ds",
                total_elapsed %/% 60L, total_elapsed %% 60L), file = pipeline_log)
log_msg(sprintf("Batches in beta : %s", paste(batches_p3, collapse = ", ")),
        file = pipeline_log)
log_msg(sprintf("Beta values     : %s", beta_csv),  file = pipeline_log)
log_msg(sprintf("mSet (NOOB)     : %s", mset_rds),  file = pipeline_log)
log_msg(sprintf("PCA plots       : %s/pca_by_<col>.pdf", pca_dir), file = pipeline_log)
log_msg(sprintf("Full log        : %s", pipeline_log), file = pipeline_log)
