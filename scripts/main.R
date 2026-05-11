#!/usr/bin/env Rscript
# main.R
#
# Orchestrates the 4-part methylation preprocessing pipeline:
#   Stage 1 — Part1_read_idat.R   : run per batch, in parallel - read IDAT files ("raw data")
#   Stage 2 — Part2_qc.R          : run per batch, in parallel (after Stage 1) - detection P-value QC on raw data
#   Stage 2b — Part2b_noob.R      : run per batch, in parellel (after Stage 1) - Noob Normalization raw data
#   Stage 2c — Part2c_qc_info.R   : single combined run (after Stage 2b) - evaluate QC info from ENmix on Noob Normalized values
#   Stage 3 — Part3_betas.R       : single combined run (after Stage 2c) - calculate Betas and some annotation filtering
#   Stage 4 — Part4_pca.R         : PCA visualization colored by batch (after Stage 3) - generate PCA plots
#
# Checkpointing: each stage checks for its output files before running. (Caching to avoid unnecessary reruns.))
# Existing outputs are skipped; rerun specific stages with --force-part{N}.
#
# Usage:
#   Rscript main.R <sample_sheet> <datadir> <outdir> [options]
#
# Options:
#   --batches 1,2,3   Only process these batch numbers (default: all from sample sheet)
#   --batch-col NAME    Sample sheet column defining plate/batch (default: Sample_Plate)
#   --cores N         Max parallel workers for Stage 1 and 2 (default: nCores - 1) 
#   --force-part1     Force rerun Stage 1 even if cached
#   --force-part2     Force rerun Stage 2 even if cached
#   --force-part2b    Force rerun Stage 2 part b even if cached
#   --force-part2c    Force rerun Stage 2 part c even if cached
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
# Required packages: minfi, parallel, matrixStats, ggplot2, ENmix


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

log_reproducibility_info <- function(file = "") { # enhanced reproducibility
  
  # R version and platform
  rv <- R.version
  log_msg("R Version:  ", rv$version.string, file = file)
  log_msg("Platform:   ", rv$platform, file = file)
  log_msg("OS:         ", utils::sessionInfo()$running, file = file)
  
  # Date/time and timezone
  log_msg("Date/Time:  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), file = file)
  log_msg("Timezone:   ", Sys.timezone(), file = file)
  
  # Locale
  log_msg("Locale:     ", Sys.getlocale("LC_ALL"), file = file)
  
  # Loaded packages and versions
  log_msg("--- Installed Packages ---", file = file)
  pkgs <- as.data.frame(installed.packages()[, c("Package", "Version")],
                        stringsAsFactors = FALSE)
  pkgs <- pkgs[order(pkgs$Package), ]
  if (nrow(pkgs) == 0) {
    log_msg("(no installed packages found)", file = file)
  } else {
    for (i in seq_len(nrow(pkgs))) {
      log_msg(sprintf("  %-25s %s", pkgs$Package[i], pkgs$Version[i]),
              file = file)
    }
  }
  
  # Random seed
  seed <- tryCatch(.Random.seed[2], error = function(e) NA)
  log_msg("RNG Seed:   ", ifelse(is.na(seed), "not yet initialized", seed), file = file)
  
  # Working directory
  log_msg("Working Dir:", getwd(), file = file)
  
  invisible(NULL)
}

# -- Argument parsing ----------------------------------------------------------

argv <- commandArgs(trailingOnly = TRUE)

if (length(argv) < 3) {
    cat(paste0(
        "Usage: Rscript main.R <sample_sheet> <datadir> <outdir> [options]\n\n",
        "Options:\n",
        "  --batches 1,2,3   Run only specific batches (default: all)\n",
        "  --batch-col NAME  Sample sheet column defining plate/batch (default: Sample_Plate)\n",
        "  --cores N         Parallel workers (default: physical cores - 1)\n",
        "  --force-part1     Force rerun Stage 1 (Part 1) even if outputs exist\n",
        "  --force-part2     Force rerun Stage 2 (Part 2) even if outputs exist\n",
        "  --force-part2b     Force rerun Stage 2b (Part 2b) even if outputs exist\n",
        "  --force-part2c     Force rerun Stage 2c (Part 2c) even if outputs exist\n",
        "  --force-part3     Force rerun Stage 3 (Part 3) even if outputs exist\n",
        "  --force-part4     Force rerun Stage 4 (Part 4) even if outputs exist\n",
        "  --meta-cols A,B   pData columns to color PCA by, comma-separated (default: Batch)\n"
    ))
    quit(status = 1)
}

sample_sheet <- argv[1]
datadir      <- argv[2]
outdir       <- argv[3]
batch_col    <- "Sample_Plate"
n_cores     <- max(1L, detectCores(logical = FALSE) - 1L)
batches_arg <- NULL
meta_cols   <- "Batch"
force1 <- force2 <- force2b <- force2c <- force3 <- force4 <- FALSE

i <- 4L
while (i <= length(argv)) {
    if (argv[i] == "--batches" && i < length(argv)) {
        batches_arg <- as.integer(strsplit(argv[i + 1L], ",")[[1]])
        i <- i + 2L
     } else if (argv[i] == "--batch-col" && i < length(argv)) {
        batch_col <- argv[i + 1L]
        i <- i + 2L
    } else if (argv[i] == "--cores" && i < length(argv)) {
        n_cores <- as.integer(argv[i + 1L])
        i <- i + 2L
    } else if (argv[i] == "--meta-cols" && i < length(argv)) {
        meta_cols <- argv[i + 1L]
        i <- i + 2L
    } else if (argv[i] == "--force-part1") { force1 <- TRUE; i <- i + 1L
    } else if (argv[i] == "--force-part2") { force2 <- TRUE; i <- i + 1L
    } else if (argv[i] == "--force-part2b") { force2b <- TRUE; i <- i + 1L
    } else if (argv[i] == "--force-part2c") { force2c <- TRUE; i <- i + 1L
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
part2b_script <- normalizePath(file.path(script_dir, "Part2b_noob.R"),    mustWork = TRUE)
part2c_script <- normalizePath(file.path(script_dir, "Part2c_qc_info.R"), mustWork = TRUE)
part3_script <- normalizePath(file.path(script_dir, "Part3_betas.R"),     mustWork = TRUE)
part4_script <- normalizePath(file.path(script_dir, "Part4_pca.R"),       mustWork = TRUE)

# -- Determine batches from sample sheet ---------------------------------------
sheet <- read.csv(sample_sheet, skip = 8, header = TRUE, stringsAsFactors = FALSE)

if (!batch_col %in% colnames(sheet)) {
    stop("Batch column '", batch_col, "' not found in sample sheet")
}

# Filter for target samples (those starting with "M") ---------------------- > This part may be unnecessary
#targets <- sheet[substr(sheet[[batch_col]], 1L, 1L) == "M", ]
targets <- sheet
                   
# Define batch IDs based on plate column
targets$Batch <- as.integer(factor(targets[[batch_col]]))
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
log_msg("Batch column : ", batch_col, file = pipeline_log)
log_msg("Max workers  : ", n_cores,               file = pipeline_log)
log_msg("Meta columns : ", meta_cols,                file = pipeline_log)
log_msg("Force flags  : part1=", force1, " part2=", force2, " part2b=", force2b, " part2c=", force2c, " part3=", force3, " part4=", force4, file = pipeline_log)

pipeline_start <- proc.time()[["elapsed"]]

log_reproducibility_info()

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
  
  tryCatch({
    
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
    
    ok <- (rc == 0L &&
             file.exists(qc_csv) &&
             file.exists(flag_p2))
    
    return(list(
      batch   = p$batch,
      status  = if (ok) "success" else "failed",
      dir     = batch_dir,
      elapsed = elapsed,
      log     = log_file
    ))
    
  }, error = function(e) {
    
    return(list(
      batch  = p$batch,
      status = "error",
      dir    = file.path(p$outdir, paste0("batch", p$batch)),
      error  = as.character(e)
    ))
    
  })
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

part2_ok   <- unlist(sapply(part2_res, function(r) r$status %in% c("success", "cached")))
batches_p3 <- batches_p2[part2_ok]
log_msg(sprintf("Stage 2 complete: %d/%d batches OK", sum(part2_ok), length(batches_p2)),
        file = pipeline_log)

# --------------------------------------------------------------------------------
# STAGE 2b — Part 2b: Per-batch NOOB normalization
#
# Cached when:  mSet_noob_batch<batch>.rds  AND  .completed  both exist
# Log file:     <outdir>/batch<batch>/part2b.log
# --------------------------------------------------------------------------------

log_msg("", file = pipeline_log)
log_msg("--- Stage 2b: NOOB normalization (per batch) ---", file = pipeline_log)

if (length(batches_p2) == 0L) {
    log_msg("ERROR: No batches passed Stage 1 + 2. Cannot run Stage 2b.",
            file = pipeline_log)
    quit(status = 1L)
}

for (batch in batches_p2) {

    batch_dir  <- file.path(outdir, paste0("batch", batch))
    log_file2b <- file.path(batch_dir, "part2b.log")
    mset_rds   <- file.path(batch_dir,
                            paste0("mSet_noob_batch", batch, ".rds"))
    flag_p2b   <- file.path(batch_dir, ".completed")

    if (!force2b && file.exists(mset_rds) && file.exists(flag_p2b)) {
        log_msg(sprintf("  Batch %d: outputs already exist — skipping (use --force-part2b to rerun)",
                        batch),
                file = pipeline_log)
        next
    }

    log_msg(sprintf("  Batch %d: running NOOB normalization", batch),
            file = pipeline_log)

    t0 <- proc.time()[["elapsed"]]

    rc <- system2(
        "Rscript",
        args   = c(part2b_script,
                   batch,
                   batch_dir,
                   batch_dir),
        stdout = log_file2b,
        stderr = log_file2b,
        wait   = TRUE
    )

    elapsed <- as.integer(round(proc.time()[["elapsed"]] - t0))

    if (rc != 0L || !file.exists(mset_rds)) {
        log_msg(sprintf("  Batch %d: FAILED (%ds)  log: %s",
                        batch, elapsed, log_file2b),
                file = pipeline_log)
        quit(status = 1L)
    }

    log_msg(sprintf("  Batch %d: SUCCESS (%ds)  log: %s",
                    batch, elapsed, log_file2b),
            file = pipeline_log)
}

# --------------------------------------------------------------------------------
# STAGE 2c — Part 2c: ENmix QCinfo (single combined run)
#
# Cached when:  QCinfo_combined.rds  AND  .completed  both exist
# Log file:     <outdir>/combined/part2c.log
# --------------------------------------------------------------------------------

log_msg("", file = pipeline_log)
log_msg("--- Stage 2c: ENmix QCinfo (combined) ---", file = pipeline_log)

qcinfo_rds  <- file.path(combined_dir, "QCinfo_combined.rds")
flag_p2c    <- file.path(combined_dir, ".completed")

if (!force2c && file.exists(qcinfo_rds) && file.exists(flag_p2c)) {
    log_msg("  QCinfo already exists — skipping (use --force-part2c to rerun)",
            file = pipeline_log)
} else {

    log_file2c <- file.path(combined_dir, "part2c.log")

    log_msg(sprintf("  Running QCinfo on %d batches: %s",
                    length(batches_p3),
                    paste(batches_p3, collapse = ", ")),
            file = pipeline_log)

    t0 <- proc.time()[["elapsed"]]

    batch_dirs_p2c <- file.path(outdir, paste0("batch", batches_p3))

    rc <- system2(
        "Rscript",
        args = c(part2c_script,
                 paste(batches_p3, collapse = ","),
                 batch_dirs_p2c,
                 combined_dir),
        stdout = log_file2c,
        stderr = log_file2c,
        wait   = TRUE
    )

    elapsed <- as.integer(round(proc.time()[["elapsed"]] - t0))

    if (rc != 0L || !file.exists(qcinfo_rds)) {
        log_msg(sprintf("  FAILED (%ds)  log: %s", elapsed, log_file2c),
                file = pipeline_log)
        quit(status = 1L)
    }

    log_msg(sprintf("  SUCCESS (%ds)  log: %s", elapsed, log_file2c),
            file = pipeline_log)
}
                            
# --------------------------------------------------------------------------------
# STAGE 3 — Part 3: Calculate Betas (single combined run)
#                        
# Cached when:  beta_values_combined_filtered.csv AND gRatioSet_combined.rds both exist
# Log file:     <outdir>/combined/part3.log
# --------------------------------------------------------------------------------

log_msg("", file = pipeline_log)
log_msg("--- Stage 3: Calculate Betas ---", file = pipeline_log)

if (length(batches_p3) == 0L) {
    log_msg("ERROR: No batches passed Stage 1 + 2. Cannot run Stage 3.", file = pipeline_log)
    quit(status = 1L)
}

# Guard to check for QC information before starting stage 3
qcinfo_rds <- file.path(combined_dir, "QCinfo_combined.rds")
if (!file.exists(qcinfo_rds)) {
    log_msg("ERROR: QCinfo missing. Stage 2c must run before Stage 3.",
            file = pipeline_log)
    quit(status = 1L)
}

beta_csv <- file.path(combined_dir, "beta_values_combined_filtered.csv")
gratio_rds <- file.path(combined_dir, "gRatioSet_combined.rds")
flag_p3 <- file.path(combined_dir, ".completed")

if (!force3 &&
    file.exists(beta_csv) &&
    file.exists(gratio_rds) &&
    file.exists(flag_p3)) {

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
