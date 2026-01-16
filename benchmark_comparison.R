#!/usr/bin/env Rscript

# Parameters
OLD_COMMIT <- "527cd54b24356f765e3829a94629cc00cde0754f"
N_ITER <- 5 # Reduced for robustness
WORKER_SCRIPT <- "_benchmark_worker.R"

suppressPackageStartupMessages(library(data.table))

# Pull in helper suite so that the loops used for comparison are accessible here as well.
helper_file <- file.path(getwd(), "benchmark_inference_loops.R")
if (file.exists(helper_file)) {
  source(helper_file)
}

# Setup Paths
wd <- getwd()
tmp_dir <- tempdir()
repo_dir <- file.path(tmp_dir, "repo_old")
lib_old <- file.path(tmp_dir, "lib_old")
lib_new <- file.path(tmp_dir, "lib_new")
dir.create(lib_old, showWarnings = FALSE)
dir.create(lib_new, showWarnings = FALSE)

# --- 1. Install CURRENT Version ---
cat("\n>>> Installing CURRENT version to temporary library...\n")
pkg_path <- file.path(wd, "SeqExpMatch")
install_cmd_new <- sprintf("R CMD INSTALL %s -l %s", shQuote(pkg_path), shQuote(lib_new))
if (system(install_cmd_new) != 0) stop("Failed to install current version")

# --- 2. Checkout and Install OLD Version ---
cat("\n>>> cloning and checking out OLD version (", OLD_COMMIT, ")...\n")
# Clean up repo dir if exists
unlink(repo_dir, recursive = TRUE)
system(sprintf("git clone %s %s", shQuote(wd), shQuote(repo_dir)))

old_wd <- getwd()
setwd(repo_dir)
if (system(sprintf("git checkout %s", OLD_COMMIT)) != 0) stop("Failed to checkout old commit")
setwd(old_wd)

cat("\n>>> Patching OLD version for build...\n")
pkg_path_old <- file.path(repo_dir, "SeqExpMatch")
src_path_old <- file.path(pkg_path_old, "src")

# Patch 0: DELETE ALL PRE-COMPILED BINARIES (Crucial fix for architecture mismatch)
cat("Cleaning up old binaries and Windows-specific files...\n")
bin_files <- list.files(src_path_old, pattern = "\\.(o|so|dll|exe|a|lib|exp|def)$", full.names = TRUE, ignore.case = TRUE)
if (length(bin_files) > 0) {
    cat(sprintf("Removing %d binary files...\n", length(bin_files)))
    unlink(bin_files)
}
unlink(file.path(src_path_old, "Makevars.win"))

# Patch 1: Copy Makevars from current (to match system build env)
if (file.exists(file.path(pkg_path, "src/Makevars"))) {
    cat("Copying Makevars from current version...\n")
    file.copy(file.path(pkg_path, "src/Makevars"), file.path(pkg_path_old, "src/Makevars"), overwrite = TRUE)
}

# Patch 2: Run compileAttributes to regenerate exports
cat("Running Rcpp::compileAttributes on old version...\n")
tryCatch({
    library(Rcpp)
    compileAttributes(pkg_path_old)
}, error = function(e) warning("compileAttributes failed: ", e$message))

cat("\n>>> Installing OLD version to temporary library...\n")
install_cmd_old <- sprintf("R CMD INSTALL %s -l %s", shQuote(pkg_path_old), shQuote(lib_old))
if (system(install_cmd_old) != 0) {
    cat("\n!!! Failed to install OLD version. Proceeding with NEW version only benchmark. !!!\n")
    SKIP_OLD <- TRUE
} else {
    SKIP_OLD <- FALSE
}

# --- 3. Run Benchmarks ---
res_file_old <- file.path(tmp_dir, "res_old.rds")
res_file_new <- file.path(tmp_dir, "res_new.rds")

if (!SKIP_OLD) {
    cat(sprintf("\n>>> Running Benchmark on OLD version (%d iterations)...\n", N_ITER))
    cmd_old <- sprintf("Rscript %s %s %s %d %s", WORKER_SCRIPT, shQuote(lib_old), shQuote(res_file_old), N_ITER, "old")
    if (system(cmd_old) != 0) {
        warning("Benchmark worker failed for OLD version")
        SKIP_OLD <- TRUE
    }
}

cat(sprintf("\n>>> Running Benchmark on NEW version (%d iterations)...\n", N_ITER))
cmd_new <- sprintf("Rscript %s %s %s %d %s", WORKER_SCRIPT, shQuote(lib_new), shQuote(res_file_new), N_ITER, "new")
if (system(cmd_new) != 0) stop("Benchmark worker failed for NEW version")

# --- 4. Analyze Results ---
summarize_results <- function(dt) {
    if (nrow(dt) == 0) {
        return(data.table())
    }
    dt[,
       .(
           mean_time = mean(elapsed, na.rm = TRUE),
           median_time = median(elapsed, na.rm = TRUE),
           sd_time = sd(elapsed, na.rm = TRUE),
           error_rate = sum(!is.na(error)) / .N,
           successes = sum(status == "success"),
           total = .N
       ),
       by = .(response_type, design, inference, metric, version_label)
    ]
}

print_top_speedups <- function(comparison_table, n_rows = 10) {
    if (nrow(comparison_table) == 0) {
        cat("No comparable metrics were recorded between versions.\n")
        return()
    }
    comparison_table[, speedup := ifelse(mean_time_new > 0,
        mean_time_old / mean_time_new,
        NA_real_
    )]
    comparison_table[!is.finite(speedup), speedup := NA_real_]
    top_rows <- comparison_table[order(-speedup)][1:min(n_rows, .N)]
    if (nrow(top_rows) == 0) {
        cat("No valid speedup entries to display.\n")
        return()
    }
    print(top_rows[, .(
        response_type,
        design,
        inference,
        metric,
        mean_time_old,
        mean_time_new,
        speedup,
        error_rate_old,
        error_rate_new
    )])
}

cat("\n>>> Analysis <<<\n")
if (file.exists(res_file_new)) {
    res_new <- readRDS(res_file_new)

    if (!SKIP_OLD && file.exists(res_file_old)) {
        res_old <- readRDS(res_file_old)

        summary_old <- summarize_results(res_old)
        summary_new <- summarize_results(res_new)

        summary_old_core <- summary_old[, setdiff(names(summary_old), "version_label"), with = FALSE]
        summary_new_core <- summary_new[, setdiff(names(summary_new), "version_label"), with = FALSE]
        comparison <- merge(
            summary_old_core,
            summary_new_core,
            by = c("response_type", "design", "inference", "metric"),
            suffixes = c("_old", "_new")
        )

        cat("\n--- Comparative timing per metric (top speedups) ---\n")
        print_top_speedups(comparison, n_rows = 20)

        cat("\n--- Specific comparison for randomization p-value (compute_two_sided_pval_for_treatment_effect_rand) ---\n")
        rand_comp <- comparison[metric == "compute_two_sided_pval_for_treatment_effect_rand"]
        if (nrow(rand_comp) > 0) {
            print(rand_comp[order(-speedup), .(
                response_type,
                design,
                inference,
                mean_time_old,
                mean_time_new,
                speedup
            )])
        } else {
            cat("Metric not found in both versions.\n")
        }
    } else {
        cat("Old version results not available. Showing New version timing summary only.\n")
        summary_new <- summarize_results(res_new)
        cat("\n--- New version timing summary (slowest metrics) ---\n")
        print(summary_new[order(-mean_time)][1:min(20, .N)][, .(
            response_type,
            design,
            inference,
            metric,
            mean_time,
            median_time,
            error_rate
        )])
    }
}

cat("\nDone.\n")
