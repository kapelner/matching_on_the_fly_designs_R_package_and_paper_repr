#' Benchmark Worker Script
#' This script is called by the main benchmark_comparison.R script.
#' It expects 3 arguments: library_path, output_file, n_iter

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Not enough arguments")

lib_path <- args[1]
out_file <- args[2]
n_iter <- as.integer(args[3])
version_label <- if (length(args) >= 4) args[4] else "unknown"

# Load library from the specific path
.libPaths(c(lib_path, .libPaths()))
helper_path <- file.path(getwd(), "benchmark_inference_loops.R")
source(helper_path)

tryCatch({
    library(SeqExpMatch)
}, error = function(e) {
    stop("Failed to load SeqExpMatch library: ", e$message)
})

cat(sprintf("Running benchmark in library: %s\n", lib_path))
cat(sprintf("Package Version: %s\n", packageVersion("SeqExpMatch")))

# --- Execution Loop Using Comprehensive Inference Suite ---
res_list <- vector("list", n_iter)
cat("Starting comprehensive inference benchmark iterations...\n")
for (i in seq_len(n_iter)) {
  seed <- 12345 + i
  res_list[[i]] <- run_comprehensive_inference_suite(
    version_label = version_label,
    iteration = i,
    n = 100,
    p = 5,
    seed = seed,
    verbose = FALSE
  )
}

benchmark_results <- data.table::rbindlist(res_list, fill = TRUE)
saveRDS(benchmark_results, file = out_file)
cat("Benchmark worker finished.\n")
