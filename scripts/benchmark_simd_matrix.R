args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(reps = 5L, out_dir = "benchmark/simd_matrix", keep_builds = FALSE)
  for (arg in args) {
    if (startsWith(arg, "--reps=")) out$reps <- as.integer(sub("^--reps=", "", arg))
    if (startsWith(arg, "--out_dir=")) out$out_dir <- sub("^--out_dir=", "", arg)
    if (identical(arg, "--keep_builds")) out$keep_builds <- TRUE
  }
  out
}

cfg <- parse_args(args)
repo_root <- normalizePath(getwd())
out_dir <- normalizePath(cfg$out_dir, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

build_root <- file.path(out_dir, "builds")
dir.create(build_root, recursive = TRUE, showWarnings = FALSE)

head_src_root <- file.path(build_root, "head_src")
dir.create(head_src_root, recursive = TRUE, showWarnings = FALSE)

archive_tar <- file.path(build_root, "head_snapshot.tar")
status <- system2("git", c("-C", repo_root, "archive", "--format=tar", paste0("--output=", archive_tar), "HEAD"))
if (!identical(status, 0L)) stop("failed to create HEAD snapshot archive")
utils::untar(archive_tar, exdir = head_src_root)

install_build <- function(label, source_root, env = character()) {
  lib <- file.path(build_root, paste0("lib_", label))
  unlink(lib, recursive = TRUE, force = TRUE)
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  pkg_dir <- file.path(source_root, "EDI")
  cmd <- c("CMD", "INSTALL", "-l", normalizePath(lib, mustWork = TRUE), normalizePath(pkg_dir))
  status <- system2("R", cmd, env = env)
  if (!identical(status, 0L)) stop("install failed for ", label)
  lib
}

run_worker <- function(label, lib) {
  out_file <- file.path(out_dir, paste0(label, ".csv"))
  cmd <- c(
    "scripts/benchmark_simd_matrix_worker.R",
    paste0("--lib=", normalizePath(lib)),
    paste0("--out=", normalizePath(out_file, mustWork = FALSE)),
    paste0("--reps=", cfg$reps)
  )
  status <- system2("Rscript", cmd)
  if (!identical(status, 0L)) stop("worker failed for ", label)
  df <- read.csv(out_file, stringsAsFactors = FALSE)
  df$build <- label
  df
}

builds <- list(
  current_normal = list(source = repo_root, env = character()),
  current_no_vectorization = list(source = repo_root, env = c("EDI_PORTABLE=1", "EDI_DISABLE_VECTORIZATION=1")),
  head_pre_optimization = list(source = head_src_root, env = character())
)

results <- list()
for (label in names(builds)) {
  spec <- builds[[label]]
  lib <- install_build(label, spec$source, spec$env)
  results[[label]] <- run_worker(label, lib)
}

raw <- do.call(rbind, results)
summary <- aggregate(
  elapsed ~ build + kernel,
  data = raw,
  FUN = function(x) c(median = median(x), mean = mean(x), min = min(x), max = max(x))
)
summary_expanded <- data.frame(
  build = summary$build,
  kernel = summary$kernel,
  median_elapsed = vapply(summary$elapsed, `[[`, numeric(1), "median"),
  mean_elapsed = vapply(summary$elapsed, `[[`, numeric(1), "mean"),
  min_elapsed = vapply(summary$elapsed, `[[`, numeric(1), "min"),
  max_elapsed = vapply(summary$elapsed, `[[`, numeric(1), "max"),
  stringsAsFactors = FALSE
)

baseline <- subset(summary_expanded, build == "head_pre_optimization", select = c("kernel", "median_elapsed"))
names(baseline)[2] <- "baseline_median_elapsed"
summary_expanded <- merge(summary_expanded, baseline, by = "kernel", all.x = TRUE, sort = FALSE)
summary_expanded$ratio_vs_head <- summary_expanded$median_elapsed / summary_expanded$baseline_median_elapsed

current_baseline <- subset(summary_expanded, build == "current_normal", select = c("kernel", "median_elapsed"))
names(current_baseline)[2] <- "current_normal_median_elapsed"
summary_expanded <- merge(summary_expanded, current_baseline, by = "kernel", all.x = TRUE, sort = FALSE)
summary_expanded$ratio_vs_current_normal <- summary_expanded$median_elapsed / summary_expanded$current_normal_median_elapsed

write.csv(raw, file.path(out_dir, "simd_matrix_raw.csv"), row.names = FALSE)
write.csv(summary_expanded[order(summary_expanded$kernel, summary_expanded$build), ],
          file.path(out_dir, "simd_matrix_summary.csv"), row.names = FALSE)

cat("Benchmark matrix written to:\n")
cat("  ", file.path(out_dir, "simd_matrix_raw.csv"), "\n", sep = "")
cat("  ", file.path(out_dir, "simd_matrix_summary.csv"), "\n", sep = "")

if (!cfg$keep_builds) {
  unlink(build_root, recursive = TRUE, force = TRUE)
}
