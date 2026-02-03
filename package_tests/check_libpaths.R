#!/usr/bin/env Rscript
rm(list = ls())
pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
file_arg = "--file="
script_path = sub(file_arg, "", commandArgs(trailingOnly = FALSE)[grep(file_arg, commandArgs(trailingOnly = FALSE))])
test_dir = dirname(normalizePath(script_path))
lib_path = normalizePath(file.path(test_dir, "..", "Rlib"), mustWork = FALSE)

cat("script_path:", script_path, "\n")
cat("test_dir:", test_dir, "\n")
cat("lib_path:", lib_path, "\n")
cat("dir.exists(lib_path):", dir.exists(lib_path), "\n")

if (dir.exists(lib_path)){
  .libPaths(c(lib_path, .libPaths()))
}

cat("\n.libPaths after modification:\n")
print(.libPaths())

cat("\nLoading SeqExpMatch...\n")
library(SeqExpMatch)

cat("\nSeqExpMatch loaded from:\n")
print(system.file(package = "SeqExpMatch"))

cat("\nChecking if mean_cpp is exported:\n")
cat("mean_cpp exists:", exists("mean_cpp"), "\n")

if (exists("mean_cpp")) {
  cat("Testing mean_cpp:\n")
  result = mean_cpp(c(1, 2, 3, 4, 5))
  cat("mean_cpp(1:5) =", result, "\n")
}
