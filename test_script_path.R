#!/usr/bin/env Rscript
rm(list = ls())

cat("Testing script path detection...\n")

file_arg = "--file="
all_args = commandArgs(trailingOnly = FALSE)
cat("All commandArgs:", paste(all_args, collapse = ", "), "\n")

grep_result = grep(file_arg, all_args)
cat("grep result:", paste(grep_result, collapse = ", "), "\n")
cat("grep result length:", length(grep_result), "\n")

if (length(grep_result) > 0) {
  script_path = sub(file_arg, "", all_args[grep_result])
  cat("script_path:", script_path, "\n")

  test_dir = dirname(normalizePath(script_path))
  cat("test_dir:", test_dir, "\n")

  lib_path = normalizePath(file.path(test_dir, "..", "Rlib"), mustWork = FALSE)
  cat("lib_path:", lib_path, "\n")
  cat("lib_path exists:", dir.exists(lib_path), "\n")
} else {
  cat("ERROR: No --file= argument found\n")
  cat("This can happen when running interactively or with certain R execution modes\n")
}

cat("Now trying to load SeqExpMatch...\n")
.libPaths(c("Rlib", .libPaths()))
library(SeqExpMatch)
cat("SUCCESS!\n")
