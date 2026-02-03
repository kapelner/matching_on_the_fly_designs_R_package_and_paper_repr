#!/usr/bin/env Rscript
rm(list = ls())

cat("Loading packages...\n")
pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)

cat("Setting libpaths...\n")
.libPaths(c("Rlib", .libPaths()))

cat("Loading SeqExpMatch...\n")
library(SeqExpMatch)

cat("After loading SeqExpMatch\n")
cat("Setting seed and max_n_dataset...\n")
set.seed(1)
max_n_dataset = 500

cat("Sourcing _dataset_load.R...\n")
source("package_tests/_dataset_load.R")

cat("SUCCESS: Dataset loading complete\n")
