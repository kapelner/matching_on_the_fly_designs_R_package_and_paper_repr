#!/usr/bin/env Rscript
cat("Step 1: Before loading packages\n")
rm(list = ls())

cat("Step 2: Loading pacman packages\n")
pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)

cat("Step 3: Setting lib paths\n")
.libPaths(c("Rlib", .libPaths()))

cat("Step 4: Loading SeqExpMatch\n")
library(SeqExpMatch)

cat("Step 5: After loading SeqExpMatch\n")
cat("Step 6: Setting seed\n")
set.seed(1)

cat("Step 7: Creating simple data\n")
n <- 20
X <- as.data.frame(matrix(rnorm(n * 3), ncol = 3))
colnames(X) <- paste0("x", 1:3)
y <- rnorm(n)

cat("Step 8: Creating design object\n")
seq_des_obj <- SeqDesignKK21$new(response_type = "continuous", n = n)

cat("Step 9: Adding subjects\n")
for (t in 1:n) {
  seq_des_obj$add_subject_to_experiment_and_assign(X[t, ])
  seq_des_obj$add_subject_response(t, y[t], 1)
}

cat("Step 10: Creating inference object\n")
seq_des_inf <- SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj)

cat("Step 11: Computing treatment estimate\n")
result <- seq_des_inf$compute_treatment_estimate()

cat("SUCCESS: Treatment estimate =", result, "\n")
