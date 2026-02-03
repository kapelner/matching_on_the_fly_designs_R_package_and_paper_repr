rm(list = ls())
pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
.libPaths(c("Rlib", .libPaths()))
library(SeqExpMatch)
set.seed(1)
max_n_dataset = 500
source("package_tests/_dataset_load.R")

D = datasets_and_response_models$boston
n = min(nrow(D$X), 200)
y = D$y_original$continuous
dead = rep(1, n)

cat("Creating KK21 design...\n")
seq_des_obj = SeqDesignKK21$new(response_type = "continuous", n = n)
for (t in 1:n) {
  seq_des_obj$add_subject_to_experiment_and_assign(D$X[t, ])
  seq_des_obj$add_subject_response(t, y[t], dead[t])
}
cat("Created design with", n, "subjects\n")

cat("Creating AllSimpleMeanDiff...\n")
seq_des_inf1 = SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj)
cat("Computing treatment estimate...\n")
seq_des_inf1$compute_treatment_estimate()
cat("Treatment estimate:", seq_des_inf1$compute_treatment_estimate(), "\n")

cat("Computing MLE pval...\n")
seq_des_inf1$compute_mle_two_sided_pval_for_treatment_effect()

cat("Computing MLE CI...\n")
seq_des_inf1$compute_mle_confidence_interval(0.05)

cat("Computing bootstrap CI (B=10)...\n")
seq_des_inf1$compute_bootstrap_confidence_interval(B = 10, na.rm = TRUE)

cat("All tests passed!\n")
