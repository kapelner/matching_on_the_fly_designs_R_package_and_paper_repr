#!/usr/bin/env Rscript
pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
library(SeqExpMatch)

set.seed(1)
max_n_dataset = 500
source("package_tests/_dataset_load.R")

D = datasets_and_response_models$boston
n = min(nrow(D$X), 200)
y = D$y_original$continuous
dead = rep(1, n)

cat("Creating SeqDesignKK21 object\n")
seq_des_obj = SeqDesignKK21$new(response_type = "continuous", n = n)

cat("Adding subjects...\n")
for (t in 1 : n){
  seq_des_obj$add_subject_to_experiment_and_assign(D$X[t, ])
  seq_des_obj$add_subject_response(t, y[t], dead[t])
}

cat("Creating SeqDesignInferenceAllKKCompoundMeanDiff...\n")
inf2 = SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj)

cat("Computing treatment estimate...\n")
te = inf2$compute_treatment_estimate()
cat("Treatment estimate:", te, "\n")

cat("\nCalling approximate_bootstrap_distribution_beta_hat_T with B=10...\n")
result = inf2$approximate_bootstrap_distribution_beta_hat_T(B = 10)
cat("Result class:", class(result), "\n")
cat("Result type:", typeof(result), "\n")
cat("Result length:", length(result), "\n")
cat("Result:\n")
print(result)

cat("\nTrying to call quantile on result...\n")
q = quantile(result, c(0.025, 0.975))
cat("Quantiles:", q, "\n")

cat("\nSUCCESS!\n")
