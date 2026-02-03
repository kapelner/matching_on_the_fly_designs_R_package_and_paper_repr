#!/usr/bin/env Rscript
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

cat("Creating SeqDesignKK21 object\n")
seq_des_obj = SeqDesignKK21$new(response_type = "continuous", n = n)

cat("Adding subjects...\n")
for (t in 1 : n){
  seq_des_obj$add_subject_to_experiment_and_assign(D$X[t, ])
  seq_des_obj$add_subject_response(t, y[t], dead[t])
}

cat("Subjects added successfully\n")

cat("Testing SeqDesignInferenceAllSimpleMeanDiff...\n")
inf1 = SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj)
cat("  Created successfully\n")
cat("  Treatment estimate:", inf1$compute_treatment_estimate(), "\n")

cat("\nNow testing SeqDesignInferenceAllKKCompoundMeanDiff...\n")
cat("  About to call $new()...\n")
flush.console()

# Try with R-level error handling
result = tryCatch({
  inf2 = SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj)
  cat("  Created successfully!\n")
  inf2
}, error = function(e) {
  cat("R ERROR:", e$message, "\n")
  traceback()
  NULL
})

if (!is.null(result)) {
  cat("  Treatment estimate:", result$compute_treatment_estimate(), "\n")
  cat("SUCCESS!\n")
} else {
  cat("FAILED - R error caught\n")
}
