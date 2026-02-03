#!/usr/bin/env Rscript
pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
library(SeqExpMatch)

set.seed(1)
max_n_dataset = 500
source("package_tests/_dataset_load.R")

D = datasets_and_response_models$boston
n = min(nrow(D$X), 200)

cat("Testing all response types with basic inference...\n\n")

# Continuous
cat("1. CONTINUOUS response: ")
y = D$y_original$continuous
seq_des = SeqDesignKK21$new(response_type = "continuous", n = n)
for (t in 1:n) {
  seq_des$add_subject_to_experiment_and_assign(D$X[t, ])
  seq_des$add_subject_response(t, y[t], 1)
}
inf = SeqDesignInferenceAllSimpleMeanDiff$new(seq_des)
te = inf$compute_treatment_estimate()
cat("✓ Treatment estimate =", te, "\n")

# Incidence
cat("2. INCIDENCE response: ")
y = D$y_original$incidence
seq_des = SeqDesignKK21$new(response_type = "incidence", n = n)
for (t in 1:n) {
  seq_des$add_subject_to_experiment_and_assign(D$X[t, ])
  seq_des$add_subject_response(t, y[t], 1)
}
inf = SeqDesignInferenceIncidUnivLogRegr$new(seq_des)
te = inf$compute_treatment_estimate()
cat("✓ Treatment estimate =", te, "\n")

# Proportion
cat("3. PROPORTION response: ")
y = D$y_original$proportion
seq_des = SeqDesignKK21$new(response_type = "proportion", n = n)
for (t in 1:n) {
  seq_des$add_subject_to_experiment_and_assign(D$X[t, ])
  seq_des$add_subject_response(t, y[t], 1)
}
inf = SeqDesignInferencePropUniBetaRegr$new(seq_des)
te = inf$compute_treatment_estimate()
cat("✓ Treatment estimate =", te, "\n")

# Count
cat("4. COUNT response: ")
y = D$y_original$count
seq_des = SeqDesignKK21$new(response_type = "count", n = n)
for (t in 1:n) {
  seq_des$add_subject_to_experiment_and_assign(D$X[t, ])
  seq_des$add_subject_response(t, y[t], 1)
}
inf = SeqDesignInferenceCountUnivNegBinRegr$new(seq_des)
te = inf$compute_treatment_estimate()
cat("✓ Treatment estimate =", te, "\n")

# Survival
cat("5. SURVIVAL response: ")
y = D$y_original$survival
seq_des = SeqDesignKK21$new(response_type = "survival", n = n)
for (t in 1:n) {
  seq_des$add_subject_to_experiment_and_assign(D$X[t, ])
  seq_des$add_subject_response(t, y[t], 1)
}
inf = SeqDesignInferenceSurvivalRestrictedMeanDiff$new(seq_des)
te = inf$compute_treatment_estimate()
cat("✓ Treatment estimate =", te, "\n")

cat("\n✓✓✓ ALL RESPONSE TYPES WORKING! ✓✓✓\n")
