library(SeqExpMatch)

set.seed(123)

# Test proportion response with delta=0.5 (should error, not return NA)
seq_des = SeqDesignCRD$new(n = 30, response_type = 'proportion')

for (i in 1:30) {
  seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[i, 2:10])
}

# Generate some proportion data
proportions = runif(30, 0.1, 0.9)
seq_des$add_all_subject_responses(proportions)

seq_des_inf = SeqDesignInferencePropUniBetaRegr$new(seq_des)

cat("\n=== Testing proportion randomization with delta=0 (should work) ===\n")
pval_delta0 = seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(delta=0, nsim_exact_test=10)
cat("p-value with delta=0:", pval_delta0, "\n")

cat("\n=== Testing proportion randomization with delta=0.5 (should error) ===\n")
tryCatch({
  pval_delta05 = seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(delta=0.5, nsim_exact_test=10)
  cat("UNEXPECTED: Got p-value:", pval_delta05, "(should have errored)\n")
}, error = function(e) {
  cat("Expected error received:\n")
  cat("  Error:", e$message, "\n")
})

cat("\n=== Test complete ===\n")
