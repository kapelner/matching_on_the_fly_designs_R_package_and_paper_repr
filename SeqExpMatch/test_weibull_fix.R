library(SeqExpMatch)

set.seed(789)

# Create a simple CRD design with survival response
seq_des = SeqDesignCRD$new(n = 30, response_type = "survival")

for (i in 1:30) {
  seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[i, 2:10])
}

# Generate survival data (time and censoring indicator)
time = rexp(30, rate = 0.1)
dead = rbinom(30, 1, 0.8)  # 80% uncensored

seq_des$add_all_subject_responses_and_censoring(time, dead)

# Create Weibull inference object
seq_des_inf = SeqDesignInferenceSurvivalUniWeibullRegr$new(seq_des)

cat("\n== Inference: SeqDesignInferenceSurvivalUniWeibullRegr\n\n")

cat("Calling compute_treatment_estimate()\n")
tryCatch({
  est = seq_des_inf$compute_treatment_estimate()
  print(est)
}, error = function(e) {
  cat("Error:", e$message, "\n")
})

cat("\nCalling compute_mle_two_sided_pval_for_treatment_effect()\n")
tryCatch({
  pval = seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
  print(pval)
}, error = function(e) {
  cat("Error:", e$message, "\n")
})

cat("\nCalling compute_mle_confidence_interval()\n")
tryCatch({
  ci = seq_des_inf$compute_mle_confidence_interval()
  print(ci)
}, error = function(e) {
  cat("Error:", e$message, "\n")
})

cat("\nCalling compute_bootstrap_confidence_interval()\n")
tryCatch({
  boot_ci = seq_des_inf$compute_bootstrap_confidence_interval(B = 50)
  print(boot_ci)
}, error = function(e) {
  cat("Error:", e$message, "\n")
})

cat("\n== Test Complete ==\n")
