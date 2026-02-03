library(SeqExpMatch)

set.seed(123)

# Create a simple KK design
seq_des = SeqDesignKK14$new(n = 30, response_type = "continuous")

for (i in 1:30) {
  seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[i, 2:10])
}

# Add continuous responses
responses = rnorm(30, mean = 5, sd = 2)
seq_des$add_all_subject_responses(responses)

# Create inference object with verbose mode
seq_des_inf = SeqDesignInferenceContinMultOLSKK$new(seq_des, verbose = TRUE)

cat("\n=== Treatment Estimate ===\n")
est = seq_des_inf$compute_treatment_estimate()
cat("Estimate:", est, "\n")

cat("\n=== Randomization Distribution (first 10 draws) ===\n")
t0s = seq_des_inf$compute_beta_hat_T_randomization_distr_under_sharp_null(
  nsim_exact_test = 10,
  delta = 0,
  transform_responses = "none",
  show_progress = FALSE
)
cat("t0s:", t0s, "\n")
cat("Range:", range(t0s, na.rm = TRUE), "\n")
cat("NAs:", sum(is.na(t0s)), "\n")

cat("\n=== P-value Computation ===\n")
pval = seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test = 10)
cat("P-value:", pval, "\n")

cat("\n=== Checking comparison ===\n")
cat("Observed estimate:", est, "\n")
cat("Number of t0s >= est:", sum(t0s >= est, na.rm = TRUE), "\n")
cat("Number of t0s <= est:", sum(t0s <= est, na.rm = TRUE), "\n")
cat("Total valid t0s:", sum(!is.na(t0s)), "\n")
