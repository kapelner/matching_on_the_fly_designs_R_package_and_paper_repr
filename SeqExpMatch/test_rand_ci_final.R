library(SeqExpMatch)

set.seed(456)  # Different seed to test robustness

# Create KK design with same parameters as user's case
seq_des = SeqDesignKK14$new(n = 30, response_type = "continuous")

for (i in 1:30) {
  seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[i, 2:10])
}

# Add continuous responses
responses = rnorm(30, mean = 5, sd = 2)
seq_des$add_all_subject_responses(responses)

# Create inference object
seq_des_inf = SeqDesignInferenceContinMultOLSKK$new(seq_des)

cat("\n== Testing Randomization CI with More Simulations ==\n\n")

cat("Estimate:", seq_des_inf$compute_treatment_estimate(), "\n\n")

cat("MLE CI:\n")
print(seq_des_inf$compute_mle_confidence_interval())

cat("\nBootstrap CI (B=200):\n")
boot_ci = seq_des_inf$compute_bootstrap_confidence_interval(B = 200)
print(boot_ci)

cat("\nRandomization CI (nsim=200):\n")
rand_ci = seq_des_inf$compute_confidence_interval_rand(
  alpha = 0.05,
  nsim_exact_test = 200,
  pval_epsilon = 0.01
)
print(rand_ci)

cat("\n== Comparison ==\n")
cat(sprintf("Bootstrap CI width:      %.2f\n", boot_ci[2] - boot_ci[1]))
cat(sprintf("Randomization CI width:  %.2f\n", rand_ci[2] - rand_ci[1]))
cat(sprintf("Ratio (Rand/Boot):       %.2f\n",
            (rand_ci[2] - rand_ci[1]) / (boot_ci[2] - boot_ci[1])))
