library(SeqExpMatch)

set.seed(123)

# Create KK design
seq_des = SeqDesignKK14$new(n = 30, response_type = "continuous")

for (i in 1:30) {
  seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[i, 2:10])
}

# Add continuous responses
responses = rnorm(30, mean = 5, sd = 2)
seq_des$add_all_subject_responses(responses)

# Create inference object
seq_des_inf = SeqDesignInferenceContinMultOLSKK$new(seq_des)

cat("\n== Inference: SeqDesignInferenceContinMultOLSKK\n\n")

cat("Calling compute_treatment_estimate()\n")
est = seq_des_inf$compute_treatment_estimate()
print(est)

cat("Calling compute_mle_confidence_interval()\n")
mle_ci = seq_des_inf$compute_mle_confidence_interval()
print(mle_ci)

cat("Calling compute_bootstrap_confidence_interval()\n")
boot_ci = seq_des_inf$compute_bootstrap_confidence_interval(B = 100)
print(boot_ci)

cat("Calling compute_confidence_interval_rand() [THIS MAY TAKE A WHILE]\n")
rand_ci = seq_des_inf$compute_confidence_interval_rand(alpha = 0.05, nsim_exact_test = 50)
print(rand_ci)

cat("\n== Summary ==\n")
cat(sprintf("Estimate:      %8.3f\n", est))
cat(sprintf("MLE CI:        [%8.3f, %8.3f] (width: %.3f)\n",
            mle_ci[1], mle_ci[2], mle_ci[2] - mle_ci[1]))
cat(sprintf("Bootstrap CI:  [%8.3f, %8.3f] (width: %.3f)\n",
            boot_ci[1], boot_ci[2], boot_ci[2] - boot_ci[1]))
cat(sprintf("Rand CI:       [%8.3f, %8.3f] (width: %.3f)\n",
            rand_ci[1], rand_ci[2], rand_ci[2] - rand_ci[1]))
