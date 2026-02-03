library(SeqExpMatch)

set.seed(789)

cat("Creating design...\n")
seq_des = SeqDesignCRD$new(n = 30, response_type = "survival")

cat("Adding subjects...\n")
for (i in 1:30) {
  seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[i, 2:10])
}

cat("Generating survival data...\n")
time = rexp(30, rate = 0.1)
dead = rbinom(30, 1, 0.8)

cat("Adding responses...\n")
seq_des$add_all_subject_responses_and_censoring(time, dead)

cat("Creating inference object...\n")
seq_des_inf = SeqDesignInferenceSurvivalUniWeibullRegr$new(seq_des, verbose = TRUE)

cat("\nAttempting to compute treatment estimate...\n")
tryCatch({
  est = seq_des_inf$compute_treatment_estimate()
  cat("Estimate:", est, "\n")
}, error = function(e) {
  cat("ERROR in compute_treatment_estimate:\n")
  cat("  Message:", e$message, "\n")
  cat("  Call:", deparse(e$call), "\n")
  traceback()
})
