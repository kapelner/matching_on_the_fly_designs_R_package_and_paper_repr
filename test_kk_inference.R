.libPaths(c("Rlib", .libPaths()))
library(SeqExpMatch)

set.seed(1)
n <- 20
X <- as.data.frame(matrix(rnorm(n * 3), ncol = 3))
colnames(X) <- paste0("x", 1:3)
y <- rnorm(n)

cat("Creating DesignSeqOneByOneKK21...\n")
des_obj <- DesignSeqOneByOneKK21$new(response_type = "continuous", n = n)

cat("Adding subjects...\n")
for (t in 1:n) {
	des_obj$add_subject_to_experiment_and_assign(X[t, ])
	des_obj$add_subject_response(t, y[t], 1)
}

cat("Creating InferenceAllKKCompoundMeanDiff...\n")
tryCatch({
	seq_des_inf <- InferenceAllKKCompoundMeanDiff$new(des_obj)
	cat("SUCCESS: Created inference object\n")

	cat("Computing treatment estimate...\n")
	result <- seq_des_inf$compute_treatment_estimate()
	cat("Treatment estimate:", result, "\n")
}, error = function(e) {
	cat("ERROR:", e$message, "\n")
	traceback()
})
