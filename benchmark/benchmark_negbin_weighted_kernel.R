# Before/after benchmark for the NegBin bootstrap-weighted refit path.
# "Before": MASS::glm.nb(weights = ...) — the engine InferenceCountNegBin$
#   compute_estimate_with_bootstrap_weights() used prior to this change.
# "After":  fast_neg_bin_weighted_cpp(...) — the new EDI C++ kernel now wired
#   into that same method.
library(EDI)

bench_one_n <- function(n, reps = 20L) {
	set.seed(20260728)
	X <- cbind(`(Intercept)` = 1, w = rbinom(n, 1L, 0.5), x1 = rnorm(n), x2 = rnorm(n))
	beta_true <- c(0.5, -0.3, 0.2, -0.1)
	theta_true <- 4
	mu <- exp(X %*% beta_true)
	y <- rnbinom(n, size = theta_true, mu = as.numeric(mu))
	weights <- sample(c(0.5, 1, 1.5, 2), n, replace = TRUE)
	df <- data.frame(y = y, w = X[, "w"], x1 = X[, "x1"], x2 = X[, "x2"])

	t_mass <- system.time(for (i in seq_len(reps)) {
		suppressWarnings(MASS::glm.nb(y ~ ., data = df, weights = weights))
	})[["elapsed"]]

	t_cpp <- system.time(for (i in seq_len(reps)) {
		fast_neg_bin_weighted_cpp(X = X, y = as.integer(y), weights = as.numeric(weights), smart_cold_start = TRUE)
	})[["elapsed"]]

	data.frame(n = n, reps = reps, mass_glm_nb_sec = t_mass, fast_neg_bin_weighted_cpp_sec = t_cpp, speedup = t_mass / t_cpp)
}

results <- do.call(rbind, lapply(c(50L, 200L, 1000L, 5000L), bench_one_n))
print(results, row.names = FALSE)
