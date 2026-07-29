# Before/after benchmark for KK robust-regression classes.
# "Before": use_rcpp = FALSE (MASS::rlm) — the only engine these classes used
#   prior to this change.
# "After":  use_rcpp = TRUE (fast_robust_regression_cpp) — now the default.
library(EDI)

make_kk_robust_design <- function(n_pairs, n_single, seed = 20260729) {
	set.seed(seed)
	n <- 2L * n_pairs + n_single
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = "continuous", verbose = FALSE)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
	}
	des$.__enclos_env__$private$m <- c(rep(seq_len(n_pairs), each = 2L), rep(0L, n_single))
	priv <- des$.__enclos_env__$private
	eta <- 0.4 - 0.5 * priv$w + 0.3 * priv$X[, 1] - 0.2 * priv$X[, 2]
	des$add_all_subject_responses(eta + rt(n, df = 4) * 0.5)
	des
}

bench_one <- function(n_pairs, n_single, reps = 10L, class_name = "InferenceContinKKRobustRegrOneLik") {
	des <- make_kk_robust_design(n_pairs, n_single)
	t_mass <- system.time(for (i in seq_len(reps)) {
		inf <- get(class_name)$new(des, method = "M", use_rcpp = FALSE, verbose = FALSE)
		inf$compute_estimate()
	})[["elapsed"]]
	t_cpp <- system.time(for (i in seq_len(reps)) {
		inf <- get(class_name)$new(des, method = "M", use_rcpp = TRUE, verbose = FALSE)
		inf$compute_estimate()
	})[["elapsed"]]
	data.frame(class = class_name, n = 2L * n_pairs + n_single, reps = reps,
		mass_rlm_sec = t_mass, fast_robust_regression_cpp_sec = t_cpp, speedup = t_mass / t_cpp)
}

results <- do.call(rbind, c(
	lapply(list(c(10L, 5L), c(40L, 20L), c(200L, 100L)), function(ns) bench_one(ns[1], ns[2], class_name = "InferenceContinKKRobustRegrOneLik")),
	lapply(list(c(10L, 5L), c(40L, 20L), c(200L, 100L)), function(ns) bench_one(ns[1], ns[2], class_name = "InferenceContinKKRobustRegrIVWC"))
))
print(results, row.names = FALSE)
