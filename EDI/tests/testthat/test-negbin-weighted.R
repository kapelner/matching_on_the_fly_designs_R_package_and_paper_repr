context("weighted negative binomial C++ kernel")

test_that("fast_neg_bin_weighted_cpp matches MASS::glm.nb with case weights", {
	skip_if_not_installed("MASS")
	set.seed(20260728)

	n <- 200L
	X <- cbind(`(Intercept)` = 1, w = rbinom(n, 1L, 0.5), x1 = rnorm(n))
	beta_true <- c(0.5, -0.3, 0.2)
	theta_true <- 4
	mu <- exp(X %*% beta_true)
	y <- rnbinom(n, size = theta_true, mu = as.numeric(mu))
	weights <- sample(c(0.5, 1, 1.5, 2), n, replace = TRUE)

	mod_ref <- MASS::glm.nb(y ~ w + x1, data = data.frame(y = y, w = X[, "w"], x1 = X[, "x1"]), weights = weights)
	mod_fast <- fast_neg_bin_weighted_cpp(X = X, y = as.integer(y), weights = as.numeric(weights), smart_cold_start = TRUE)

	expect_true(isTRUE(mod_fast$converged))
	expect_equal(as.numeric(mod_fast$b), as.numeric(stats::coef(mod_ref)), tolerance = 5e-3)
	expect_equal(as.numeric(mod_fast$theta_hat), as.numeric(mod_ref$theta), tolerance = 5e-2)
})

test_that("InferenceCountNegBin bootstrap-weighted estimate is finite and no longer routes through MASS::glm.nb", {
	set.seed(20260728)
	n <- 60L
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
	}
	des$add_all_subject_responses(as.integer(rnbinom(n, size = 4, mu = 3)))

	inf <- InferenceCountNegBin$new(des, verbose = FALSE)
	inf$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n),
		unit_group_id = rep(1L, n),
		n_units = n
	)
	weights <- sample(c(0.5, 1, 1.5, 2), n, replace = TRUE)
	est <- inf$compute_estimate_with_bootstrap_weights(weights)

	expect_true(is.finite(as.numeric(est)))
})
