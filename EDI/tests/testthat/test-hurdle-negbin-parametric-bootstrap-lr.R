library(testthat)
library(EDI)

make_param_boot_hurdle_negbin_design <- function(seed = 20260519L, n = 140L){
	set.seed(seed)
	x1 <- rnorm(n)
	x2 <- rnorm(n)
	w <- rep(c(1, 0), length.out = n)
	p_pos <- plogis(-0.35 + 0.75 * w + 0.30 * x1 - 0.20 * x2)
	mu <- exp(0.45 + 0.30 * w + 0.20 * x1 - 0.15 * x2)
	theta <- 2.5
	u_pos <- rbinom(n, 1L, p_pos)
	cdf0 <- pnbinom(0, size = theta, mu = mu)
	u_trunc <- cdf0 + (1 - cdf0) * runif(n)
	y_pos <- pmax(1L, as.integer(qnbinom(u_trunc, size = theta, mu = mu)))
	y <- ifelse(u_pos == 1L, y_pos, 0L)

	des <- DesignFixedBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1, x2 = x2))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)
	des
}

test_that("InferenceCountHurdleNegBin parametric bootstrap LR completes with low attrition", {
	des <- make_param_boot_hurdle_negbin_design()

	inf <- InferenceCountHurdleNegBin$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	priv <- inf$.__enclos_env__$private

	expect_true(isTRUE(priv$supports_lik_ratio_param_bootstrap()))

	for (seed in c(4101L, 4102L, 4103L)) {
		inf$set_seed(seed)
		inf$num_cores <- 1L
		p_boot <- inf$compute_lik_ratio_bootstrap_two_sided_pval(
			delta = 0,
			B = 21L,
			show_progress = FALSE,
			min_number_usable_samples = 15L,
			max_attempts_per_replicate = 2L
		)
		diag <- inf$get_last_param_bootstrap_diagnostics()

		expect_true(is.finite(p_boot))
		expect_gte(p_boot, 0)
		expect_lte(p_boot, 1)
		expect_true(is.list(diag))
		expect_equal(diag$B, 21L)
		expect_equal(diag$n_success + diag$n_failure, diag$B)
		expect_equal(sum(unlist(diag$reason_counts, use.names = FALSE)), diag$B)
		expect_gte(diag$n_success, 15L)
		expect_gte(diag$success_fraction, 0.70)
		expect_lte(diag$prop_full_refit_failure, 0.25)
		expect_lte(diag$prop_null_refit_failure, 0.25)
		expect_lte(diag$prop_non_finite_lr, 0.10)
	}
})

test_that("InferenceCountHurdleNegBin bootstrap LR confidence interval returns finite ordered bounds", {
	des <- make_param_boot_hurdle_negbin_design(seed = 20260520L, n = 120L)
	inf <- InferenceCountHurdleNegBin$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf$set_seed(5151)
	inf$num_cores <- 1L

	ci <- inf$compute_lik_ratio_bootstrap_confidence_interval(
		alpha = 0.20,
		B = 9L,
		show_progress = FALSE,
		min_number_usable_samples = 5L,
		max_attempts_per_replicate = 2L
	)
	diag <- inf$get_last_param_bootstrap_diagnostics()

	expect_true(is.numeric(ci))
	expect_length(ci, 2L)
	expect_true(all(is.finite(ci)))
	expect_lte(ci[1], ci[2])
	expect_true(is.list(diag))
	expect_gte(diag$n_success, 5L)
})
