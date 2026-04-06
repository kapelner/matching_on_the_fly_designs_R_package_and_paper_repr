test_that("Modified Poisson reusable bootstrap worker matches generic bootstrap path", {
	SlowInferenceIncidMultiModifiedPoisson = R6::R6Class(
		"SlowInferenceIncidMultiModifiedPoisson",
		inherit = InferenceIncidMultiModifiedPoisson,
		private = list(
			supports_reusable_bootstrap_worker = function(){
				FALSE
			}
		)
	)

	set.seed(20260406)
	n = 48
	p = 4
	X = as.data.frame(matrix(rnorm(n * p), nrow = n, ncol = p))
	colnames(X) = paste0("x", seq_len(p))

	des = FixedDesignBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	linpred = -0.5 + 0.7 * w + 0.4 * X$x1 - 0.2 * X$x2
	y = stats::rbinom(n, 1, stats::plogis(linpred))
	des$add_all_subject_responses(y)

	fast_inf = InferenceIncidMultiModifiedPoisson$new(des, verbose = FALSE)
	slow_inf = SlowInferenceIncidMultiModifiedPoisson$new(des, verbose = FALSE)
	fast_inf$num_cores = 1L
	slow_inf$num_cores = 1L

	set.seed(99)
	fast_boot = fast_inf$approximate_bootstrap_distribution_beta_hat_T(B = 21, show_progress = FALSE)
	set.seed(99)
	slow_boot = slow_inf$approximate_bootstrap_distribution_beta_hat_T(B = 21, show_progress = FALSE)

	expect_equal(fast_boot, slow_boot, tolerance = 1e-12)
})
