test_that("univariate logit reusable bootstrap worker matches generic bootstrap path", {
	SlowInferenceIncidUnivLogRegr = R6::R6Class(
		"SlowInferenceIncidUnivLogRegr",
		inherit = InferenceIncidUnivLogRegr,
		private = list(
			supports_reusable_bootstrap_worker = function(){
				FALSE
			}
		)
	)

	set.seed(20260407)
	n = 52
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))

	des = FixedDesignBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	y = stats::rbinom(n, 1, stats::plogis(-0.4 + 0.8 * w))
	des$add_all_subject_responses(y)

	fast_inf = InferenceIncidUnivLogRegr$new(des, verbose = FALSE)
	slow_inf = SlowInferenceIncidUnivLogRegr$new(des, verbose = FALSE)
	fast_inf$num_cores = 1L
	slow_inf$num_cores = 1L

	set.seed(123)
	fast_boot = fast_inf$approximate_bootstrap_distribution_beta_hat_T(B = 21, show_progress = FALSE)
	set.seed(123)
	slow_boot = slow_inf$approximate_bootstrap_distribution_beta_hat_T(B = 21, show_progress = FALSE)

	expect_equal(fast_boot, slow_boot, tolerance = 1e-12)
})

test_that("multivariate logit reusable bootstrap worker matches generic bootstrap path", {
	SlowInferenceIncidMultiLogRegr = R6::R6Class(
		"SlowInferenceIncidMultiLogRegr",
		inherit = InferenceIncidMultiLogRegr,
		private = list(
			supports_reusable_bootstrap_worker = function(){
				FALSE
			}
		)
	)

	set.seed(20260408)
	n = 56
	p = 4
	X = as.data.frame(matrix(rnorm(n * p), nrow = n, ncol = p))
	colnames(X) = paste0("x", seq_len(p))

	des = FixedDesignBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	linpred = -0.3 + 0.7 * w + 0.5 * X$x1 - 0.3 * X$x2 + 0.1 * X$x3
	y = stats::rbinom(n, 1, stats::plogis(linpred))
	des$add_all_subject_responses(y)

	fast_inf = InferenceIncidMultiLogRegr$new(des, verbose = FALSE)
	slow_inf = SlowInferenceIncidMultiLogRegr$new(des, verbose = FALSE)
	fast_inf$num_cores = 1L
	slow_inf$num_cores = 1L

	set.seed(456)
	fast_boot = fast_inf$approximate_bootstrap_distribution_beta_hat_T(B = 21, show_progress = FALSE)
	set.seed(456)
	slow_boot = slow_inf$approximate_bootstrap_distribution_beta_hat_T(B = 21, show_progress = FALSE)

	expect_equal(fast_boot, slow_boot, tolerance = 1e-12)
})
