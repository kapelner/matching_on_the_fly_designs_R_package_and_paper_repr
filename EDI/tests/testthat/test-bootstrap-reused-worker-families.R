compare_bootstrap_fast_slow <- function(fast_inf, slow_inf, B = 11L, seed = 1L, tolerance = 1e-10){
	fast_inf$num_cores = 1L
	slow_inf$num_cores = 1L
	set.seed(seed)
	fast_boot = fast_inf$approximate_bootstrap_distribution_beta_hat_T(B = B, show_progress = FALSE)
	set.seed(seed)
	slow_boot = slow_inf$approximate_bootstrap_distribution_beta_hat_T(B = B, show_progress = FALSE)
	expect_equal(unname(fast_boot), unname(slow_boot), tolerance = tolerance)
}

test_that("log-binomial reusable bootstrap worker matches generic path", {
	SlowInferenceIncidUnivLogBinomial = R6::R6Class(
		"SlowInferenceIncidUnivLogBinomial",
		inherit = InferenceIncidUnivLogBinomial,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidMultiLogBinomial = R6::R6Class(
		"SlowInferenceIncidMultiLogBinomial",
		inherit = InferenceIncidMultiLogBinomial,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260410)
	n = 60
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
	des = FixedDesignBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	p = pmin(0.95, pmax(0.02, exp(-1.4 + 0.5 * w + 0.2 * X$x1 - 0.1 * X$x2)))
	y = stats::rbinom(n, 1, p)
	des$add_all_subject_responses(y)

	compare_bootstrap_fast_slow(
		InferenceIncidUnivLogBinomial$new(des, verbose = FALSE),
		SlowInferenceIncidUnivLogBinomial$new(des, verbose = FALSE),
		seed = 101
	)
	compare_bootstrap_fast_slow(
		InferenceIncidMultiLogBinomial$new(des, verbose = FALSE),
		SlowInferenceIncidMultiLogBinomial$new(des, verbose = FALSE),
		seed = 102
	)
})

test_that("binomial identity reusable bootstrap worker matches generic path", {
	SlowInferenceIncidUnivBinomialIdentityRiskDiff = R6::R6Class(
		"SlowInferenceIncidUnivBinomialIdentityRiskDiff",
		inherit = InferenceIncidUnivBinomialIdentityRiskDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidMultiBinomialIdentityRiskDiff = R6::R6Class(
		"SlowInferenceIncidMultiBinomialIdentityRiskDiff",
		inherit = InferenceIncidMultiBinomialIdentityRiskDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260411)
	n = 64
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
	des = FixedDesignBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	p = pmin(0.95, pmax(0.02, 0.25 + 0.18 * w + 0.05 * X$x1 - 0.03 * X$x2))
	y = stats::rbinom(n, 1, p)
	des$add_all_subject_responses(y)

	compare_bootstrap_fast_slow(
		InferenceIncidUnivBinomialIdentityRiskDiff$new(des, verbose = FALSE),
		SlowInferenceIncidUnivBinomialIdentityRiskDiff$new(des, verbose = FALSE),
		seed = 103
	)
	compare_bootstrap_fast_slow(
		InferenceIncidMultiBinomialIdentityRiskDiff$new(des, verbose = FALSE),
		SlowInferenceIncidMultiBinomialIdentityRiskDiff$new(des, verbose = FALSE),
		seed = 104
	)
})

test_that("count reusable bootstrap workers match generic paths", {
	SlowInferenceCountUnivPoissonRegr = R6::R6Class(
		"SlowInferenceCountUnivPoissonRegr",
		inherit = InferenceCountUnivPoissonRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountMultiPoissonRegr = R6::R6Class(
		"SlowInferenceCountMultiPoissonRegr",
		inherit = InferenceCountMultiPoissonRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountUnivRobustPoissonRegr = R6::R6Class(
		"SlowInferenceCountUnivRobustPoissonRegr",
		inherit = InferenceCountUnivRobustPoissonRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountMultiRobustPoissonRegr = R6::R6Class(
		"SlowInferenceCountMultiRobustPoissonRegr",
		inherit = InferenceCountMultiRobustPoissonRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountUnivQuasiPoissonRegr = R6::R6Class(
		"SlowInferenceCountUnivQuasiPoissonRegr",
		inherit = InferenceCountUnivQuasiPoissonRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountMultiQuasiPoissonRegr = R6::R6Class(
		"SlowInferenceCountMultiQuasiPoissonRegr",
		inherit = InferenceCountMultiQuasiPoissonRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260412)
	n = 68
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
	des = FixedDesignBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	lambda = exp(0.2 + 0.35 * w + 0.15 * X$x1 - 0.1 * X$x2)
	y = stats::rpois(n, lambda)
	des$add_all_subject_responses(y)

	compare_bootstrap_fast_slow(InferenceCountUnivPoissonRegr$new(des, verbose = FALSE), SlowInferenceCountUnivPoissonRegr$new(des, verbose = FALSE), seed = 105)
	compare_bootstrap_fast_slow(InferenceCountMultiPoissonRegr$new(des, verbose = FALSE), SlowInferenceCountMultiPoissonRegr$new(des, verbose = FALSE), seed = 106)
	compare_bootstrap_fast_slow(InferenceCountUnivRobustPoissonRegr$new(des, verbose = FALSE), SlowInferenceCountUnivRobustPoissonRegr$new(des, verbose = FALSE), seed = 107)
	compare_bootstrap_fast_slow(InferenceCountMultiRobustPoissonRegr$new(des, verbose = FALSE), SlowInferenceCountMultiRobustPoissonRegr$new(des, verbose = FALSE), seed = 108)
	compare_bootstrap_fast_slow(InferenceCountUnivQuasiPoissonRegr$new(des, verbose = FALSE), SlowInferenceCountUnivQuasiPoissonRegr$new(des, verbose = FALSE), seed = 109)
	compare_bootstrap_fast_slow(InferenceCountMultiQuasiPoissonRegr$new(des, verbose = FALSE), SlowInferenceCountMultiQuasiPoissonRegr$new(des, verbose = FALSE), seed = 110)
})

test_that("continuous robust reusable bootstrap worker matches generic path", {
	SlowInferenceContinUnivRobustRegr = R6::R6Class(
		"SlowInferenceContinUnivRobustRegr",
		inherit = InferenceContinUnivRobustRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceContinMultiRobustRegr = R6::R6Class(
		"SlowInferenceContinMultiRobustRegr",
		inherit = InferenceContinMultiRobustRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260413)
	n = 60
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
	des = FixedDesignBernoulli$new(n = n, response_type = "continuous", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	y = 0.4 * w + 0.3 * X$x1 - 0.2 * X$x2 + stats::rt(n, df = 5)
	des$add_all_subject_responses(y)

	compare_bootstrap_fast_slow(InferenceContinUnivRobustRegr$new(des, method = "M", verbose = FALSE), SlowInferenceContinUnivRobustRegr$new(des, method = "M", verbose = FALSE), seed = 111)
	compare_bootstrap_fast_slow(InferenceContinMultiRobustRegr$new(des, method = "M", verbose = FALSE), SlowInferenceContinMultiRobustRegr$new(des, method = "M", verbose = FALSE), seed = 112)
})

test_that("continuous quantile reusable bootstrap worker matches generic path", {
	skip_if_not_installed("quantreg")

	SlowInferenceContinUnivQuantileRegr = R6::R6Class(
		"SlowInferenceContinUnivQuantileRegr",
		inherit = InferenceContinUnivQuantileRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceContinMultiQuantileRegr = R6::R6Class(
		"SlowInferenceContinMultiQuantileRegr",
		inherit = InferenceContinMultiQuantileRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260414)
	n = 56
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
	des = FixedDesignBernoulli$new(n = n, response_type = "continuous", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	y = 0.5 * w + 0.25 * X$x1 - 0.15 * X$x2 + stats::rt(n, df = 6)
	des$add_all_subject_responses(y)

	compare_bootstrap_fast_slow(InferenceContinUnivQuantileRegr$new(des, tau = 0.5, verbose = FALSE), SlowInferenceContinUnivQuantileRegr$new(des, tau = 0.5, verbose = FALSE), seed = 113)
	compare_bootstrap_fast_slow(InferenceContinMultiQuantileRegr$new(des, tau = 0.5, verbose = FALSE), SlowInferenceContinMultiQuantileRegr$new(des, tau = 0.5, verbose = FALSE), seed = 114)
})
