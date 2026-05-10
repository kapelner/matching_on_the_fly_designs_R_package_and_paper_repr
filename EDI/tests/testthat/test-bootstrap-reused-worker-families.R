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
	SlowInferenceIncidLogBinomial = R6::R6Class(
		"SlowInferenceIncidLogBinomial",
		inherit = InferenceIncidLogBinomial,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidLogBinomial = R6::R6Class(
		"SlowInferenceIncidLogBinomial",
		inherit = InferenceIncidLogBinomial,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260410)
	n = 60
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
	des = DesignFixedBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	p = pmin(0.95, pmax(0.02, exp(-1.4 + 0.5 * w + 0.2 * X$x1 - 0.1 * X$x2)))
	y = stats::rbinom(n, 1, p)
	des$add_all_subject_responses(y)

	compare_bootstrap_fast_slow(
		InferenceIncidLogBinomial$new(des, verbose = FALSE),
		SlowInferenceIncidLogBinomial$new(des, verbose = FALSE),
		seed = 101
	)
	compare_bootstrap_fast_slow(
		InferenceIncidLogBinomial$new(des, verbose = FALSE),
		SlowInferenceIncidLogBinomial$new(des, verbose = FALSE),
		seed = 102
	)
})

test_that("binomial identity reusable bootstrap worker matches generic path", {
	SlowInferenceIncidBinomialIdentityRiskDiff = R6::R6Class(
		"SlowInferenceIncidBinomialIdentityRiskDiff",
		inherit = InferenceIncidBinomialIdentityRiskDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceIncidBinomialIdentityRiskDiff = R6::R6Class(
		"SlowInferenceIncidBinomialIdentityRiskDiff",
		inherit = InferenceIncidBinomialIdentityRiskDiff,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260411)
	n = 64
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
	des = DesignFixedBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	p = pmin(0.95, pmax(0.02, 0.25 + 0.18 * w + 0.05 * X$x1 - 0.03 * X$x2))
	y = stats::rbinom(n, 1, p)
	des$add_all_subject_responses(y)

	compare_bootstrap_fast_slow(
		InferenceIncidBinomialIdentityRiskDiff$new(des, verbose = FALSE),
		SlowInferenceIncidBinomialIdentityRiskDiff$new(des, verbose = FALSE),
		seed = 103
	)
	compare_bootstrap_fast_slow(
		InferenceIncidBinomialIdentityRiskDiff$new(des, verbose = FALSE),
		SlowInferenceIncidBinomialIdentityRiskDiff$new(des, verbose = FALSE),
		seed = 104
	)
})

test_that("count reusable bootstrap workers match generic paths", {
	SlowInferenceCountPoisson = R6::R6Class(
		"SlowInferenceCountPoisson",
		inherit = InferenceCountPoisson,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountPoisson = R6::R6Class(
		"SlowInferenceCountPoisson",
		inherit = InferenceCountPoisson,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountRobustPoisson = R6::R6Class(
		"SlowInferenceCountRobustPoisson",
		inherit = InferenceCountRobustPoisson,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountRobustPoisson = R6::R6Class(
		"SlowInferenceCountRobustPoisson",
		inherit = InferenceCountRobustPoisson,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountQuasiPoisson = R6::R6Class(
		"SlowInferenceCountQuasiPoisson",
		inherit = InferenceCountQuasiPoisson,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceCountQuasiPoisson = R6::R6Class(
		"SlowInferenceCountQuasiPoisson",
		inherit = InferenceCountQuasiPoisson,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260412)
	n = 68
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
	des = DesignFixedBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	lambda = exp(0.2 + 0.35 * w + 0.15 * X$x1 - 0.1 * X$x2)
	y = stats::rpois(n, lambda)
	des$add_all_subject_responses(y)

	compare_bootstrap_fast_slow(InferenceCountPoisson$new(des, verbose = FALSE), SlowInferenceCountPoisson$new(des, verbose = FALSE), seed = 105)
	compare_bootstrap_fast_slow(InferenceCountPoisson$new(des, verbose = FALSE), SlowInferenceCountPoisson$new(des, verbose = FALSE), seed = 106)
	compare_bootstrap_fast_slow(InferenceCountRobustPoisson$new(des, verbose = FALSE), SlowInferenceCountRobustPoisson$new(des, verbose = FALSE), seed = 107)
	compare_bootstrap_fast_slow(InferenceCountRobustPoisson$new(des, verbose = FALSE), SlowInferenceCountRobustPoisson$new(des, verbose = FALSE), seed = 108)
	compare_bootstrap_fast_slow(InferenceCountQuasiPoisson$new(des, verbose = FALSE), SlowInferenceCountQuasiPoisson$new(des, verbose = FALSE), seed = 109)
	compare_bootstrap_fast_slow(InferenceCountQuasiPoisson$new(des, verbose = FALSE), SlowInferenceCountQuasiPoisson$new(des, verbose = FALSE), seed = 110)
})

test_that("count likelihood phase 1 and 2 classes share the count likelihood base", {
	set.seed(20260415)
	n = 40
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = DesignFixedBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	y = stats::rpois(n, exp(0.3 + 0.25 * w + 0.1 * X$x1))
	des$add_all_subject_responses(y)

	expect_true(inherits(InferenceCountPoisson$new(des, verbose = FALSE), "InferenceCountLikelihood"))
	expect_true(inherits(InferenceCountNegBin$new(des, verbose = FALSE), "InferenceCountLikelihood"))
	expect_true(inherits(InferenceCountZeroInflatedPoisson$new(des, verbose = FALSE), "InferenceCountLikelihood"))
	expect_true(inherits(InferenceCountZeroInflatedNegBin$new(des, verbose = FALSE), "InferenceCountLikelihood"))
	expect_true(inherits(InferenceCountHurdlePoisson$new(des, verbose = FALSE), "InferenceCountLikelihood"))
	expect_true(inherits(InferenceCountHurdleNegBin$new(des, verbose = FALSE), "InferenceCountLikelihood"))
	expect_true(inherits(InferenceCountRobustPoisson$new(des, verbose = FALSE), "InferenceCountCompositeLikelihood"))
	expect_true(inherits(InferenceCountQuasiPoisson$new(des, verbose = FALSE), "InferenceCountCompositeLikelihood"))
})

test_that("continuous robust reusable bootstrap worker matches generic path", {
	SlowInferenceContinRobustRegr = R6::R6Class(
		"SlowInferenceContinRobustRegr",
		inherit = InferenceContinRobustRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceContinRobustRegr = R6::R6Class(
		"SlowInferenceContinRobustRegr",
		inherit = InferenceContinRobustRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260413)
	n = 60
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
	des = DesignFixedBernoulli$new(n = n, response_type = "continuous", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	y = 0.4 * w + 0.3 * X$x1 - 0.2 * X$x2 + stats::rt(n, df = 5)
	des$add_all_subject_responses(y)

	compare_bootstrap_fast_slow(InferenceContinRobustRegr$new(des, method = "M", verbose = FALSE), SlowInferenceContinRobustRegr$new(des, method = "M", verbose = FALSE), seed = 111)
	compare_bootstrap_fast_slow(InferenceContinRobustRegr$new(des, method = "M", verbose = FALSE), SlowInferenceContinRobustRegr$new(des, method = "M", verbose = FALSE), seed = 112)
})

test_that("continuous quantile reusable bootstrap worker matches generic path", {
	skip_if_not_installed("quantreg")

	SlowInferenceContinQuantileRegr = R6::R6Class(
		"SlowInferenceContinQuantileRegr",
		inherit = InferenceContinQuantileRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)
	SlowInferenceContinQuantileRegr = R6::R6Class(
		"SlowInferenceContinQuantileRegr",
		inherit = InferenceContinQuantileRegr,
		private = list(supports_reusable_bootstrap_worker = function() FALSE)
	)

	set.seed(20260414)
	n = 56
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
	des = DesignFixedBernoulli$new(n = n, response_type = "continuous", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()
	w = des$get_w()
	y = 0.5 * w + 0.25 * X$x1 - 0.15 * X$x2 + stats::rt(n, df = 6)
	des$add_all_subject_responses(y)

	compare_bootstrap_fast_slow(InferenceContinQuantileRegr$new(des, tau = 0.5, verbose = FALSE), SlowInferenceContinQuantileRegr$new(des, tau = 0.5, verbose = FALSE), seed = 113)
	compare_bootstrap_fast_slow(InferenceContinQuantileRegr$new(des, tau = 0.5, verbose = FALSE), SlowInferenceContinQuantileRegr$new(des, tau = 0.5, verbose = FALSE), seed = 114)
})
