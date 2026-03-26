.libPaths(c("Rlib", .libPaths()))
library(testthat)
library(SeqExpMatch)
library(data.table)

test_that("compute_confidence_interval_rand works for continuous response", {
	set.seed(123)
	n <- 40
	des <- SeqDesignBernoulli$new(n = n, response_type = "continuous", verbose = TRUE)
	for (i in 1:n) {
	des$add_subject_to_experiment_and_assign(data.table(x1 = rnorm(1), x2 = rnorm(1)))
	}
	# Add treatment effect of 1.0
	treatment <- des$.__enclos_env__$private$w
	y <- rnorm(n) + treatment * 1.0
	des$add_all_subject_responses(y)

	inf <- DesignInferenceAllSimpleMeanDiff$new(des, verbose = TRUE)

	# Compute randomization CI
	# Using small nsim for speed in tests
	ci <- inf$compute_confidence_interval_rand(alpha = 0.05, r = 100, pval_epsilon = 0.05)

	expect_equal(length(ci), 2)
	expect_true(ci[1] < ci[2])
	# The estimate should be within the CI
	est <- inf$compute_treatment_estimate()
	expect_true(est >= ci[1] && est <= ci[2])

	message("Continuous Rand CI: [", ci[1], ", ", ci[2], "] Est: ", est)
})

test_that("compute_confidence_interval_rand works for proportion response", {
	set.seed(123)
	n <- 100
	des <- SeqDesignBernoulli$new(n = n, response_type = "proportion", verbose = FALSE)
	for (i in 1:n) {
	des$add_subject_to_experiment_and_assign(data.table(x1 = rnorm(1)))
	}
	treatment <- des$.__enclos_env__$private$w
	# Simulate proportions
	mu <- plogis(-0.5 + treatment * 1.0)
	y <- rbeta(n, shape1 = mu * 10, shape2 = (1 - mu) * 10)
	des$add_all_subject_responses(y)

	inf <- DesignInferencePropUniBetaRegr$new(des, verbose = FALSE)

	# Compute randomization CI
	ci <- inf$compute_confidence_interval_rand(alpha = 0.05, r = 100, pval_epsilon = 0.05)

	expect_equal(length(ci), 2)
	expect_true(ci[1] < ci[2])
	expect_true(all(is.finite(ci)))

	est <- inf$compute_treatment_estimate()
	expect_true(est >= ci[1] && est <= ci[2])

	message("Proportion Rand CI: [", ci[1], ", ", ci[2], "] Est: ", est)
})

test_that("compute_confidence_interval_rand works for survival response (uncensored)", {
	set.seed(123)
	n <- 50
	des <- SeqDesignBernoulli$new(n = n, response_type = "survival", verbose = FALSE)
	for (i in 1:n) {
	des$add_subject_to_experiment_and_assign(data.table(x1 = rnorm(1)))
	}
	treatment <- des$.__enclos_env__$private$w
	# Simulate survival times (log-normal)
	y <- exp(1.0 + treatment * 0.8 + rnorm(n, 0, 0.5))
	des$add_all_subject_responses(y, deads = rep(1, n)) # All events, no censoring

	inf <- DesignInferenceSurvivalUniWeibullRegr$new(des, verbose = FALSE)

	# Compute randomization CI
	ci <- inf$compute_confidence_interval_rand(alpha = 0.05, r = 100, pval_epsilon = 0.05)

	expect_equal(length(ci), 2)
	expect_true(ci[1] < ci[2])
	# Survival times (ratios) should be positive
	expect_true(all(ci > 0))

	est <- inf$compute_treatment_estimate()
	expect_true(est >= ci[1] && est <= ci[2])

	message("Survival Rand CI: [", ci[1], ", ", ci[2], "] Est: ", est)
})

test_that("compute_confidence_interval_rand works for ordinal response (cumulative logit)", {
	set.seed(456)
	n <- 40
	des <- SeqDesignBernoulli$new(n = n, response_type = "ordinal", verbose = FALSE)
	for (i in 1:n) {
		des$add_subject_to_experiment_and_assign(data.table(x1 = rnorm(1)))
	}
	treatment <- des$.__enclos_env__$private$w
	score <- rnorm(n, mean = treatment * 0.5)
	y <- 1L + (score > 0) + (score > 1)
	des$add_all_subject_responses(y)

	inf <- DesignInferenceOrdinalUniPropOddsRegr$new(des, verbose = FALSE)

	ci <- inf$compute_confidence_interval_rand(alpha = 0.05, r = 100, pval_epsilon = 0.05, show_progress = FALSE)

	expect_equal(length(ci), 2)
	expect_true(ci[1] < ci[2])
	est <- inf$compute_treatment_estimate()
	expect_true(est >= ci[1] && est <= ci[2])
	expect_true(all(is.finite(ci)))

	message("Ordinal Rand CI: [", ci[1], ", ", ci[2], "] Est: ", est)
})

test_that("compute_confidence_interval_rand throws error for unsupported types", {
	n <- 20
	des_incid <- SeqDesignBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	for (i in 1:n) des_incid$add_subject_to_experiment_and_assign(data.table(x=1))
	des_incid$add_all_subject_responses(rbinom(n, 1, 0.5))
	inf_incid <- DesignInferenceIncidUnivLogRegr$new(des_incid)
	expect_error(inf_incid$compute_confidence_interval_rand(), "Confidence intervals are not supported for randomization tests for incidence outcomes")

	des_count <- SeqDesignBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	for (i in 1:n) des_count$add_subject_to_experiment_and_assign(data.table(x=1))
	des_count$add_all_subject_responses(rpois(n, 5))
	inf_count <- DesignInferenceCountUnivNegBinRegr$new(des_count)
	ci_count <- inf_count$compute_confidence_interval_rand(alpha = 0.05, r = 100, pval_epsilon = 0.05)
	expect_equal(length(ci_count), 2)
	expect_true(ci_count[1] < ci_count[2])
	expect_true(all(is.finite(ci_count)))
})

test_that("Zhang exact incidence CI matches across serial and multicore", {
	set.seed(321)
	n <- 24
	des <- SeqDesignBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	for (i in 1:n) {
		des$add_subject_to_experiment_and_assign(data.table(x1 = rnorm(1), x2 = rnorm(1)))
	}
	treatment <- des$.__enclos_env__$private$w
	prob <- plogis(-0.2 + 0.8 * treatment)
	des$add_all_subject_responses(rbinom(n, 1, prob))

	inf_serial <- DesignInferenceIncidExactZhang$new(des, num_cores = 1, verbose = FALSE)
	inf_parallel <- DesignInferenceIncidExactZhang$new(des, num_cores = 4, verbose = FALSE)

	ci_serial <- inf_serial$compute_exact_confidence_interval(alpha = 0.10, pval_epsilon = 0.01)
	ci_parallel <- inf_parallel$compute_exact_confidence_interval(alpha = 0.10, pval_epsilon = 0.01)

	expect_length(ci_serial, 2)
	expect_true(all(is.finite(ci_serial)))
	expect_equal(ci_parallel, ci_serial, tolerance = 1e-8)
})
