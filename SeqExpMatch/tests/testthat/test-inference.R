test_that("Inference works for continuous", {
	n <- 10
	des <- SeqDesignCRD$new(n = n, response_type = "continuous", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	des$add_all_subject_responses(rnorm(n))

	# Simple Mean Diff
	inf <- SeqDesignInferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	expect_true(is.numeric(est))

	# OLS
	inf_ols <- SeqDesignInferenceContinMultOLS$new(des, verbose = FALSE)
	est_ols <- inf_ols$compute_treatment_estimate()
	expect_true(is.numeric(est_ols))
})

test_that("Inference works for incidence", {
	n <- 30
	des <- SeqDesignCRD$new(n = n, response_type = "incidence", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	des$add_all_subject_responses(rbinom(n, 1, 0.5))

	inf <- SeqDesignInferenceIncidUnivLogRegr$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	expect_true(is.numeric(est))
})

test_that("Inference works for count", {
	n <- 20
	des <- SeqDesignCRD$new(n = n, response_type = "count", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	des$add_all_subject_responses(rpois(n, 5))

	inf <- SeqDesignInferenceCountUnivNegBinRegr$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	expect_true(is.numeric(est))

	inf_robust_pois <- SeqDesignInferenceCountUnivRobustPoissonRegr$new(des, verbose = FALSE)
	est_robust_pois <- inf_robust_pois$compute_treatment_estimate()
	expect_true(is.numeric(est_robust_pois))

	inf_quasi_pois <- SeqDesignInferenceCountUnivQuasiPoissonRegr$new(des, verbose = FALSE)
	est_quasi_pois <- inf_quasi_pois$compute_treatment_estimate()
	expect_true(is.numeric(est_quasi_pois))

	if (requireNamespace("glmmTMB", quietly = TRUE)) {
		inf_zinb <- SeqDesignInferenceCountUnivZeroInflatedNegBinRegr$new(des, verbose = FALSE)
		est_zinb <- inf_zinb$compute_treatment_estimate()
		expect_true(is.numeric(est_zinb))
	}

	inf_hurdle_nb <- SeqDesignInferenceCountUnivHurdleNegBinRegr$new(des, verbose = FALSE)
	est_hurdle_nb <- inf_hurdle_nb$compute_treatment_estimate()
	expect_true(is.numeric(est_hurdle_nb))
})

test_that("Inference works for proportion", {
	n <- 10
	des <- SeqDesignCRD$new(n = n, response_type = "proportion", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	des$add_all_subject_responses(runif(n))

	inf <- SeqDesignInferencePropUniBetaRegr$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	expect_true(is.numeric(est))
})

test_that("Inference works for survival", {
	n <- 20
	des <- SeqDesignCRD$new(n = n, response_type = "survival")
	set.seed(1)
	for (i in 1:n) {
	des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	des$add_all_subject_responses(rexp(n), dead = rbinom(n, 1, 0.8))

	# KM Diff
	inf <- SeqDesignInferenceSurvivalKMDiff$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	expect_true(is.numeric(est))

	# Log-rank
	inf_logrank <- SeqDesignInferenceSurvivalLogRank$new(des, verbose = FALSE)
	est_logrank <- inf_logrank$compute_treatment_estimate()
	expect_true(is.numeric(est_logrank))

	# Cox PH
	inf_cox <- SeqDesignInferenceSurvivalUniCoxPHRegr$new(des, verbose = FALSE)
	est_cox <- inf_cox$compute_treatment_estimate()
	expect_true(is.numeric(est_cox))

	# Stratified Cox PH
	inf_strat_cox <- SeqDesignInferenceSurvivalUniStratCoxPHRegr$new(des, verbose = FALSE)
	est_strat_cox <- inf_strat_cox$compute_treatment_estimate()
	expect_true(is.numeric(est_strat_cox))
})

test_that("Inference works for ordinal partial proportional odds", {
	n <- 20
	des <- SeqDesignCRD$new(n = n, response_type = "ordinal")
	set.seed(10)
	for (i in 1:n) {
		des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	y_levels <- sample(1:4, n, replace = TRUE)
	des$add_all_subject_responses(y_levels)

	inf_ppod <- SeqDesignInferenceOrdinalPartialProportionalOdds$new(des, nonparallel = c("x"), verbose = FALSE)
	est_ppod <- inf_ppod$compute_treatment_estimate()
	expect_true(is.numeric(est_ppod))
	pval_ppod <- inf_ppod$compute_mle_two_sided_pval_for_treatment_effect()
	expect_true(is.numeric(pval_ppod))
})
