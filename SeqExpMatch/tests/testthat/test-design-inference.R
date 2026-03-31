test_that("Inference works for continuous", {
	n <- 10
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	des$add_all_subject_responses(rnorm(n))

	# Simple Mean Diff
	inf <- InferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	expect_true(is.numeric(est))

	# OLS
	inf_ols <- InferenceContinMultOLS$new(des, verbose = FALSE)
	est_ols <- inf_ols$compute_treatment_estimate()
	expect_true(is.numeric(est_ols))
})

test_that("Inference works for incidence", {
	n <- 30
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "incidence", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	des$add_all_subject_responses(rbinom(n, 1, 0.5))

	inf <- InferenceIncidUnivLogRegr$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	expect_true(is.numeric(est))
})

test_that("Simple incidence proportion difference uses pooled-variance t inference", {
	des <- FixedDesign$new(n = 10, response_type = "incidence", verbose = FALSE)
	for (i in 1:10) {
		des$add_subject(data.frame(x = i))
	}
	des$add_all_subject_assignments(c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0))
	des$add_all_subject_responses(c(1, 1, 0, 1, 0, 0, 1, 0, 0, 0))

	inf <- InferenceIncidenceSimplePropDiffPooled$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	ci <- inf$compute_asymp_confidence_interval(alpha = 0.05)
	pval <- inf$compute_asymp_two_sided_pval_for_treatment_effect()

	y_t <- c(1, 1, 0, 1, 0)
	y_c <- c(0, 1, 0, 0, 0)
	p_t <- mean(y_t)
	p_c <- mean(y_c)
	s2_t <- stats::var(y_t)
	s2_c <- stats::var(y_c)
	df <- length(y_t) + length(y_c) - 2
	s2_pooled <- ((length(y_t) - 1) * s2_t + (length(y_c) - 1) * s2_c) / df
	se <- sqrt(s2_pooled * (1 / length(y_t) + 1 / length(y_c)))
	crit <- stats::qt(0.975, df = df)
	expected_est <- p_t - p_c
	expected_ci <- c(expected_est - crit * se, expected_est + crit * se)
	expected_pval <- 2 * stats::pt(-abs(expected_est / se), df = df)

	expect_equal(est, expected_est, tolerance = 1e-12)
	expect_equal(ci, expected_ci, tolerance = 1e-12)
	expect_equal(pval, expected_pval, tolerance = 1e-12)

	des_bad <- DesignSeqOneByOneBernoulli$new(n = 6, response_type = "continuous", verbose = FALSE)
	for (i in 1:6) {
		des_bad$add_subject_to_experiment_and_assign(data.frame(x = i))
	}
	des_bad$add_all_subject_responses(rnorm(6))
	expect_error(InferenceIncidenceSimplePropDiffPooled$new(des_bad, verbose = FALSE), "incidence")
})

test_that("Azriel inference is gated to blocked incidence designs", {
	des <- FixedDesignBlocking$new(
		strata_cols = "stratum",
		n = 8,
		response_type = "incidence",
		verbose = FALSE
	)
	for (stratum in c("A", "A", "A", "A", "B", "B", "B", "B")) {
		des$add_subject(data.frame(stratum = stratum))
	}
	des$add_all_subject_assignments(c(1, 0, 1, 0, 1, 0, 1, 0))
	des$add_all_subject_responses(c(1, 0, 1, 0, 1, 1, 0, 0))

	inf <- InferenceIncidAzriel$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	ci <- inf$compute_asymp_confidence_interval()
	pval <- inf$compute_asymp_two_sided_pval_for_treatment_effect()

	y <- c(1, 0, 1, 0, 1, 1, 0, 0)
	w <- c(1, 0, 1, 0, 1, 0, 1, 0)
	m <- c(1, 1, 1, 1, 2, 2, 2, 2)
	expected_est <- mean(y[w == 1]) - mean(y[w == 0])
	var_cmh <- 0
	for (m_j in unique(m)) {
		i_m <- m == m_j
		y_m <- y[i_m]
		n_m <- sum(i_m)
		Sigma_m <- matrix(-1 / (n_m - 1), nrow = n_m, ncol = n_m)
		diag(Sigma_m) <- 1
		var_cmh <- var_cmh + as.numeric(t(y_m) %*% Sigma_m %*% y_m)
	}
	expected_se <- 2 / length(y) * sqrt(var_cmh)
	z_crit <- stats::qnorm(0.975)
	expected_ci <- c(expected_est - z_crit * expected_se, expected_est + z_crit * expected_se)
	expected_pval <- 2 * stats::pnorm(-abs(expected_est / expected_se))

	expect_equal(est, expected_est, tolerance = 1e-12)
	expect_equal(ci, expected_ci, tolerance = 1e-12)
	expect_equal(as.numeric(pval), expected_pval, tolerance = 1e-12)

	des_bad <- DesignSeqOneByOneBernoulli$new(n = 8, response_type = "incidence", verbose = FALSE)
	for (i in 1:8) {
		des_bad$add_subject_to_experiment_and_assign(data.frame(x = i))
	}
	des_bad$add_all_subject_responses(rbinom(8, 1, 0.5))
	expect_error(InferenceIncidAzriel$new(des_bad, verbose = FALSE), "blocking design")
})

test_that("Inference works for count", {
	n <- 20
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	des$add_all_subject_responses(rpois(n, 5))

	inf <- InferenceCountUnivNegBinRegr$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	expect_true(is.numeric(est))

	inf_robust_pois <- InferenceCountUnivRobustPoissonRegr$new(des, verbose = FALSE)
	est_robust_pois <- inf_robust_pois$compute_treatment_estimate()
	expect_true(is.numeric(est_robust_pois))

	inf_quasi_pois <- InferenceCountUnivQuasiPoissonRegr$new(des, verbose = FALSE)
	est_quasi_pois <- inf_quasi_pois$compute_treatment_estimate()
	expect_true(is.numeric(est_quasi_pois))

	if (requireNamespace("glmmTMB", quietly = TRUE)) {
		inf_zinb <- InferenceCountUnivZeroInflatedNegBinRegr$new(des, verbose = FALSE)
		est_zinb <- inf_zinb$compute_treatment_estimate()
		expect_true(is.numeric(est_zinb))
	}

	inf_hurdle_nb <- InferenceCountUnivHurdleNegBinRegr$new(des, verbose = FALSE)
	est_hurdle_nb <- inf_hurdle_nb$compute_treatment_estimate()
	expect_true(is.numeric(est_hurdle_nb))
})

test_that("Inference works for proportion", {
	n <- 10
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "proportion", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	des$add_all_subject_responses(runif(n))

	inf <- InferencePropUniBetaRegr$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	expect_true(is.numeric(est))
})

test_that("Inference works for survival", {
	n <- 20
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "survival")
	set.seed(1)
	for (i in 1:n) {
	des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	des$add_all_subject_responses(rexp(n), dead = rbinom(n, 1, 0.8))

	# KM Diff
	inf <- InferenceSurvivalKMDiff$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	expect_true(is.numeric(est))

	# Log-rank
	inf_logrank <- InferenceSurvivalLogRank$new(des, verbose = FALSE)
	est_logrank <- inf_logrank$compute_treatment_estimate()
	expect_true(is.numeric(est_logrank))

	# Cox PH
	inf_cox <- InferenceSurvivalUniCoxPHRegr$new(des, verbose = FALSE)
	est_cox <- inf_cox$compute_treatment_estimate()
	expect_true(is.numeric(est_cox))

	# Stratified Cox PH
	inf_strat_cox <- InferenceSurvivalUniStratCoxPHRegr$new(des, verbose = FALSE)
	est_strat_cox <- inf_strat_cox$compute_treatment_estimate()
	expect_true(is.numeric(est_strat_cox))
})

test_that("Inference works for ordinal partial proportional odds", {
	n <- 20
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "ordinal")
	set.seed(10)
	for (i in 1:n) {
		des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	y_levels <- sample(1:4, n, replace = TRUE)
	des$add_all_subject_responses(y_levels)

	inf_ppod <- InferenceOrdinalPartialProportionalOdds$new(des, nonparallel = c("x"), verbose = FALSE)
	est_ppod <- inf_ppod$compute_treatment_estimate()
	expect_true(is.numeric(est_ppod))
	pval_ppod <- inf_ppod$compute_asymp_two_sided_pval_for_treatment_effect()
	expect_true(is.numeric(pval_ppod))
})

test_that("Inference works for incidence KK Newcombe IVWC", {
	n <- 20
	set.seed(1)
	x_dat <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1, 0.5))
	y <- rbinom(n, 1, 0.5)

	seq_des <- DesignSeqOneByOneKK14$new(n = n, response_type = "incidence", verbose = FALSE)
	for (i in 1:n) {
		seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
	}
	seq_des$add_all_subject_responses(y)

	inf <- InferenceIncidUnivKKNewcombeRiskDiff$new(seq_des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	expect_true(is.numeric(est))
	
	ci <- inf$compute_asymp_confidence_interval()
	expect_true(is.numeric(ci))
	expect_length(ci, 2)
	
	pval <- inf$compute_asymp_two_sided_pval_for_treatment_effect()
	expect_true(is.numeric(pval))
})
