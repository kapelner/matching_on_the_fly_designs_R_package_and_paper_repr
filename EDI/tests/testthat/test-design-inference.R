test_that("Inference works for continuous", {
	n <- 10
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	add_all_subject_responses_seq(des, rnorm(n))

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
	des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	add_all_subject_responses_seq(des, rbinom(n, 1, 0.5))

	inf <- InferenceIncidUnivLogRegr$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	expect_true(is.numeric(est))
})

test_that("Simple incidence proportion difference uses pooled-variance t inference", {
	des <- FixedDesign$new(n = 10, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = 1:10))
	des$overwrite_all_subject_assignments(c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0))
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
	expect_equal(unname(ci), expected_ci, tolerance = 1e-12)
	expect_equal(pval, expected_pval, tolerance = 1e-12)

	des_bad <- DesignSeqOneByOneBernoulli$new(n = 6, response_type = "continuous", verbose = FALSE)
	for (i in 1:6) {
		des_bad$add_one_subject_to_experiment_and_assign(data.frame(x = i))
	}
	add_all_subject_responses_seq(des_bad, rnorm(6))
	expect_error(InferenceIncidenceSimplePropDiffPooled$new(des_bad, verbose = FALSE), "incidence")
})

test_that("Azriel inference is gated to blocked incidence designs", {
	des <- FixedDesignBlocking$new(
		strata_cols = "stratum",
		n = 8,
		response_type = "incidence",
		verbose = FALSE
	)
	des$add_all_subjects_to_experiment(data.frame(stratum = c("A", "A", "A", "A", "B", "B", "B", "B")))
	des$overwrite_all_subject_assignments(c(1, 0, 1, 0, 1, 0, 1, 0))
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
	expected_ci <- c(-0.4989475801457, 1.4989475801457)
	expected_pval <- 0.2665697

	expect_equal(est, expected_est, tolerance = 1e-12)
	expect_equal(unname(ci), expected_ci, tolerance = 1e-12)
	expect_equal(as.numeric(pval), expected_pval, tolerance = 1e-8)

	des_bad <- DesignSeqOneByOneBernoulli$new(n = 8, response_type = "incidence", verbose = FALSE)
	for (i in 1:8) {
		des_bad$add_one_subject_to_experiment_and_assign(data.frame(x = i))
	}
	add_all_subject_responses_seq(des_bad, rbinom(8, 1, 0.5))
	expect_error(InferenceIncidAzriel$new(des_bad, verbose = FALSE), "blocking design")
})

test_that("Azriel inference requires even treatment allocation", {
	des <- FixedDesignBlocking$new(
		strata_cols = "stratum",
		n = 8,
		response_type = "incidence",
		prob_T = 0.25,
		verbose = FALSE
	)
	des$add_all_subjects_to_experiment(data.frame(stratum = c("A", "A", "A", "A", "B", "B", "B", "B")))
	des$overwrite_all_subject_assignments(c(1, 0, 1, 0, 1, 0, 1, 0))
	des$add_all_subject_responses(c(1, 0, 1, 0, 1, 1, 0, 0))

	expect_error(
		InferenceIncidAzriel$new(des, verbose = FALSE),
		"even treatment allocation"
	)
})

test_that("Azriel inference requires equal block sizes", {
	des <- FixedDesignBlocking$new(
		strata_cols = "stratum",
		n = 8,
		response_type = "incidence",
		verbose = FALSE
	)
	des$add_all_subjects_to_experiment(data.frame(stratum = c(rep("A", 3), rep("B", 5))))
	des$overwrite_all_subject_assignments(c(1, 0, 1, 0, 1, 0, 1, 0))
	des$add_all_subject_responses(c(1, 0, 1, 0, 1, 1, 0, 0))

	expect_error(
		InferenceIncidAzriel$new(des, verbose = FALSE),
		"same number of subjects"
	)
})

test_that("Azriel and Extended Robins standard errors match a fixed blocked simulation", {
	set.seed(20260405)
	expected <- matrix(
		c(
			0.3818813079129866, 0.6059599821770412,
			0.3535533905932738, 0.4330127018922193,
			0.2886751345948129, 0.5448623679425842,
			0.2500000000000000, 0.3423265984407288,
			0.3818813079129866, 0.6059599821770412
		),
		ncol = 2L,
		byrow = TRUE,
		dimnames = list(
			NULL,
			c("azriel", "robins")
		)
	)

	ys <- replicate(5, rbinom(8, 1, 0.5), simplify = FALSE)
	se_pairs <- t(vapply(ys, function(y) {
		des <- FixedDesignBlocking$new(
			strata_cols = "stratum",
			n = 8,
			response_type = "incidence",
			verbose = FALSE
		)
		des$add_all_subjects_to_experiment(data.frame(stratum = c(rep("A", 4), rep("B", 4))))
		des$overwrite_all_subject_assignments(c(1, 0, 1, 0, 1, 0, 1, 0))
		des$add_all_subject_responses(y)

		inf_azriel <- InferenceIncidAzriel$new(des, verbose = FALSE)
		inf_robins <- InferenceIncidExtendedRobins$new(des, verbose = FALSE)
		c(
			azriel = inf_azriel$.__enclos_env__$private$get_standard_error(),
			robins = inf_robins$.__enclos_env__$private$get_standard_error()
		)
	}, numeric(2)))

	expect_equal(se_pairs, expected, tolerance = 1e-12)
})

test_that("Extended Robins standard error matches the blockwise formula", {
	des <- FixedDesignBlocking$new(
		strata_cols = "stratum",
		n = 8,
		response_type = "incidence",
		verbose = FALSE
	)
	des$add_all_subjects_to_experiment(data.frame(stratum = c(rep("A", 4), rep("B", 4))))
	des$overwrite_all_subject_assignments(c(1, 0, 1, 0, 1, 0, 1, 0))
	des$add_all_subject_responses(c(1, 0, 1, 0, 1, 1, 0, 0))

	inf <- InferenceIncidExtendedRobins$new(des, verbose = FALSE)
	se_cpp <- inf$.__enclos_env__$private$get_standard_error()

	m <- des$get_block_ids()
	B <- length(unique(m))
	n_B <- sum(m == 1L)
	n_B_over_two <- n_B / 2
	variance_tot <- 0
	for (b in 1:B) {
		y_b <- des$get_y()[m == b]
		w_b <- des$get_w()[m == b]
		p_hat_T_b <- sum(y_b[w_b == 1]) / n_B_over_two
		p_hat_C_b <- sum(y_b[w_b == 0]) / n_B_over_two
		m_1_b <- max(p_hat_T_b, p_hat_C_b)
		m_0_b <- min(p_hat_T_b, p_hat_C_b)
		variance_tot <- variance_tot +
			m_1_b * (1 - m_1_b) / n_B_over_two +
			m_0_b * (1 - m_0_b) / n_B_over_two +
			((2 * m_0_b - m_1_b) * (1 - m_1_b) - m_0_b * (1 - m_0_b)) / n_B
	}
	p_hat_T <- mean(des$get_y()[des$get_w() == 1])
	p_hat_C <- mean(des$get_y()[des$get_w() == 0])
	var_robbins_ext <- 1 / des$get_n() * (
		p_hat_T * (1 - p_hat_T) + p_hat_C * (1 - p_hat_C)
	)
	se_r <- sqrt(variance_tot + var_robbins_ext)

	expect_equal(se_cpp, se_r, tolerance = 1e-12)
})

test_that("G-computation risk-ratio intervals error when log-scale bounds overflow", {
	des <- FixedDesigniBCRD$new(n = 4, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = 1:4))
	des$overwrite_all_subject_assignments(c(1, 0, 1, 0))
	des$add_all_subject_responses(c(1, 0, 1, 0))

	inf <- InferenceIncidMultiGCompRiskRatio$new(des, verbose = FALSE)
	priv <- inf$.__enclos_env__$private
	priv$cached_values$rr <- 1
	priv$cached_values$s_beta_hat_T <- 1
	priv$cached_values$log_rr <- 1000
	priv$cached_values$se_log_rr <- 100

	expect_error(
		inf$compute_asymp_confidence_interval(),
		"could not compute a finite delta-method confidence interval"
	)
})

test_that("Inference works for count", {
	n <- 20
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "count", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	add_all_subject_responses_seq(des, rpois(n, 5))

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

test_that("KK count combined-likelihood multi inference handles full-width covariates", {
	skip_if_not_installed("MASS")
	skip_if_not_installed("dplyr")

	set.seed(1)
	cars_subset <- MASS::Cars93 |>
		stats::na.omit() |>
		dplyr::slice_sample(n = 48, replace = TRUE) |>
		dplyr::select(-Make, -Model)
	X <- model.matrix(Price ~ 0 + ., data = cars_subset)
	X <- apply(X, 2, scale)
	X <- X / ncol(X)
	X <- as.data.frame(X)
	X <- dplyr::select(X, dplyr::where(~ !any(is.na(.))))
	y <- round(cars_subset$Price - min(cars_subset$Price))

	des <- DesignSeqOneByOneKK14$new(n = nrow(X), response_type = "count", verbose = FALSE)
	for (i in seq_len(nrow(X))) {
		w_i <- des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
		y_i <- round(y[i] * exp(rnorm(1, mean = 0, sd = 0.1) * w_i))
		des$add_one_subject_response(i, y_i, 1)
	}

	inf <- InferenceCountPoissonMultiKKCPoissonCombinedLikelihood$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()

	expect_true(is.numeric(est))
	expect_length(est, 1L)
	expect_true(is.finite(est))
})

test_that("Inference works for proportion", {
	n <- 10
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "proportion", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
	des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	add_all_subject_responses_seq(des, runif(n))

	inf <- InferencePropUniBetaRegr$new(des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	expect_true(is.numeric(est))
})

test_that("Inference works for survival", {
	n <- 20
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "survival")
	set.seed(1)
	for (i in 1:n) {
	des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	add_all_subject_responses_seq(des, rexp(n), deads = rbinom(n, 1, 0.8))

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
		des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	y_levels <- sample(1:4, n, replace = TRUE)
	add_all_subject_responses_seq(des, y_levels)

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
		seq_des$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
	}
	add_all_subject_responses_seq(seq_des, y)

	inf <- InferenceIncidUnivKKNewcombeRiskDiff$new(seq_des, verbose = FALSE)
	est <- inf$compute_treatment_estimate()
	expect_true(is.numeric(est))
	
	ci <- inf$compute_asymp_confidence_interval()
	expect_true(is.numeric(ci))
	expect_length(ci, 2)
	
	pval <- inf$compute_asymp_two_sided_pval_for_treatment_effect()
	expect_true(is.numeric(pval))
})
