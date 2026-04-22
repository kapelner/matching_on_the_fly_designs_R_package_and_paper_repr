test_that("Weibull Frailty Inference works for KK designs", {
	skip_if_not_installed("parfm")
	
	n <- 40
	# Use KK14 design to ensure we have matches
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = "survival", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
	}
	# Generate some survival data with many events to ensure convergence
	y <- rexp(n, rate = 0.1)
	dead <- rep(1L, n)
	add_all_subject_responses_seq(des, y, deads = dead)

	# 1. Univariate IVWC (model_formula = ~ 1)
	inf_univ_ivwc <- InferenceSurvivalKKWeibullFrailtyIVWC$new(des, model_formula = ~ 1, verbose = FALSE)
	est_univ_ivwc <- inf_univ_ivwc$compute_estimate()
	expect_true(is.numeric(est_univ_ivwc))
	expect_true(is.finite(est_univ_ivwc))
	
	ci_univ_ivwc <- inf_univ_ivwc$compute_asymp_confidence_interval()
	expect_length(ci_univ_ivwc, 2)
	expect_true(all(is.finite(ci_univ_ivwc)))

	# 2. Multivariate IVWC (default)
	inf_multi_ivwc <- InferenceSurvivalKKWeibullFrailtyIVWC$new(des, verbose = FALSE)
	est_multi_ivwc <- inf_multi_ivwc$compute_estimate()
	expect_true(is.numeric(est_multi_ivwc))
	expect_true(is.finite(est_multi_ivwc))

	# 3. Univariate OneLik (model_formula = ~ 1)
	inf_univ_onelik <- InferenceSurvivalKKWeibullFrailtyOneLik$new(des, model_formula = ~ 1, verbose = FALSE)
	est_univ_onelik <- inf_univ_onelik$compute_estimate()
	expect_true(is.numeric(est_univ_onelik))
	expect_true(is.finite(est_univ_onelik))
	
	pv_univ_onelik <- inf_univ_onelik$compute_asymp_two_sided_pval()
	expect_true(is.numeric(pv_univ_onelik))
	expect_true(pv_univ_onelik >= 0 && pv_univ_onelik <= 1)

	# 4. Multivariate OneLik (default)
	inf_multi_onelik <- InferenceSurvivalKKWeibullFrailtyOneLik$new(des, verbose = FALSE)
	est_multi_onelik <- inf_multi_onelik$compute_estimate()
	expect_true(is.numeric(est_multi_onelik))
	expect_true(is.finite(est_multi_onelik))
})
