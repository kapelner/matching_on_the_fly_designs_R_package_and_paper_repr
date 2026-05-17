library(testthat)
library(EDI)

test_that("InferenceIncidProbitRegr matches stats::glm on the treatment slope", {
	set.seed(20260510)
	n <- 160
	x <- rnorm(n)
	z <- rnorm(n)
	w <- rep(c(1, 0), length.out = n)
	p <- pnorm(-0.25 + 0.9 * w + 0.45 * x - 0.35 * z)
	y <- rbinom(n, 1, p)

	des <- DesignFixed$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = x, z = z))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)

	inf <- InferenceIncidProbitRegr$new(des, verbose = FALSE, smart_cold_start_default = TRUE)
	fit_ref <- stats::glm(y ~ w + x + z, family = stats::binomial(link = "probit"))

	expect_equal(
		inf$compute_estimate(),
		as.numeric(stats::coef(fit_ref)[["w"]]),
		tolerance = 1e-4
	)
})

test_that("InferenceIncidProbitRegr supports asymptotic and likelihood-based inference paths", {
	set.seed(20260511)
	n <- 120
	x <- rnorm(n)
	w <- rep(c(1, 0), length.out = n)
	p <- pnorm(-0.1 + 0.75 * w + 0.3 * x)
	y <- rbinom(n, 1, p)

	des <- DesignFixed$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = x))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)

	inf <- InferenceIncidProbitRegr$new(des, verbose = FALSE, smart_cold_start_default = TRUE)
	est <- inf$compute_estimate()
	ci <- inf$compute_asymp_confidence_interval(alpha = 0.1)
	p_asymp <- inf$compute_asymp_two_sided_pval()

	expect_true(is.finite(est))
	expect_equal(length(ci), 2L)
	expect_true(all(is.finite(ci)))
	expect_true(is.finite(p_asymp))
	expect_gte(p_asymp, 0)
	expect_lte(p_asymp, 1)

	for (testing_type in c("score", "gradient", "lik_ratio")) {
		inf$set_testing_type(testing_type)
		pval <- inf$compute_asymp_two_sided_pval(0)
		expect_true(is.finite(pval))
		expect_gte(pval, 0)
		expect_lte(pval, 1)
	}
})
