library(testthat)
library(EDI)

test_that("KK OLS one-lik uses HC2 post-fit standard errors and finite df", {
	des = FixedDesignBinaryMatch$new(response_type = "continuous", n = 6, verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = c(-2, -1, 0, 1, 2, 3)))
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(c(-1, 0, 1, 1, 4, 9))

	inf = InferenceContinKKOLSOneLik$new(des, verbose = FALSE)
	priv = inf$.__enclos_env__$private

	X = rbind(
		cbind(0, 1, c(-2, 0, 2)),
		cbind(1, c(0, 1, 0), c(-1, 1, 3))
	)
	y = c(1, 1, 5, 0, 1, 9)
	j_treat = 2L

	fit = priv$fit_ols(X, y, j_treat = j_treat, estimate_only = FALSE)
	post_fit = EDI:::ols_hc2_post_fit_cpp(X, y, fit$b, j_treat)

	expect_equal(fit$se, as.numeric(post_fit$std_err[j_treat]), tolerance = 1e-10)
	expect_equal(fit$vcov, post_fit$vcov, tolerance = 1e-10)
	expect_equal(fit$df, nrow(X) - ncol(X))

	inf$compute_estimate()
	expect_true(is.finite(priv$cached_values$df))
	expect_gt(priv$cached_values$df, 0)
})

test_that("KK OLS IVWC uses HC2 component standard errors and finite df", {
	des = FixedDesignBinaryMatch$new(response_type = "continuous", n = 6, verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = c(-2, -1, 0, 1, 2, 3)))
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(c(-1, 0, 1, 1, 4, 9))

	inf = InferenceContinKKOLSIVWC$new(des, verbose = FALSE)
	priv = inf$.__enclos_env__$private

	X = cbind(1, c(0, 1, 0, 1, 0, 1), c(-2, -1, 0, 1, 2, 3))
	y = c(-1, 0, 1, 1, 4, 9)
	j_treat = 2L

	fit = priv$fit_ols_with_treatment(X, y, j_treat = j_treat, estimate_only = FALSE)
	post_fit = EDI:::ols_hc2_post_fit_cpp(X, y, as.numeric(stats::coef(stats::lm.fit(X, y))), j_treat)

	expect_equal(sqrt(fit$ssq), as.numeric(post_fit$std_err[j_treat]), tolerance = 1e-10)
	expect_equal(fit$df, nrow(X) - ncol(X))
})
