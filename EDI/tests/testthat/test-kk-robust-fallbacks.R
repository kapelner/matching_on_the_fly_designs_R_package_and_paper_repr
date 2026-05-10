library(testthat)
library(EDI)

test_that("KK robust regression muffles MM non-convergence and falls back internally", {
	skip_if_not_installed("MASS")

	des = DesignFixedBinaryMatch$new(response_type = "continuous", n = 4, verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = c(-1, -0.5, 0.5, 1)))
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(c(0, 0, 100, 100))

	inf_one = InferenceContinKKRobustRegrOneLik$new(des, method = "MM", maxit = 1, acc = 1e-12, verbose = FALSE)
	inf_ivwc = InferenceContinKKRobustRegrIVWC$new(des, method = "MM", maxit = 1, acc = 1e-12, verbose = FALSE)

	X = cbind(1, stats::rnorm(20))
	y = c(rep(0, 10), rep(100, 10))

	fit_one = expect_no_warning(inf_one$.__enclos_env__$private$fit_rlm(X, y, j_treat = 2L, estimate_only = FALSE))
	fit_ivwc = expect_no_warning(inf_ivwc$.__enclos_env__$private$fit_rlm_with_treatment(X, y, j_treat = 2L, estimate_only = FALSE))

	expect_true(is.list(fit_one))
	expect_true(is.finite(fit_one$beta))
	expect_true(is.finite(fit_one$se))
	expect_true(is.list(fit_ivwc))
	expect_true(is.finite(fit_ivwc$beta))
	expect_true(is.finite(fit_ivwc$ssq))
})

test_that("KK hurdle-Poisson combined-likelihood warns once before bootstrap fallback", {
	des = DesignFixedBinaryMatch$new(response_type = "count", n = 4, verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = c(-1, -0.5, 0.5, 1)))
	des$assign_w_to_all_subjects()
	des$add_all_subject_responses(c(0L, 1L, 2L, 0L))

	inf = InferenceCountKKHurdlePoissonOneLik$new(des, verbose = FALSE, use_rcpp = TRUE)
	inf$.__enclos_env__$private$cached_values$warned_bootstrap_se_unavailable = FALSE

	expect_warning(inf$.__enclos_env__$private$warn_bootstrap_fallback_once(), "falling back to bootstrap")
	expect_no_warning(inf$.__enclos_env__$private$warn_bootstrap_fallback_once())
	expect_true(isTRUE(inf$.__enclos_env__$private$cached_values$warned_bootstrap_se_unavailable))
})
