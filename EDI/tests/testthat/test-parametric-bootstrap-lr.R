library(testthat)
library(EDI)

make_param_boot_logit_design <- function(seed = 20260518L, n = 120L){
	set.seed(seed)
	x1 <- rnorm(n)
	x2 <- rnorm(n)
	w <- rep(c(1, 0), length.out = n)
	p <- plogis(-0.4 + 0.8 * w + 0.35 * x1 - 0.25 * x2)
	y <- rbinom(n, 1, p)

	des <- DesignFixed$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1, x2 = x2))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)
	des
}

test_that("parametric bootstrap LR reusable worker path matches generic path", {
	des <- make_param_boot_logit_design()

	inf_worker <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf_worker$set_seed(4242)
	inf_worker$num_cores <- 1L

	inf_generic <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf_generic$set_seed(4242)
	inf_generic$num_cores <- 1L
	inf_generic$.__enclos_env__$private$reusable_bootstrap_worker_enabled <- FALSE

	p_worker <- inf_worker$compute_lik_ratio_bootstrap_two_sided_pval(B = 21L, show_progress = FALSE)
	p_generic <- inf_generic$compute_lik_ratio_bootstrap_two_sided_pval(B = 21L, show_progress = FALSE)
	ci_worker <- inf_worker$compute_lik_ratio_bootstrap_confidence_interval(alpha = 0.2, B = 9L, show_progress = FALSE)
	ci_generic <- inf_generic$compute_lik_ratio_bootstrap_confidence_interval(alpha = 0.2, B = 9L, show_progress = FALSE)

	expect_equal(p_worker, p_generic, tolerance = 1e-12)
	expect_equal(ci_worker, ci_generic, tolerance = 1e-12)
	expect_true(is.finite(p_worker))
	expect_true(all(is.finite(ci_worker)))

	diag_worker <- inf_worker$get_last_param_bootstrap_diagnostics()
	expect_true(is.list(diag_worker))
	expect_equal(diag_worker$B, 9L)
	expect_true(isTRUE(diag_worker$used_reusable_worker))
	expect_equal(diag_worker$n_success + diag_worker$n_failure, diag_worker$B)
	expect_equal(sum(unlist(diag_worker$reason_counts, use.names = FALSE)), diag_worker$B)
	expect_equal(length(diag_worker$replicate_results), diag_worker$B)
	expect_true(all(vapply(diag_worker$replicate_results, function(x) is.list(x) && !is.null(x$reason), logical(1))))
})

test_that("parametric bootstrap LR enforces a minimum usable replicate threshold", {
	des <- make_param_boot_logit_design(seed = 20260519L, n = 80L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf$set_seed(99)

	expect_error(
		inf$compute_lik_ratio_bootstrap_two_sided_pval(B = 4L, min_number_usable_samples = 5L, show_progress = FALSE),
		"B must be at least min_number_usable_samples"
	)
})
