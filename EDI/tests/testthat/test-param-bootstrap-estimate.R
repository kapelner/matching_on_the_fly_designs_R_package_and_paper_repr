library(testthat)
library(EDI)

make_param_boot_logit_design <- function(seed = 20260723L, n = 150L){
	set.seed(seed)
	x1 <- rnorm(n)
	x2 <- rnorm(n)
	w <- rep(c(1, 0), length.out = n)
	p <- plogis(-0.4 + 0.8 * w + 0.35 * x1 - 0.25 * x2)
	y <- rbinom(n, 1, p)

	des <- EDI:::DesignFixed$new(n = n, response_type = "incidence", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1, x2 = x2))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)
	des
}

make_param_boot_cont_ratio_design <- function(seed = 20260724L, n = 96L){
	set.seed(seed)
	x1 <- rnorm(n)
	x2 <- rnorm(n)
	w <- rep(c(1, 0), length.out = n)
	eta <- 0.45 * w + 0.25 * x1 - 0.15 * x2
	cut_1 <- plogis(-1.1 - eta)
	cut_2 <- plogis(-0.1 - eta)
	cut_3 <- plogis(0.9 - eta)
	u <- runif(n)
	y <- ifelse(u <= cut_1, 1L, ifelse(u <= cut_2, 2L, ifelse(u <= cut_3, 3L, 4L)))

	des <- EDI:::DesignFixed$new(n = n, response_type = "ordinal", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = x1, x2 = x2))
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)
	des
}

test_that("compute_param_bootstrap_estimate is finite and deterministic under a fixed seed", {
	des <- make_param_boot_logit_design()

	inf1 <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf1$set_seed(555)
	inf1$num_cores <- 1L

	inf2 <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf2$set_seed(555)
	inf2$num_cores <- 1L

	inf3 <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf3$set_seed(556)
	inf3$num_cores <- 1L

	est1 <- inf1$compute_param_bootstrap_estimate(B = 25L, show_progress = FALSE)
	est2 <- inf2$compute_param_bootstrap_estimate(B = 25L, show_progress = FALSE)
	est3 <- inf3$compute_param_bootstrap_estimate(B = 25L, show_progress = FALSE)

	expect_true(is.finite(est1))
	expect_equal(est1, est2, tolerance = 0)
	expect_false(isTRUE(all.equal(est1, est3)))
})

test_that("continuation-ratio parametric bootstrap exposes full-parameter refits", {
	des <- make_param_boot_cont_ratio_design()
	inf <- InferenceOrdinalContRatioRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf$set_seed(42)
	inf$num_cores <- 1L

	est <- inf$compute_param_bootstrap_estimate(
		B = 7L,
		min_number_usable_samples = 1L,
		max_attempts_per_replicate = 5L,
		show_progress = FALSE
	)
	pval <- inf$compute_param_bootstrap_pval(
		delta = 0,
		B = 7L,
		min_number_usable_samples = 1L,
		max_attempts_per_replicate = 5L,
		show_progress = FALSE
	)
	ci <- inf$compute_param_bootstrap_confidence_interval(
		alpha = 0.2,
		B = 7L,
		min_number_usable_samples = 1L,
		max_attempts_per_replicate = 5L,
		show_progress = FALSE
	)
	diag <- inf$get_last_param_bootstrap_estimate_diagnostics()

	expect_true(is.finite(est))
	expect_true(is.finite(pval))
	expect_true(pval >= 0 && pval <= 1)
	expect_true(all(is.finite(ci)))
	expect_true(ci[[1]] <= ci[[2]])
	expect_equal(diag$n_failure, 0L)
	expect_equal(diag$n_success, 7L)
})

test_that("compute_param_bootstrap_estimate satisfies the bias-correction reconciliation identity", {
	des <- make_param_boot_logit_design(seed = 20260724L, n = 130L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf$set_seed(909)
	inf$num_cores <- 1L

	est <- inf$compute_param_bootstrap_estimate(B = 31L, show_progress = FALSE)
	diag <- inf$get_last_param_bootstrap_estimate_diagnostics()

	expect_true(is.list(diag))
	expect_equal(diag$B, 31L)
	expect_true(is.finite(diag$raw_estimate))
	expect_equal(diag$n_success + diag$n_failure, diag$B)
	expect_true(diag$n_success >= 1L)

	finite_reps <- diag$replicate_estimates[is.finite(diag$replicate_estimates)]
	expected <- 2 * diag$raw_estimate - mean(finite_reps)
	expect_equal(est, expected, tolerance = 1e-12)
})

test_that("parametric bootstrap estimate drops extreme finite refits", {
	des <- make_param_boot_logit_design(seed = 20260725L, n = 120L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf$num_cores <- 1L
	priv_env <- inf$.__enclos_env__$private
	priv_env$param_bootstrap_extreme_estimate_threshold <- 10
	draws <- c(0.2, 1e9, 0.4)
	i <- 0L
	unlockBinding("simulate_under_lik_null", priv_env)
	assign("simulate_under_lik_null", function(spec, delta, null_fit){
		i <<- i + 1L
		b <- rep.int(0, length(spec$full_fit$b))
		b[spec$j] <- draws[[i]]
		list(full_fit = list(b = b))
	}, envir = priv_env)
	lockBinding("simulate_under_lik_null", priv_env)

	est <- inf$compute_param_bootstrap_estimate(
		B = length(draws),
		min_number_usable_samples = 2L,
		max_attempts_per_replicate = 1L,
		show_progress = FALSE
	)
	diag <- inf$get_last_param_bootstrap_estimate_diagnostics()

	expect_true(is.finite(est))
	expect_equal(diag$n_extreme, 1L)
	expect_equal(diag$n_success, 2L)
	expect_equal(diag$n_failure, 1L)
	expect_equal(diag$replicate_estimates, c(0.2, NA_real_, 0.4))
	expect_equal(est, 2 * diag$raw_estimate - mean(c(0.2, 0.4)), tolerance = 1e-12)
})

test_that("compute_param_bootstrap_estimate errors when the family does not support it", {
	des <- make_param_boot_logit_design(seed = 20260725L, n = 90L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	priv_env <- inf$.__enclos_env__$private
	unlockBinding("supports_lik_ratio_param_bootstrap", priv_env)
	assign("supports_lik_ratio_param_bootstrap", function() FALSE, envir = priv_env)
	lockBinding("supports_lik_ratio_param_bootstrap", priv_env)

	expect_error(
		inf$compute_param_bootstrap_estimate(B = 10L, show_progress = FALSE),
		"does not support parametric-bootstrap"
	)
})

test_that("compute_param_bootstrap_confidence_interval returns a finite, ordered, reproducible interval", {
	des <- make_param_boot_logit_design(seed = 20260726L, n = 140L)

	inf1 <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf1$set_seed(321)
	inf1$num_cores <- 1L

	inf2 <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf2$set_seed(321)
	inf2$num_cores <- 1L

	ci1 <- inf1$compute_param_bootstrap_confidence_interval(alpha = 0.2, B = 41L, show_progress = FALSE)
	ci2 <- inf2$compute_param_bootstrap_confidence_interval(alpha = 0.2, B = 41L, show_progress = FALSE)

	expect_true(all(is.finite(ci1)))
	expect_length(ci1, 2)
	expect_true(ci1[[1]] <= ci1[[2]])
	expect_equal(ci1, ci2, tolerance = 0)
})

test_that("compute_param_bootstrap_confidence_interval satisfies the reflected-quantile reconciliation identity", {
	des <- make_param_boot_logit_design(seed = 20260727L, n = 130L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf$set_seed(707)
	inf$num_cores <- 1L

	ci <- inf$compute_param_bootstrap_confidence_interval(alpha = 0.1, B = 51L, show_progress = FALSE)
	diag <- inf$get_last_param_bootstrap_estimate_diagnostics()

	expect_true(is.list(diag))
	expect_equal(diag$B, 51L)

	finite_reps <- diag$replicate_estimates[is.finite(diag$replicate_estimates)]
	expected_lo <- 2 * diag$raw_estimate - stats::quantile(finite_reps, 0.95, names = FALSE)
	expected_hi <- 2 * diag$raw_estimate - stats::quantile(finite_reps, 0.05, names = FALSE)

	expect_equal(unname(ci[[1]]), expected_lo, tolerance = 1e-12)
	expect_equal(unname(ci[[2]]), expected_hi, tolerance = 1e-12)
})

test_that("compute_param_bootstrap_confidence_interval errors when the family does not support it", {
	des <- make_param_boot_logit_design(seed = 20260728L, n = 90L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	priv_env <- inf$.__enclos_env__$private
	unlockBinding("supports_lik_ratio_param_bootstrap", priv_env)
	assign("supports_lik_ratio_param_bootstrap", function() FALSE, envir = priv_env)
	lockBinding("supports_lik_ratio_param_bootstrap", priv_env)

	expect_error(
		inf$compute_param_bootstrap_confidence_interval(alpha = 0.05, B = 10L, show_progress = FALSE),
		"does not support parametric-bootstrap"
	)
})

test_that("compute_param_bootstrap_pval is a finite probability, deterministic under a fixed seed", {
	des <- make_param_boot_logit_design(seed = 20260729L, n = 140L)

	inf1 <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf1$set_seed(111)
	inf1$num_cores <- 1L

	inf2 <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf2$set_seed(111)
	inf2$num_cores <- 1L

	p1 <- inf1$compute_param_bootstrap_pval(delta = 0, B = 41L, show_progress = FALSE)
	p2 <- inf2$compute_param_bootstrap_pval(delta = 0, B = 41L, show_progress = FALSE)

	expect_true(is.finite(p1))
	expect_true(p1 >= 0 && p1 <= 1)
	expect_equal(p1, p2, tolerance = 0)
})

test_that("compute_param_bootstrap_pval satisfies the reflected-empirical-CDF reconciliation identity", {
	des <- make_param_boot_logit_design(seed = 20260730L, n = 130L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	inf$set_seed(222)
	inf$num_cores <- 1L

	delta <- 0.1
	p <- inf$compute_param_bootstrap_pval(delta = delta, B = 51L, show_progress = FALSE)
	diag <- inf$get_last_param_bootstrap_estimate_diagnostics()

	expect_true(is.list(diag))
	expect_equal(diag$B, 51L)

	finite_reps <- diag$replicate_estimates[is.finite(diag$replicate_estimates)]
	t_delta <- 2 * diag$raw_estimate - delta
	n <- length(finite_reps)
	left_tail <- (1 + sum(finite_reps <= t_delta)) / (1 + n)
	right_tail <- (1 + sum(finite_reps >= t_delta)) / (1 + n)
	expected <- min(1, 2 * min(left_tail, right_tail))

	expect_equal(p, expected, tolerance = 1e-12)
})

test_that("compute_param_bootstrap_pval errors when the family does not support it", {
	des <- make_param_boot_logit_design(seed = 20260731L, n = 90L)
	inf <- InferenceIncidLogRegr$new(des, model_formula = ~ x1 + x2, verbose = FALSE)
	priv_env <- inf$.__enclos_env__$private
	unlockBinding("supports_lik_ratio_param_bootstrap", priv_env)
	assign("supports_lik_ratio_param_bootstrap", function() FALSE, envir = priv_env)
	lockBinding("supports_lik_ratio_param_bootstrap", priv_env)

	expect_error(
		inf$compute_param_bootstrap_pval(delta = 0, B = 10L, show_progress = FALSE),
		"does not support parametric-bootstrap"
	)
})
