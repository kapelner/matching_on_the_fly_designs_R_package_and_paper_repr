context("use_rcpp wiring for KK robust regression classes")

make_kk_robust_design <- function(n_pairs = 40L, n_single = 20L, seed = 20260729) {
	set.seed(seed)
	n <- 2L * n_pairs + n_single
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = "continuous", verbose = FALSE)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
	}
	des$.__enclos_env__$private$m <- c(rep(seq_len(n_pairs), each = 2L), rep(0L, n_single))
	priv <- des$.__enclos_env__$private
	w <- priv$w
	X <- as.data.frame(priv$X)
	eta <- 0.4 - 0.5 * w + 0.3 * X[[1]] - 0.2 * X[[2]]
	y <- eta + rt(n, df = 4) * 0.5
	# a few outliers to actually exercise robust weighting
	y[1:3] <- y[1:3] + 15
	des$add_all_subject_responses(y)
	des
}

test_that("InferenceContinKKRobustRegrOneLik: use_rcpp matches MASS::rlm fallback", {
	skip_if_not_installed("MASS")
	des <- make_kk_robust_design()

	inf_fast <- InferenceContinKKRobustRegrOneLik$new(des, method = "M", use_rcpp = TRUE, verbose = FALSE)
	inf_ref  <- InferenceContinKKRobustRegrOneLik$new(des, method = "M", use_rcpp = FALSE, verbose = FALSE)

	est_fast <- inf_fast$compute_estimate()
	est_ref  <- inf_ref$compute_estimate()
	expect_equal(est_fast, est_ref, tolerance = 1e-3)

	se_fast <- inf_fast$.__enclos_env__$private$get_standard_error()
	se_ref  <- inf_ref$.__enclos_env__$private$get_standard_error()
	if (is.finite(se_ref)) {
		expect_true(is.finite(se_fast))
		expect_equal(se_fast, se_ref, tolerance = 1e-2)
	} else {
		expect_true(!is.finite(se_fast))
	}
})

test_that("InferenceContinKKRobustRegrIVWC: use_rcpp matches MASS::rlm fallback", {
	skip_if_not_installed("MASS")
	des <- make_kk_robust_design()

	inf_fast <- InferenceContinKKRobustRegrIVWC$new(des, method = "M", use_rcpp = TRUE, verbose = FALSE)
	inf_ref  <- InferenceContinKKRobustRegrIVWC$new(des, method = "M", use_rcpp = FALSE, verbose = FALSE)

	est_fast <- inf_fast$compute_estimate()
	est_ref  <- inf_ref$compute_estimate()
	expect_true(is.finite(est_fast))
	expect_true(is.finite(est_ref))
	expect_equal(est_fast, est_ref, tolerance = 1e-3)
})

test_that("InferenceContinKKRobustRegrOneLik: bootstrap-weighted estimate matches between use_rcpp TRUE/FALSE", {
	skip_if_not_installed("MASS")
	des <- make_kk_robust_design()
	n <- des$get_n()

	inf_fast <- InferenceContinKKRobustRegrOneLik$new(des, method = "M", use_rcpp = TRUE, verbose = FALSE)
	inf_ref  <- InferenceContinKKRobustRegrOneLik$new(des, method = "M", use_rcpp = FALSE, verbose = FALSE)
	boot_ctx <- list(row_to_unit = seq_len(n), unit_group_id = rep(1L, n), n_units = n)
	inf_fast$.__enclos_env__$private$current_bayesian_bootstrap_context <- boot_ctx
	inf_ref$.__enclos_env__$private$current_bayesian_bootstrap_context <- boot_ctx

	set.seed(1)
	weights <- runif(n, 0.5, 1.5)

	est_fast <- inf_fast$compute_estimate_with_bootstrap_weights(weights, estimate_only = TRUE)
	est_ref  <- inf_ref$compute_estimate_with_bootstrap_weights(weights, estimate_only = TRUE)

	expect_true(is.finite(est_fast))
	expect_true(is.finite(est_ref))
	expect_equal(est_fast, est_ref, tolerance = 3e-2)
})

test_that("InferenceContinRobustRegr (non-KK): bootstrap-weighted estimate matches between use_rcpp TRUE/FALSE", {
	skip_if_not_installed("MASS")
	set.seed(20260729)
	n <- 80L
	x1 <- rnorm(n)
	w <- rep(c(0, 1), length.out = n)
	eta <- 0.3 - 0.4 * w + 0.2 * x1
	y <- eta + rt(n, df = 4) * 0.5

	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous", verbose = FALSE)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = x1[i]))
	}
	des$overwrite_all_subject_assignments(w)
	des$add_all_subject_responses(y)

	inf_fast <- InferenceContinRobustRegr$new(des, method = "M", use_rcpp = TRUE, verbose = FALSE)
	inf_ref  <- InferenceContinRobustRegr$new(des, method = "M", use_rcpp = FALSE, verbose = FALSE)
	boot_ctx <- list(row_to_unit = seq_len(n), unit_group_id = rep(1L, n), n_units = n)
	inf_fast$.__enclos_env__$private$current_bayesian_bootstrap_context <- boot_ctx
	inf_ref$.__enclos_env__$private$current_bayesian_bootstrap_context <- boot_ctx

	weights <- runif(n, 0.5, 1.5)
	est_fast <- inf_fast$compute_estimate_with_bootstrap_weights(weights)
	est_ref  <- inf_ref$compute_estimate_with_bootstrap_weights(weights)

	expect_true(is.finite(est_fast))
	expect_true(is.finite(est_ref))
	expect_equal(est_fast, est_ref, tolerance = 1e-3)
})
