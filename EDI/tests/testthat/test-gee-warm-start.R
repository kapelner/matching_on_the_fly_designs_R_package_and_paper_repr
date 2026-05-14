context("GEE Warm Start")
library(EDI)

test_that("GEE Rcpp solver utilizes warm starts correctly", {
	set.seed(20260513)
	n_pairs <- 100L
	n_single <- 50L
	group_id <- c(rep(seq_len(n_pairs), each = 2L), rep(n_pairs + seq_len(n_single), each = 1L))
	n <- length(group_id)
	X <- cbind("(Intercept)" = 1, w = rnorm(n), x1 = rnorm(n), x2 = rnorm(n))
	b <- rnorm(max(group_id), sd = 0.45)
	eta <- drop(X %*% c(0.4, -0.35, 0.2, -0.15) + b[group_id])
	y <- rbinom(n, size = 1L, prob = plogis(eta))

	# First fit: No warm start
	fit1 <- EDI:::gee_pairs_singletons_cpp(X, y, group_id, family_str = "binomial")
	expect_true(fit1$converged)
	expect_gt(fit1$niter, 1L)

	# Second fit: Warm start with same data
	fit2 <- EDI:::gee_pairs_singletons_cpp(X, y, group_id, family_str = "binomial", 
	                                      warm_start_beta = fit1$beta, 
	                                      warm_start_fisher_info = fit1$fisher_information)
	
	expect_true(fit2$converged)
	# Should converge in 0 or 1 iteration if beta is exactly the same and Bread is warm
	expect_lte(fit2$niter, 1L)
	expect_equal(as.numeric(fit2$beta), as.numeric(fit1$beta), tolerance = 1e-10)

	# Third fit: Perturbed warm start
	fit3 <- EDI:::gee_pairs_singletons_cpp(X, y, group_id, family_str = "binomial", 
	                                      warm_start_beta = fit1$beta + 0.1)
	expect_true(fit3$converged)
	expect_equal(as.numeric(fit3$beta), as.numeric(fit1$beta), tolerance = 1e-6)
})

test_that("InferenceAbstractKKGEE class handles warm starts", {
	set.seed(20260513)
	n <- 100
	X <- data.frame(x1 = rnorm(n))
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = "incidence", verbose = FALSE)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	}
	des$add_all_subject_responses(rbinom(n, 1, 0.5))
	
	inf <- InferenceIncidKKGEE$new(des, use_rcpp = TRUE)
	
	# First compute: No warm start cached
	inf$compute_estimate()
	# Warm start should now be cached in private
	warm_beta <- inf$.__enclos_env__$private$fit_warm_start
	expect_false(is.null(warm_beta))
	
	# Re-running to ensure it doesn't crash and returns same result
	est1 <- inf$compute_estimate()
	est2 <- inf$compute_estimate()
	expect_equal(est1, est2)
})
