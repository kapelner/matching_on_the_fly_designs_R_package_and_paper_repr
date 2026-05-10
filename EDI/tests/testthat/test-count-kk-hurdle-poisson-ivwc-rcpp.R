test_that("KK hurdle-Poisson IVWC uses the Rcpp matched-pair GLMM path", {
	skip_if_not_installed("GreedyExperimentalDesign")

	set.seed(3)
	n = 40
	X = data.frame(x = rnorm(n))

	des = DesignFixedBinaryMatch$new(response_type = "count", n = n, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()

	y = stats::rpois(n, lambda = exp(0.4 * des$get_w() + 0.1 * X$x))
	y[seq(1, n, 5)] = 0L
	des$add_all_subject_responses(y)

	inf = InferenceCountKKHurdlePoissonIVWC$new(des, verbose = FALSE, use_rcpp = TRUE)

	expect_identical(inf$get_optimization_alg(), "lbfgs")
	expect_true(is.finite(inf$compute_estimate()))
	expect_true(all(is.finite(inf$compute_asymp_confidence_interval())))
})
