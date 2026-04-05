test_that("Different designs work with continuous response", {
	n <- 20
	designs <- list(
	Bernoulli = DesignSeqOneByOneBernoulli$new(response_type = "continuous", n = n),
	Efron = DesignSeqOneByOneEfron$new(response_type = "continuous", n = n),
	Atkinson = DesignSeqOneByOneAtkinson$new(response_type = "continuous", n = n),
	iBCRD = DesignSeqOneByOneiBCRD$new(response_type = "continuous", n = n),
	KK14 = DesignSeqOneByOneKK14$new(response_type = "continuous", n = n),
	KK21 = DesignSeqOneByOneKK21$new(response_type = "continuous", n = n),
	KK21stepwise = DesignSeqOneByOneKK21stepwise$new(response_type = "continuous", n = n)
	)

	for (name in names(designs)) {
	des <- designs[[name]]
	set.seed(1)
	for (i in 1:n) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
	}
	add_all_subject_responses_seq(des, rnorm(n))

	# Check if we can run a simple inference
	inf <- InferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
	expect_true(is.numeric(inf$compute_treatment_estimate()), info = paste("Design:", name))
	}
})

test_that("KK14 with Morrison works", {
	n <- 20
	p <- 2
	des <- DesignSeqOneByOneKK14$new(response_type = "continuous", n = n, morrison = TRUE, p = p)
	set.seed(1)
	for (i in 1:n) {
	des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
	}
	add_all_subject_responses_seq(des, rnorm(n))

	inf <- InferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
	expect_true(is.numeric(inf$compute_treatment_estimate()))
})

test_that("Multivariate inference works", {
	n <- 40
	des <- DesignSeqOneByOneBernoulli$new(response_type = "continuous", n = n)
	set.seed(1)
	X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
	for (i in 1:n) {
	des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	}
	add_all_subject_responses_seq(des, rnorm(n))

	# Continuous MultOLS
	inf_cont <- InferenceContinMultOLS$new(des, verbose = FALSE)
	expect_true(is.numeric(inf_cont$compute_treatment_estimate()))

	# Incidence MultiLogRegr
	des_inc <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "incidence")
	for (i in 1:n) {
	des_inc$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	}
	add_all_subject_responses_seq(des_inc, rbinom(n, 1, 0.5))
	inf_inc <- InferenceIncidMultiLogRegr$new(des_inc, verbose = FALSE)
	expect_true(is.numeric(inf_inc$compute_treatment_estimate()))

	# Count MultiNegBinRegr
	des_count <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "count")
	for (i in 1:n) {
	des_count$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	}
	add_all_subject_responses_seq(des_count, rpois(n, 5))
	inf_count <- InferenceCountMultiNegBinRegr$new(des_count, verbose = FALSE)
	expect_true(is.numeric(inf_count$compute_treatment_estimate()))
})

test_that("Generated permutations support randomization p-values", {
	n <- 30
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous", verbose = FALSE)
	set.seed(1)
	for (i in 1:n) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
	}
	add_all_subject_responses_seq(des, rnorm(n))

	inf <- InferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
	perms <- inf$.__enclos_env__$private$generate_permutations(64)

	expect_true(is.character(attr(perms, "sig")))
	expect_length(attr(perms, "sig"), 1)

	pval <- inf$compute_two_sided_pval_for_treatment_effect_rand(
		r = 64,
		show_progress = FALSE,
		permutations = perms
	)
	expect_true(is.numeric(pval))
	expect_true(is.finite(pval))
	expect_true(pval >= 0 && pval <= 1)
})

test_that("KK designs with KK inference work", {
	n <- 30
	designs <- list(
	KK14 = DesignSeqOneByOneKK14$new(response_type = "continuous", n = n),
	KK21 = DesignSeqOneByOneKK21$new(response_type = "continuous", n = n),
	KK21stepwise = DesignSeqOneByOneKK21stepwise$new(response_type = "continuous", n = n)
	)

	inf_classes <- list(
	KK14 = InferenceBaiAdjustedTKK14,
	KK21 = InferenceBaiAdjustedTKK21,
	KK21stepwise = InferenceBaiAdjustedTKK21
	)

	for (name in names(designs)) {
	des <- designs[[name]]
	set.seed(1)
	for (i in 1:n) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
	}
	add_all_subject_responses_seq(des, rnorm(n))

	inf <- inf_classes[[name]]$new(des, verbose = FALSE)
	expect_true(is.numeric(inf$compute_treatment_estimate()), info = paste("KK Design/Inf:", name))
	}
})

test_that("Debug randomization does not poison the cached p-value path", {
	n <- 24
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous", verbose = FALSE)
	set.seed(11)
	for (i in 1:n) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
	}
	add_all_subject_responses_seq(des, rnorm(n))

	inf <- InferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
	perms <- inf$.__enclos_env__$private$generate_permutations(32)

	debug_rand <- inf$approximate_randomization_distribution_beta_hat_T(
		r = 32,
		show_progress = FALSE,
		permutations = perms,
		debug = TRUE
	)
	pval_after_debug <- inf$compute_two_sided_pval_for_treatment_effect_rand(
		r = 32,
		show_progress = FALSE,
		permutations = perms
	)

	inf_fresh <- InferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
	pval_fresh <- inf_fresh$compute_two_sided_pval_for_treatment_effect_rand(
		r = 32,
		show_progress = FALSE,
		permutations = perms
	)

	expect_true(is.list(debug_rand))
	expect_true(all(is.finite(debug_rand$values)))
	expect_true(is.numeric(pval_after_debug))
	expect_true(is.finite(pval_after_debug))
	expect_true(pval_after_debug >= 0 && pval_after_debug <= 1)
	expect_equal(pval_after_debug, pval_fresh, tolerance = 1e-12)
})



test_that("Fast randomization worker does not mutate shared proportion design state", {
	skip_if_not_installed("betareg")
	n <- 30
	des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "proportion", verbose = FALSE)
	set.seed(19)
	X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
	for (i in 1:n) {
		des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	}
	y <- pmin(0.98, pmax(0.02, plogis(0.5 * X$x1 - 0.3 * X$x2 + rnorm(n, sd = 0.2))))
	add_all_subject_responses_seq(des, y)

	inf <- InferencePropMultiBetaRegr$new(des, verbose = FALSE)
	perms <- inf$.__enclos_env__$private$generate_permutations(19)
	y_before <- des$.__enclos_env__$private$y

	pval <- inf$compute_two_sided_pval_for_treatment_effect_rand(
		r = 19,
		delta = 0.5,
		transform_responses = "logit",
		show_progress = FALSE,
		permutations = perms
	)
	y_after <- des$.__enclos_env__$private$y

	expect_true(is.numeric(pval))
	expect_true(is.finite(pval))
	expect_equal(as.numeric(y_after), as.numeric(y_before), tolerance = 0)
	expect_true(all(y_after >= 0 & y_after <= 1))
})
