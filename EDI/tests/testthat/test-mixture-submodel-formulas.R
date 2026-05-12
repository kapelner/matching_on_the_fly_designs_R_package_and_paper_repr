test_that("auxiliary formulas control zero-inflated and hurdle count submodels", {
	make_count_design = function() {
		set.seed(101)
		des = DesignSeqOneByOneBernoulli$new(
			n = 60L,
			response_type = "count"
		)
		x1 = rnorm(60)
		x2 = rep(c(0, 1), 30)
		for (i in seq_len(60)) {
			des$add_one_subject_to_experiment_and_assign(data.frame(x1 = x1[i], x2 = x2[i]))
		}
		y = rpois(60, lambda = exp(0.2 + 0.3 * x1 + 0.4 * x2))
		y[sample.int(60, 20)] = 0
		des$add_all_subject_responses(y)
		des
	}

	des = make_count_design()

	zip_treat_only = InferenceCountZeroInflatedPoisson$new(
		des,
		model_formula = ~ x1,
		model_formula_zero = ~ 1,
		use_rcpp = TRUE,
		verbose = FALSE,
		optimization_alg = "lbfgs"
	)
	zip_treat_only$compute_estimate(estimate_only = TRUE)
	zip_ctx_1 = zip_treat_only$.__enclos_env__$private$cached_values$likelihood_test_context
	expect_equal(ncol(zip_ctx_1$Xzi), 2L)

	zip_full = InferenceCountZeroInflatedPoisson$new(
		des,
		model_formula = ~ x1,
		model_formula_zero = ~ x1 + x2,
		use_rcpp = TRUE,
		verbose = FALSE,
		optimization_alg = "lbfgs"
	)
	zip_full$compute_estimate(estimate_only = TRUE)
	zip_ctx_2 = zip_full$.__enclos_env__$private$cached_values$likelihood_test_context
	expect_gt(ncol(zip_ctx_2$Xzi), ncol(zip_ctx_1$Xzi))

	hurdle_treat_only = InferenceCountHurdlePoisson$new(
		des,
		model_formula = ~ x1,
		model_formula_hurdle = ~ 1,
		use_rcpp = TRUE,
		verbose = FALSE,
		optimization_alg = "lbfgs"
	)
	hurdle_treat_only$compute_estimate(estimate_only = TRUE)
	hurdle_ctx_1 = hurdle_treat_only$.__enclos_env__$private$cached_values$likelihood_test_context
	expect_equal(ncol(hurdle_ctx_1$Xzi), 2L)

	hurdle_full = InferenceCountHurdlePoisson$new(
		des,
		model_formula = ~ x1,
		model_formula_hurdle = ~ x1 + x2,
		use_rcpp = TRUE,
		verbose = FALSE,
		optimization_alg = "lbfgs"
	)
	hurdle_full$compute_estimate(estimate_only = TRUE)
	hurdle_ctx_2 = hurdle_full$.__enclos_env__$private$cached_values$likelihood_test_context
	expect_gt(ncol(hurdle_ctx_2$Xzi), ncol(hurdle_ctx_1$Xzi))

	hnb_treat_only = InferenceCountHurdleNegBin$new(
		des,
		model_formula = ~ x1,
		model_formula_hurdle = ~ 1,
		verbose = FALSE
	)
	hnb_treat_only$compute_estimate(estimate_only = TRUE)
	hnb_ctx_1 = hnb_treat_only$.__enclos_env__$private$cached_values$count_likelihood_context
	expect_equal(ncol(hnb_ctx_1$X_hurdle), 2L)

	hnb_full = InferenceCountHurdleNegBin$new(
		des,
		model_formula = ~ x1,
		model_formula_hurdle = ~ x1 + x2,
		verbose = FALSE
	)
	hnb_full$compute_estimate(estimate_only = TRUE)
	hnb_ctx_2 = hnb_full$.__enclos_env__$private$cached_values$count_likelihood_context
	expect_gt(ncol(hnb_ctx_2$X_hurdle), ncol(hnb_ctx_1$X_hurdle))
})

test_that("zero-one inflated beta uses the separate zero-one formula", {
	set.seed(202)
	des = DesignSeqOneByOneBernoulli$new(
		n = 60L,
		response_type = "proportion"
	)
	x1 = rnorm(60)
	x2 = rep(c(0, 1), 30)
	for (i in seq_len(60)) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = x1[i], x2 = x2[i]))
	}
	y = plogis(-0.2 + 0.4 * x1)
	y[seq(1, 60, by = 6)] = 0
	y[seq(4, 60, by = 6)] = 1
	des$add_all_subject_responses(y)

	zoib_treat_only = InferencePropZeroOneInflatedBetaRegr$new(
		des,
		model_formula = ~ x1,
		model_formula_zero_one = ~ 1,
		verbose = FALSE
	)
	zoib_treat_only$compute_estimate(estimate_only = TRUE)
	ctx_1 = zoib_treat_only$.__enclos_env__$private$cached_values$likelihood_test_context
	expect_equal(ncol(ctx_1$X_zero_one), 2L)

	zoib_full = InferencePropZeroOneInflatedBetaRegr$new(
		des,
		model_formula = ~ x1,
		model_formula_zero_one = ~ x1 + x2,
		verbose = FALSE
	)
	zoib_full$compute_estimate(estimate_only = TRUE)
	ctx_2 = zoib_full$.__enclos_env__$private$cached_values$likelihood_test_context
	expect_gt(ncol(ctx_2$X_zero_one), ncol(ctx_1$X_zero_one))
})
