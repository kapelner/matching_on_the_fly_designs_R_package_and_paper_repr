test_that("SimulationFramework accepts merged design classes and params", {
	set.seed(20260429)
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(
			FixedDesignBernoulli,
			FixedDesignBlocking = list(B_preferred = 2L)
		),
		inference_classes_and_params = list(
			InferenceAllSimpleMeanDiff = list(max_resample_attempts = 10L)
		),
		n = 8L,
		p = 2L,
		Nrep = 1L,
		betaT = 0,
		inference_types_and_params = list(asymp_pval = list(delta = 0)),
		turn_off_asserts_for_speed = FALSE
	)

	sim$run()
	raw <- sim$get_results()
	sm <- sim$summarize()

	expect_true("inference_type" %in% names(raw))
	expect_false("method" %in% names(raw))
	expect_true("inference_type" %in% names(sm))
	expect_false("method" %in% names(sm))
	expect_true("asymp_pval" %in% raw$inference_type)
	expect_true("asymp_pval" %in% sm$inference_type)
	expect_equal(sort(sm$design), sort(c(
		"FixedDesignBernoulli",
		"FixedDesignBlocking (B_preferred=2L)"
	)))
	expect_true("design_params" %in% names(sm))
	expect_true("B_preferred=2L" %in% sm$design_params)
	expect_true("max_resample_attempts=10L" %in% sm$inference_params)
	expect_true("asymp_pval(delta=0)" %in% sm$inference_type_params)
})

test_that("SimulationFramework validates merged inference constructor params", {
	expect_error(
		SimulationFramework$new(
			response_type = "continuous",
			design_classes_and_params = list(FixedDesignBernoulli),
			inference_classes_and_params = list(
				InferenceAllSimpleMeanDiff = list(not_a_constructor_arg = TRUE)
			)
		),
		"not accepted by initialize"
	)
})

test_that("SimulationFramework validates inference type params against function arguments", {
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(FixedDesignBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(asymp_pval = list(not_an_arg = TRUE)),
		n = 8L,
		p = 2L,
		Nrep = 1L,
		betaT = 0
	)

	expect_error(sim$run(), "not accepted by compute_asymp_two_sided_pval")
})
