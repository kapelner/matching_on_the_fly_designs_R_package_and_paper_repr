test_that("blocked incidence inference works for CMH and Extended Robins", {
	des <- DesignFixedBlocking$new(
		n = 20,
		response_type = "incidence",
		strata_cols = "x1",
		verbose = FALSE
	)

	des$add_all_subjects_to_experiment(
		data.frame(x1 = factor(rep(1:2, 10), levels = 1:2))
	)
	des$overwrite_all_subject_assignments(rep(c(1, 0), 10))
	des$add_all_subject_responses(rep(c(1, 0, 0, 1), 5))

	inf_cmh <- InferenceIncidCMH$new(des, verbose = FALSE)
	inf_robins <- InferenceIncidExtendedRobins$new(des, verbose = FALSE)

	expect_true(is.finite(inf_cmh$compute_estimate()))
	expect_true(is.finite(inf_robins$compute_estimate()))
	expect_equal(
		inf_cmh$compute_estimate(),
		inf_robins$compute_estimate(),
		tolerance = 1e-12
	)

	expect_equal(length(inf_cmh$compute_asymp_confidence_interval()), 2L)
	expect_equal(length(inf_robins$compute_asymp_confidence_interval()), 2L)
	expect_true(all(is.finite(inf_cmh$compute_asymp_confidence_interval())))
	expect_true(all(is.finite(inf_robins$compute_asymp_confidence_interval())))
	expect_true(is.finite(inf_cmh$compute_asymp_two_sided_pval()))
	expect_true(is.finite(inf_robins$compute_asymp_two_sided_pval()))
})
