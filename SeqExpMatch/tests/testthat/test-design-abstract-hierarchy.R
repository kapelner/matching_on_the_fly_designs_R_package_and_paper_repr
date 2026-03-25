test_that("Design hierarchy supports both fixed and sequential designs", {
	seq_des = SeqDesignBernoulli$new(n = 4, response_type = "continuous", verbose = FALSE)
	fixed_des = FixedDesignBernoulli$new(n = 4, response_type = "continuous", verbose = FALSE)

	expect_true(is(seq_des, "Design"))
	expect_true(is(fixed_des, "Design"))
	expect_false(is(seq_des, "FixedDesign"))
	expect_true(is(fixed_des, "FixedDesign"))
})

test_that("plain FixedDesign supports analysis but not redraw-based resampling", {
	des = FixedDesign$new(n = 4, response_type = "continuous", verbose = FALSE)
	for (i in 1:4) {
		des$add_subject(data.frame(x1 = i))
	}
	des$add_all_subject_assignments(c(0, 1, 0, 1))
	des$add_all_subject_responses(c(1, 3, 2, 4))

	expect_false(des$supports_resampling())
	expect_error(des$randomize(), "Plain FixedDesign objects do not support randomization")

	inf = DesignInferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
	expect_equal(inf$compute_treatment_estimate(), 2)
	expect_length(inf$compute_asymp_confidence_interval(), 2)
	expect_true(is.finite(inf$compute_asymp_two_sided_pval_for_treatment_effect()))
	expect_error(
		inf$compute_bootstrap_two_sided_pval(B = 11),
		"Bootstrap inference is not available for plain FixedDesign objects"
	)
	expect_error(
		inf$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test = 11),
		"Randomization inference is not available for plain FixedDesign objects"
	)

	# Test matched pair ID setter
	des$add_all_subject_matched_pair_ids(c(1, 1, 2, 2))
	expect_equal(des$get_t(), 4)
	expect_true(des$is_blocking_design())
	expect_true(des$is_complete_blocking_design())
})
