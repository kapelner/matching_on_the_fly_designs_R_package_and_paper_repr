test_that("GreedyExperimentalDesign-backed fixed designs randomize as expected", {
	if (!check_package_installed("GreedyExperimentalDesign")) skip("GreedyExperimentalDesign (or Java) unavailable")

	set.seed(1)
	X = data.frame(
		x1 = rnorm(8),
		x2 = rnorm(8)
	)

	designs = list(
		DesignFixedBinaryMatch$new(
			response_type = "continuous",
			n = 8,
			verbose = FALSE
		),
		DesignFixedRerandomization$new(
			response_type = "continuous",
			n = 8,
			obj_val_cutoff = 100,
			objective = "abs_sum_diff",
			verbose = FALSE
		),
		DesignFixedMatchingGreedyPairSwitching$new(
			response_type = "continuous",
			n = 8,
			objective = "abs_sum_diff",
			max_designs = 4,
			verbose = FALSE
		)
	)

	set_num_cores(1)
	on.exit(unset_num_cores())

	for (des in designs) {
		des$add_all_subjects_to_experiment(X)
		des$assign_w_to_all_subjects()
		w = des$get_w()

		expect_length(w, 8)
		expect_true(all(w %in% 0:1))
		expect_equal(sum(w), 4)
	}
})

test_that("GreedyExperimentalDesign-backed fixed designs honor multicore randomization path", {
	if (!check_package_installed("GreedyExperimentalDesign")) skip("GreedyExperimentalDesign (or Java) unavailable")

	set.seed(1)
	X = data.frame(
		x1 = rnorm(8),
		x2 = rnorm(8)
	)

	designs = list(
		DesignFixedBinaryMatch$new(
			response_type = "continuous",
			n = 8,
			verbose = FALSE
		),
		DesignFixedRerandomization$new(
			response_type = "continuous",
			n = 8,
			obj_val_cutoff = 100,
			objective = "abs_sum_diff",
			verbose = FALSE
		),
		DesignFixedMatchingGreedyPairSwitching$new(
			response_type = "continuous",
			n = 8,
			objective = "abs_sum_diff",
			max_designs = 4,
			verbose = FALSE
		)
	)

	set_num_cores(2)
	on.exit(unset_num_cores())

	for (des in designs) {
		des$add_all_subjects_to_experiment(X)
		expect_silent(des$assign_w_to_all_subjects())
		w = des$get_w()

		expect_length(w, 8)
		expect_true(all(w %in% 0:1))
		expect_equal(sum(w), 4)
	}
})
