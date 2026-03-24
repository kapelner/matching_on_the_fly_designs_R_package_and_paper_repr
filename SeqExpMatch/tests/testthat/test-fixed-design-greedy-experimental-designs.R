test_that("GreedyExperimentalDesign-backed fixed designs randomize as expected", {
	skip_if_not_installed("GreedyExperimentalDesign")

	set.seed(1)
	X = data.frame(
		x1 = rnorm(8),
		x2 = rnorm(8)
	)

	designs = list(
		FixedDesignBinaryMatch$new(
			n = 8,
			num_cores = 1,
			verbose = FALSE
		),
		FixedDesignRerandomization$new(
			n = 8,
			obj_val_cutoff = 100,
			objective = "abs_sum_diff",
			num_cores = 1,
			verbose = FALSE
		),
		FixedDesignMatchingGreedyPairSwitching$new(
			n = 8,
			objective = "abs_sum_diff",
			max_designs = 4,
			num_cores = 1,
			verbose = FALSE
		)
	)

	for (des in designs) {
		for (i in 1:nrow(X)) {
			des$add_subject(X[i, , drop = FALSE])
		}
		des$randomize()
		w = des$get_w()

		expect_length(w, 8)
		expect_true(all(w %in% 0:1))
		expect_equal(sum(w), 4)
	}
})

test_that("GreedyExperimentalDesign-backed fixed designs honor multicore randomization path", {
	skip_if_not_installed("GreedyExperimentalDesign")

	set.seed(1)
	X = data.frame(
		x1 = rnorm(8),
		x2 = rnorm(8)
	)

	designs = list(
		FixedDesignBinaryMatch$new(
			n = 8,
			num_cores = 2,
			verbose = FALSE
		),
		FixedDesignRerandomization$new(
			n = 8,
			obj_val_cutoff = 100,
			objective = "abs_sum_diff",
			num_cores = 2,
			verbose = FALSE
		),
		FixedDesignMatchingGreedyPairSwitching$new(
			n = 8,
			objective = "abs_sum_diff",
			max_designs = 4,
			num_cores = 2,
			verbose = FALSE
		)
	)

	for (des in designs) {
		for (i in 1:nrow(X)) {
			des$add_subject(X[i, , drop = FALSE])
		}
		expect_silent(des$randomize())
		w = des$get_w()

		expect_length(w, 8)
		expect_true(all(w %in% 0:1))
		expect_equal(sum(w), 4)
	}
})
