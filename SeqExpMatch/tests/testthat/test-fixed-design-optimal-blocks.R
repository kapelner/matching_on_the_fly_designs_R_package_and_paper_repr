test_that("FixedDesignOptimalBlocks is gated behind required libraries", {
	required = c("ompr", "ompr.roi", "ROI.plugin.glpk", "randomizr")
	all_installed = all(vapply(required, requireNamespace, logical(1), quietly = TRUE))
	if (!all_installed) {
		expect_error(
			FixedDesignOptimalBlocks$new(response_type = "continuous", K = 2, n = 6, verbose = FALSE),
			"are required for FixedDesignOptimalBlocks"
		)
		return(invisible(NULL))
	}

	n = 8
	X = data.frame(x1 = c(-3, -2.9, -2.8, -2.7, 2.7, 2.8, 2.9, 3), x2 = c(0, 0.1, -0.1, 0.05, 5, 5.1, 4.9, 5.05))
	for (dist_name in c("euclidean", "sum_abs_diff", "mahal")) {
		des = FixedDesignOptimalBlocks$new(response_type = "continuous", K = 2, dist = dist_name, n = n, verbose = FALSE)
		des$add_all_subjects_to_experiment(X)
		des$assign_w_to_all_subjects()
		w = des$get_w()
		expect_length(w, n)
		expect_equal(sum(w), n / 2)
		W = des$draw_ws_according_to_design(r = 3)
		expect_equal(dim(W), c(n, 3))
	}
})

test_that("FixedDesignOptimalBlocks errors when roughly equal block sizes are infeasible", {
	required = c("ompr", "ompr.roi", "ROI.plugin.glpk", "randomizr")
	skip_if_not(all(vapply(required, requireNamespace, logical(1), quietly = TRUE)))

	expect_error(
		FixedDesignOptimalBlocks$new(response_type = "continuous", K = 7, n = 6, verbose = FALSE),
		"Cannot partition 6 subjects into 7 roughly equal blocks"
	)
})
