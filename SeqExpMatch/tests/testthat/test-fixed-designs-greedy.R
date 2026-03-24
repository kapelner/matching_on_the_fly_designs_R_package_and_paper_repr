test_that("FixedDesignBinaryMatch works", {
	n = 10
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = FixedDesignBinaryMatch$new(n = n, response_type = "continuous", verbose = FALSE)
	for (i in 1:n) {
		des$add_subject(X[i, , drop = FALSE])
	}
	
	des$randomize()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n/2)
	
	# Check it's actually matching (simple heuristic: distances within pairs should be small)
	# This is hard to test deterministically, but we can check it runs.
})

test_that("FixedDesignRerandomization works", {
	n = 10
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	# Use a very high cutoff so it finds it quickly
	des = FixedDesignRerandomization$new(n = n, obj_val_cutoff_to_include = 100, verbose = FALSE)
	for (i in 1:n) {
		des$add_subject(X[i, , drop = FALSE])
	}
	
	des$randomize()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n/2)
})

test_that("FixedDesignGreedy works", {
	n = 10
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = FixedDesignGreedy$new(n = n, verbose = FALSE)
	for (i in 1:n) {
		des$add_subject(X[i, , drop = FALSE])
	}
	
	des$randomize()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n/2)
})
