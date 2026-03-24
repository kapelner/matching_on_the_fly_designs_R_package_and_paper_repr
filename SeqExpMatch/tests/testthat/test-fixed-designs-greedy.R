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
	
	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n/2))
})

test_that("FixedDesignRerandomization works", {
	n = 10
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	# Use a very high cutoff so it finds it quickly
	des = FixedDesignRerandomization$new(n = n, obj_val_cutoff = 100, verbose = FALSE)
	for (i in 1:n) {
		des$add_subject(X[i, , drop = FALSE])
	}
	
	des$randomize()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n/2)

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n/2))
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

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n/2))
})

test_that("FixedDesigniBCRD works", {
	n = 10
	des = FixedDesigniBCRD$new(n = n, verbose = FALSE)
	for (i in 1:n) {
		des$add_subject(data.frame(x1 = i))
	}
	
	des$randomize()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n/2)

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n/2))
})

test_that("FixedDesignBlocking works", {
	n = 12
	X = data.frame(strata = rep(c("A", "B"), each = 6), x1 = rnorm(n))
	des = FixedDesignBlocking$new(strata_cols = "strata", n = n, verbose = FALSE)
	for (i in 1:n) {
		des$add_subject(X[i, , drop = FALSE])
	}
	
	des$randomize()
	w = des$get_w()
	expect_length(w, n)
	# Check balance within strata
	expect_equal(sum(w[1:6]), 3)
	expect_equal(sum(w[7:12]), 3)

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W[1:6, ]) == 3))
	expect_true(all(colSums(W[7:12, ]) == 3))
})
