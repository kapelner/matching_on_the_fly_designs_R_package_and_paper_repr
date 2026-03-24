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

test_that("FixedDesignCluster works", {
	n = 12
	X = data.frame(cluster = rep(1:4, each = 3), x1 = rnorm(n))
	des = FixedDesignCluster$new(cluster_col = "cluster", n = n, verbose = FALSE)
	for (i in 1:n) {
		des$add_subject(X[i, , drop = FALSE])
	}
	
	des$randomize()
	w = des$get_w()
	expect_length(w, n)
	
	# Check all members of a cluster have same assignment
	for (c in 1:4){
		idxs = which(X$cluster == c)
		expect_true(all(w[idxs] == w[idxs[1]]))
	}

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	for (j in 1:5){
		wj = W[, j]
		for (c in 1:4){
			idxs = which(X$cluster == c)
			expect_true(all(wj[idxs] == wj[idxs[1]]))
		}
	}
})

test_that("FixedDesignFactorial works", {
	n = 12
	# 2x2 factorial: 4 combinations
	des = FixedDesignFactorial$new(factors = list(A=2, B=2), n = n, verbose = FALSE)
	for (i in 1:n) {
		des$add_subject(data.frame(x1 = i))
	}
	
	des$randomize()
	w = des$get_w()
	expect_length(w, n)
	
	# Check balance: each combination (1, 2, 3, 4) should appear 3 times
	expect_equal(as.numeric(table(w)), rep(3, 4))

	# Check data frame output
	wf = des$get_w_factorial()
	expect_equal(nrow(wf), n)
	expect_equal(ncol(wf), 2)
	expect_true(all(wf$Var1 %in% 1:2))
	expect_true(all(wf$Var2 %in% 1:2))

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	for (j in 1:5){
		expect_equal(as.numeric(table(W[, j])), rep(3, 4))
	}
})

test_that("FixedDesignDOptimal works", {
	n = 20
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = FixedDesignDOptimal$new(n = n, verbose = FALSE)
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

test_that("FixedDesignAOptimal works", {
	n = 20
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = FixedDesignAOptimal$new(n = n, verbose = FALSE)
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
