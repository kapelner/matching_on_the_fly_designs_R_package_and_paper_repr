test_that("DesignFixedBinaryMatch works", {
	n = 10
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = DesignFixedBinaryMatch$new(n = n, response_type = "continuous", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	
	des$assign_w_to_all_subjects()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n/2)
	
	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n/2))
})

test_that("DesignFixedRerandomization works", {
	n = 10
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	# Use a very high cutoff so it finds it quickly
	des = DesignFixedRerandomization$new(response_type = "continuous", n = n, obj_val_cutoff = 100, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	
	des$assign_w_to_all_subjects()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n/2)

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n/2))
})

test_that("DesignFixedGreedy works", {
	n = 10
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = DesignFixedGreedy$new(response_type = "continuous", n = n, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	
	des$assign_w_to_all_subjects()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n/2)

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n/2))
})

test_that("DesignFixediBCRD works", {
	n = 10
	des = DesignFixediBCRD$new(response_type = "continuous", n = n, verbose = FALSE)
	expect_identical(des$get_block_ids(), rep(1L, n))
	des$add_all_subjects_to_experiment(data.frame(x1 = 1:n))
	
	des$assign_w_to_all_subjects()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n/2)

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n/2))
})

test_that("DesignFixedBlocking works", {
	n = 12
	X = data.frame(strata = rep(c("A", "B"), each = 6), x1 = rnorm(n))
	des = DesignFixedBlocking$new(response_type = "continuous", strata_cols = "strata", n = n, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	
	des$assign_w_to_all_subjects()
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

test_that("DesignFixedCluster works", {
	n = 12
	X = data.frame(cluster = rep(1:4, each = 3), x1 = rnorm(n))
	des = DesignFixedCluster$new(response_type = "continuous", cluster_col = "cluster", n = n, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	
	des$assign_w_to_all_subjects()
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

test_that("DesignFixedFactorial works", {
	n = 12
	# 2x2 factorial: 4 combinations
	des = DesignFixedFactorial$new(response_type = "continuous", factors = list(A=2, B=2), n = n, verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = 1:n))
	
	des$assign_w_to_all_subjects()
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

test_that("DesignFixedDOptimal works", {
	n = 20
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = DesignFixedDOptimal$new(response_type = "continuous", n = n, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	
	des$assign_w_to_all_subjects()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n/2)

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n/2))
})

test_that("DesignFixedAOptimal works", {
	n = 20
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = DesignFixedAOptimal$new(response_type = "continuous", n = n, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	
	des$assign_w_to_all_subjects()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n/2)

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n/2))
})
