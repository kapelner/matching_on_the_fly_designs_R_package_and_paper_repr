test_that("DesignFixedBinaryMatch works", {
	n = 10
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = DesignFixedBinaryMatch$new(n = n, response_type = "continuous", verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	
	des$assign_w_to_all_subjects()
	w = des$get_w()
	expect_length(w, n)
	# {0,1} convention (1 = treated, 0 = control; see get_w()'s docstring in
	# design_abstract.R): a balanced n/2-n/2 allocation sums to n/2, not 0.
	expect_equal(sum(w), n / 2)

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n / 2))
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
	expect_equal(sum(w), n / 2)

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n / 2))
})

test_that("DesignFixedRerandomization rejects both obj_val_cutoff and prop_acceptable", {
	expect_error(
		DesignFixedRerandomization$new(
			response_type = "continuous", n = 10,
			obj_val_cutoff = 1.0, prop_acceptable = 0.1, verbose = FALSE
		),
		"Cannot specify both"
	)
})

test_that("DesignFixedRerandomization with obj_val_cutoff returns designs satisfying cutoff", {
	set.seed(1)
	n = 11  # odd n forces pure-R fallback so we control the covariance exactly
	r = 10
	cutoff = 1.0
	X_df = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = DesignFixedRerandomization$new(
		response_type = "continuous", n = n, obj_val_cutoff = cutoff, verbose = FALSE
	)
	des$add_all_subjects_to_experiment(X_df)

	W = des$draw_ws_according_to_design(r = r)
	expect_equal(dim(W), c(n, r))
	expect_true(all(W %in% c(0, 1)))

	X_mat = as.matrix(X_df)
	S_inv = solve(var(X_mat))
	obj_vals = apply(W, 2, function(w) {
		d = colMeans(X_mat[w == 1, , drop = FALSE]) - colMeans(X_mat[w == 0, , drop = FALSE])
		as.numeric(d %*% S_inv %*% d)
	})
	expect_true(all(obj_vals <= cutoff))
})

test_that("DesignFixedRerandomization with prop_acceptable returns correct shape and values", {
	n = 20  # even + balanced -> exercises the GED C++ fast path
	r = 10
	X_df = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = DesignFixedRerandomization$new(
		response_type = "continuous", n = n, prop_acceptable = 0.1, verbose = FALSE
	)
	des$add_all_subjects_to_experiment(X_df)

	W = des$draw_ws_according_to_design(r = r)
	expect_equal(dim(W), c(n, r))
	expect_true(all(W %in% c(0, 1)))
	expect_true(all(colSums(W) == n / 2))
})

test_that("DesignFixedRerandomization with prop_acceptable returns better-balanced designs than pure random", {
	# prop_acceptable = 0.05 with r = 20 draws 400 designs and returns the 20 best.
	# Their worst mahal dist should be well below the median of pure random designs.
	n = 20
	r = 20
	X_df = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = DesignFixedRerandomization$new(
		response_type = "continuous", n = n, prop_acceptable = 0.05, verbose = FALSE
	)
	des$add_all_subjects_to_experiment(X_df)

	W = des$draw_ws_according_to_design(r = r)

	X_mat = as.matrix(X_df)
	S_inv = solve(var(X_mat))
	compute_mahal = function(w) {
		d = colMeans(X_mat[w == 1, , drop = FALSE]) - colMeans(X_mat[w == 0, , drop = FALSE])
		as.numeric(d %*% S_inv %*% d)
	}
	obj_returned = apply(W, 2, compute_mahal)

	n_T = n / 2L
	W_rand = replicate(1000L, sample(c(rep(1L, n_T), rep(0L, n_T))))
	obj_rand = apply(W_rand, 2, compute_mahal)

	# The best 5% should comfortably sit below the random median
	expect_true(max(obj_returned) < median(obj_rand))
})

test_that("DesignFixedGreedy works", {
	n = 10
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = DesignFixedGreedy$new(response_type = "continuous", n = n, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	
	des$assign_w_to_all_subjects()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n / 2)

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n / 2))
})

test_that("DesignFixediBCRD works", {
	n = 10
	des = DesignFixediBCRD$new(response_type = "continuous", n = n, verbose = FALSE)
	expect_identical(des$get_block_ids(), rep(1L, n))
	des$add_all_subjects_to_experiment(data.frame(x1 = 1:n))

	des$assign_w_to_all_subjects()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n / 2)

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n / 2))
})

test_that("DesignFixedBlocking works", {
	n = 12
	X = data.frame(strata = rep(c("A", "B"), each = 6), x1 = rnorm(n))
	des = DesignFixedBlocking$new(response_type = "continuous", strata_cols = "strata", n = n, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	
	des$assign_w_to_all_subjects()
	w = des$get_w()
	expect_length(w, n)
	# Check balance within strata (each stratum has 6 subjects, so a balanced
	# 3-3 split sums to 3 under the {0,1} convention)
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

test_that("DesignFixedFactorial works (two-arm)", {
	n = 12
	# Single two-level factor: 2 combinations (this design is currently two-arm only)
	des = DesignFixedFactorial$new(response_type = "continuous", factors = list(A=2), n = n, verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x1 = 1:n))

	des$assign_w_to_all_subjects()
	w = des$get_w()
	expect_length(w, n)

	# Follows the standard {0,1} public convention: balanced allocation sums to n/2
	expect_true(all(w %in% c(0, 1)))
	expect_equal(sum(w), n / 2)
	expect_equal(as.numeric(table(w)), rep(n / 2, 2))

	# Check data frame output (1-indexed factor levels, independent of the {0,1} w recoding)
	wf = des$get_w_factorial()
	expect_equal(nrow(wf), n)
	expect_equal(ncol(wf), 1)
	expect_true(all(wf$A %in% 1:2))
	# get_w_factorial()'s levels agree with get_w(): 1 <-> level 2, 0 <-> level 1
	expect_equal(wf$A, ifelse(w == 1, 2, 1))

	# Test draw_ws: also follows the {0,1} convention
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	for (j in 1:5){
		expect_true(all(W[, j] %in% c(0, 1)))
		expect_equal(sum(W[, j]), n / 2)
	}
})

test_that("DesignFixedFactorial rejects more than two total combinations", {
	n = 12
	expect_error(
		DesignFixedFactorial$new(response_type = "continuous", factors = list(A=2, B=2), n = n, verbose = FALSE),
		"exactly two total factor-level"
	)
	expect_error(
		DesignFixedFactorial$new(response_type = "continuous", factors = list(A=3), n = n, verbose = FALSE),
		"exactly two total factor-level"
	)
})

test_that("DesignFixedFactorial works correctly with Inference objects", {
	n = 40
	set.seed(1)
	X = data.frame(x1 = rnorm(n))
	des = DesignFixedFactorial$new(response_type = "continuous", factors = list(treatment=2), n = n, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	des$assign_w_to_all_subjects()

	true_effect = 3
	y = rnorm(n) + true_effect * (des$get_w() == 1)
	des$add_all_subject_responses(y)

	inf = InferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
	est = inf$compute_estimate()
	expect_true(is.finite(est))
	expect_equal(est, mean(y[des$get_w() == 1]) - mean(y[des$get_w() == 0]))
	expect_equal(est, true_effect, tolerance = 1.5)

	inf_ols = InferenceContinOLS$new(des, verbose = FALSE)
	est_ols = inf_ols$compute_estimate()
	expect_true(is.finite(est_ols))
	expect_equal(est_ols, true_effect, tolerance = 1.5)
})

test_that("DesignFixedDOptimal works", {
	n = 20
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = DesignFixedDOptimal$new(response_type = "continuous", n = n, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	
	des$assign_w_to_all_subjects()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n / 2)

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n / 2))
})

test_that("DesignFixedAOptimal works", {
	n = 20
	X = data.frame(x1 = rnorm(n), x2 = rnorm(n))
	des = DesignFixedAOptimal$new(response_type = "continuous", n = n, verbose = FALSE)
	des$add_all_subjects_to_experiment(X)
	
	des$assign_w_to_all_subjects()
	w = des$get_w()
	expect_length(w, n)
	expect_equal(sum(w), n / 2)

	# Test draw_ws
	W = des$draw_ws_according_to_design(r = 5)
	expect_equal(dim(W), c(n, 5))
	expect_true(all(colSums(W) == n / 2))
})
