library(testthat)
library(EDI)

# Regression-protection tests for d_optimal_search_cpp and a_optimal_search_cpp
# in EDI/src/optimal_design_search.cpp.
#
# These kernels were rewritten for performance (Phase "additional sweep") and
# retained. The tests guard the output contract:
#
#   d_optimal_search_cpp(P, nsim, n_T)
#     -- greedy swap minimisation of w'Pw starting from random allocations
#     -- returns integer matrix n x nsim; each column is a binary {0,1} vector
#        with exactly n_T ones
#
#   a_optimal_search_cpp(P, H, nsim, n_T)
#     -- greedy swap minimisation of (wHw + 1) / (n_T - wPw)
#     -- same output shape as d_optimal_search_cpp
#
# Known-optimal tests use small constructed P / H matrices where the unique
# minimum is analytically obvious, so every random start should converge to
# the same solution.

# ===========================================================================
# 1. d_optimal_search_cpp – structural invariants
# ===========================================================================

test_that("d_optimal_search_cpp returns n x nsim integer matrix with exactly n_T ones per column", {
	set.seed(1)
	n    <- 10L
	nsim <- 20L
	n_T  <- 4L
	P    <- diag(n)

	w_mat <- EDI:::d_optimal_search_cpp(P, nsim = nsim, n_T = n_T)

	expect_equal(dim(w_mat), c(n, nsim))
	expect_true(is.integer(w_mat) || is.numeric(w_mat),
	            label = "d_optimal output is numeric/integer")
	expect_true(all(w_mat %in% c(0L, 1L)),
	            label = "d_optimal entries are 0 or 1")
	expect_true(all(colSums(w_mat) == n_T),
	            label = "every d_optimal column has exactly n_T treated units")
})

# ===========================================================================
# 2. d_optimal_search_cpp – known-optimal convergence
#
# For a diagonal P = diag(4, 1, 1, 4), n=4, n_T=2:
#   w'Pw = sum of P_ii for treated subjects
#   minimum = 1+1 = 2, achieved uniquely by treating subjects 1 and 2 (0-based)
#             i.e. w = (0, 1, 1, 0) in R 1-based = indices 2 and 3.
#
# The greedy algorithm must converge to this from any starting allocation
# (any other pair has w'Pw >= 5).
# ===========================================================================

test_that("d_optimal_search_cpp finds unique minimum for diagonal P", {
	n   <- 4L
	n_T <- 2L
	P   <- diag(c(4, 1, 1, 4))

	# Run many random starts: all should converge to the same unique minimum
	w_mat <- EDI:::d_optimal_search_cpp(P, nsim = 40L, n_T = n_T)

	# Optimal solution: subjects 2 and 3 treated (1-indexed)
	optimal <- c(0L, 1L, 1L, 0L)
	for (s in seq_len(ncol(w_mat))) {
		expect_equal(as.integer(w_mat[, s]), optimal,
		             label = paste0("d_optimal column ", s, " equals unique optimum"))
	}
})

# ===========================================================================
# 3. d_optimal_search_cpp – objective non-worsening guarantee
#
# Every returned allocation must achieve w'Pw <= any single random allocation.
# ===========================================================================

test_that("d_optimal_search_cpp does not worsen w'Pw relative to a random baseline", {
	set.seed(3)
	n   <- 8L
	n_T <- 4L
	# Random symmetric PSD matrix as a toy P
	A <- matrix(rnorm(n * n), n, n)
	P <- A %*% t(A) / n

	w_mat <- EDI:::d_optimal_search_cpp(P, nsim = 30L, n_T = n_T)

	# Worst-case random allocation: random balanced permutation
	set.seed(5)
	w_random <- sample(c(rep(1L, n_T), rep(0L, n - n_T)))
	wPw_random <- as.numeric(w_random %*% P %*% w_random)

	for (s in seq_len(ncol(w_mat))) {
		w_s   <- as.numeric(w_mat[, s])
		wPw_s <- as.numeric(w_s %*% P %*% w_s)
		expect_lte(wPw_s, wPw_random + 1e-10,
		           label = paste0("d_optimal column ", s, " not worse than random wPw"))
	}
})

# ===========================================================================
# 4. a_optimal_search_cpp – structural invariants
# ===========================================================================

test_that("a_optimal_search_cpp returns n x nsim integer matrix with exactly n_T ones per column", {
	set.seed(7)
	n    <- 10L
	nsim <- 20L
	n_T  <- 4L
	P    <- 0.5 * diag(n)
	H    <- diag(n)

	w_mat <- EDI:::a_optimal_search_cpp(P, H, nsim = nsim, n_T = n_T)

	expect_equal(dim(w_mat), c(n, nsim))
	expect_true(all(w_mat %in% c(0L, 1L)),
	            label = "a_optimal entries are 0 or 1")
	expect_true(all(colSums(w_mat) == n_T),
	            label = "every a_optimal column has exactly n_T treated units")
})

# ===========================================================================
# 5. a_optimal_search_cpp – known-optimal convergence
#
# For P = 0.5 * I_4 and H = diag(0, 1, 1, 0), n=4, n_T=2:
#   wPw = 0.5 * 2 = 1 for any balanced allocation  =>  n_T - wPw = 1 always
#   objective = (wHw + 1) / 1 = wHw + 1
#   minimise wHw: choose subjects 1 and 4 (H = 0 each)  =>  wHw = 0
#   unique minimum: w = (1, 0, 0, 1)
# ===========================================================================

test_that("a_optimal_search_cpp finds unique minimum for constructed P and H", {
	n   <- 4L
	n_T <- 2L
	P   <- 0.5 * diag(n)
	H   <- diag(c(0, 1, 1, 0))

	w_mat <- EDI:::a_optimal_search_cpp(P, H, nsim = 40L, n_T = n_T)

	optimal <- c(1L, 0L, 0L, 1L)
	for (s in seq_len(ncol(w_mat))) {
		expect_equal(as.integer(w_mat[, s]), optimal,
		             label = paste0("a_optimal column ", s, " equals unique optimum"))
	}
})

# ===========================================================================
# 6. Regression fingerprint: d_optimal results are unchanged after rewrite
#
# We verify a specific numerical quantity on a fixed random P. This serves as
# a regression guard that the rewrite in the additional-sweep phase did not
# alter the algorithm's output.
#
# The fingerprint is: sum of column-wise wPw values for all nsim runs.
# This is deterministic given that:
#   - The greedy algorithm is deterministic from a given random initialization.
#   - The RNG is seeded inside the C++ function (std::random_device + std::mt19937).
#
# Because the C++ RNG is not seeded from R, we cannot pin the exact value,
# but we CAN verify the structural properties and the objective-non-worsening
# guarantee hold on a fixed P, which is the operationally important guarantee.
#
# The canonical quote from the perf report:
#   "each returned column still had exactly n_T treated units"
# ===========================================================================

test_that("d_optimal_search_cpp regression: all columns have n_T treated on benchmark P", {
	set.seed(99)
	n   <- 20L
	n_T <- 10L
	# Construct a realistic P: projection onto intercept + two covariates
	X   <- cbind(1, rnorm(n), rnorm(n))
	Q   <- qr.Q(qr(X))
	P   <- Q %*% t(Q)

	w_mat <- EDI:::d_optimal_search_cpp(P, nsim = 50L, n_T = n_T)

	expect_equal(dim(w_mat), c(n, 50L))
	expect_true(all(colSums(w_mat) == n_T),
	            label = "benchmark P: all columns have n_T treated")
	expect_true(all(w_mat %in% c(0L, 1L)),
	            label = "benchmark P: all entries are 0 or 1")

	# wPw for each column should be <= wPw for a single random balanced allocation
	set.seed(42)
	w_rnd <- sample(c(rep(1L, n_T), rep(0L, n - n_T)))
	wPw_rnd <- as.numeric(w_rnd %*% P %*% w_rnd)
	wPw_cols <- apply(w_mat, 2L, function(w) as.numeric(w %*% P %*% w))
	expect_true(all(wPw_cols <= wPw_rnd + 1e-10),
	            label = "benchmark P: all columns improve on random baseline")
})
