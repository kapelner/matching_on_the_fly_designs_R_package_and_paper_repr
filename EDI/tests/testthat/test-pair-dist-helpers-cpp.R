library(testthat)
library(EDI)

# Regression-protection tests for the C++ helpers in EDI/src/pair_dist_helpers.cpp.
#
# These helpers were rewritten from Rcpp accessor loops to raw-pointer loops
# in the additional-sweep phase and retained (4x speedup on pair_distance_matrix).
# The tests guard correctness of the output against pure-R reference implementations.
#
# Functions tested:
#
#   compute_pair_averages_cpp(X, m_vec, m)
#     -- group-wise column means of X, one row per group ID
#
#   compute_pair_distance_matrix_cpp(pair_avg, weights)
#     -- weighted squared Euclidean distance matrix over pair averages
#     -- diagonal = Inf; symmetric

# ---------------------------------------------------------------------------
# Pure-R reference implementations
# ---------------------------------------------------------------------------

pair_avg_r <- function(X, m_vec, m) {
	out <- matrix(NA_real_, nrow = m, ncol = ncol(X))
	for (id in seq_len(m)) {
		rows <- which(m_vec == id)
		if (length(rows) > 0L) out[id, ] <- colMeans(X[rows, , drop = FALSE])
	}
	out
}

pair_dist_r <- function(pair_avg, weights) {
	m <- nrow(pair_avg)
	p <- ncol(pair_avg)
	use_w <- length(weights) == p
	D <- matrix(0, m, m)
	diag(D) <- Inf
	for (i in seq_len(m - 1L)) {
		for (j in seq(i + 1L, m)) {
			d <- sum((pair_avg[i, ] - pair_avg[j, ])^2 *
			             if (use_w) weights else rep(1, p))
			D[i, j] <- d
			D[j, i] <- d
		}
	}
	D
}

# ===========================================================================
# compute_pair_averages_cpp
# ===========================================================================

test_that("compute_pair_averages_cpp matches R reference on small exact case", {
	X     <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4L, ncol = 2L)
	m_vec <- c(1L, 1L, 2L, 2L)
	m     <- 2L

	cpp_out <- EDI:::compute_pair_averages_cpp(X, m_vec, m)
	r_out   <- pair_avg_r(X, m_vec, m)

	expect_equal(cpp_out, r_out, tolerance = 1e-14,
	             label = "compute_pair_averages_cpp exact 2-group case")
})

test_that("compute_pair_averages_cpp handles unequal group sizes", {
	set.seed(1)
	n  <- 10L
	p  <- 3L
	m  <- 4L
	X  <- matrix(rnorm(n * p), n, p)
	# Groups of size 1, 2, 3, 4
	m_vec <- c(1L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L, 4L)

	cpp_out <- EDI:::compute_pair_averages_cpp(X, m_vec, m)
	r_out   <- pair_avg_r(X, m_vec, m)

	expect_equal(cpp_out, r_out, tolerance = 1e-12)
})

test_that("compute_pair_averages_cpp returns NA for empty groups", {
	X     <- matrix(c(10, 20, 30, 40), nrow = 2L, ncol = 2L)
	m_vec <- c(1L, 3L)    # group 2 is empty
	m     <- 3L

	cpp_out <- EDI:::compute_pair_averages_cpp(X, m_vec, m)

	expect_true(all(is.na(cpp_out[2L, ])),
	            label = "empty group produces NA row")
	expect_equal(as.numeric(cpp_out[1L, ]), c(10, 30), tolerance = 1e-14)
	expect_equal(as.numeric(cpp_out[3L, ]), c(20, 40), tolerance = 1e-14)
})

test_that("compute_pair_averages_cpp is consistent with R on random medium-sized input", {
	set.seed(42)
	n  <- 60L
	p  <- 5L
	m  <- 20L
	X  <- matrix(rnorm(n * p), n, p)
	m_vec <- sample(seq_len(m), n, replace = TRUE)
	# Ensure every group has at least one observation (re-assign singletons if needed)
	for (id in seq_len(m)) {
		if (!any(m_vec == id)) m_vec[sample(n, 1L)] <- as.integer(id)
	}

	cpp_out <- EDI:::compute_pair_averages_cpp(X, m_vec, m)
	r_out   <- pair_avg_r(X, m_vec, m)

	expect_equal(cpp_out, r_out, tolerance = 1e-10)
})

# ===========================================================================
# compute_pair_distance_matrix_cpp
# ===========================================================================

test_that("compute_pair_distance_matrix_cpp matches R reference on simple 3-pair case", {
	pair_avg <- matrix(c(0, 1, 2,
	                     0, 0, 1), nrow = 3L, ncol = 2L)
	weights  <- c(1, 1)

	cpp_out <- EDI:::compute_pair_distance_matrix_cpp(pair_avg, weights)
	r_out   <- pair_dist_r(pair_avg, weights)

	expect_equal(cpp_out, r_out, tolerance = 1e-14,
	             label = "3-pair unweighted distance matrix")
})

test_that("compute_pair_distance_matrix_cpp diagonal is Inf", {
	set.seed(2)
	pair_avg <- matrix(rnorm(10L * 3L), 10L, 3L)
	weights  <- runif(3L)

	D <- EDI:::compute_pair_distance_matrix_cpp(pair_avg, weights)

	expect_true(all(is.infinite(diag(D))), label = "diagonal entries are Inf")
})

test_that("compute_pair_distance_matrix_cpp output is symmetric", {
	set.seed(3)
	pair_avg <- matrix(rnorm(8L * 4L), 8L, 4L)
	weights  <- runif(4L)

	D <- EDI:::compute_pair_distance_matrix_cpp(pair_avg, weights)

	expect_equal(D, t(D), tolerance = 1e-14, label = "distance matrix is symmetric")
})

test_that("compute_pair_distance_matrix_cpp with weight length != p uses unweighted distances", {
	pair_avg    <- matrix(c(0, 3,
	                        0, 4), nrow = 2L, ncol = 2L)
	weights_bad <- c(1)    # length 1 != p=2, so ignored

	D_no_weight <- EDI:::compute_pair_distance_matrix_cpp(pair_avg, weights_bad)

	# Expected: Euclidean squared distance = 3^2 + 4^2 = 25
	expect_equal(D_no_weight[1L, 2L], 25, tolerance = 1e-14)
	expect_equal(D_no_weight[2L, 1L], 25, tolerance = 1e-14)
})

test_that("compute_pair_distance_matrix_cpp applies weights correctly", {
	pair_avg <- matrix(c(0, 1,
	                     0, 1), nrow = 2L, ncol = 2L)   # diff = (1, 1)
	weights  <- c(2, 3)    # weighted: 2*1^2 + 3*1^2 = 5

	D <- EDI:::compute_pair_distance_matrix_cpp(pair_avg, weights)

	expect_equal(D[1L, 2L], 5, tolerance = 1e-14)
})

test_that("compute_pair_distance_matrix_cpp matches R reference on random medium input", {
	set.seed(77)
	m        <- 15L
	p        <- 4L
	pair_avg <- matrix(rnorm(m * p), m, p)
	weights  <- runif(p)

	cpp_out <- EDI:::compute_pair_distance_matrix_cpp(pair_avg, weights)
	r_out   <- pair_dist_r(pair_avg, weights)

	expect_equal(cpp_out, r_out, tolerance = 1e-10,
	             label = "random medium input: cpp == R reference")
})

# ===========================================================================
# Round-trip: compute_pair_averages_cpp -> compute_pair_distance_matrix_cpp
# ===========================================================================

test_that("pair averages fed into pair distance matrix match pure-R round-trip", {
	set.seed(55)
	n     <- 40L
	p     <- 3L
	m     <- 12L
	X     <- matrix(rnorm(n * p), n, p)
	m_vec <- sample(seq_len(m), n, replace = TRUE)
	for (id in seq_len(m)) {
		if (!any(m_vec == id)) m_vec[sample(n, 1L)] <- as.integer(id)
	}
	weights <- runif(p)

	avg_cpp <- EDI:::compute_pair_averages_cpp(X, m_vec, m)
	D_cpp   <- EDI:::compute_pair_distance_matrix_cpp(avg_cpp, weights)

	avg_r   <- pair_avg_r(X, m_vec, m)
	D_r     <- pair_dist_r(avg_r, weights)

	# NA rows (empty groups) will produce NA distances; ignore those entries
	finite_mask <- is.finite(D_r)
	expect_equal(D_cpp[finite_mask], D_r[finite_mask], tolerance = 1e-10,
	             label = "round-trip: cpp and R agree on non-empty pairs")
})

# ===========================================================================
# Regression fingerprint
#
# The perf report (additional sweep outcome) verified:
#   "upper-triangle sum before: 2087702.3513132515"
#   "upper-triangle sum after:  2087702.3513132515"
#
# We reproduce the same invariant on a locally generated dataset to confirm
# the rewrite did not alter output values.
# ===========================================================================

test_that("compute_pair_distance_matrix_cpp upper-triangle sum is stable (regression fingerprint)", {
	set.seed(2026)
	m        <- 50L
	p        <- 6L
	pair_avg <- matrix(rnorm(m * p, mean = 0, sd = 10), m, p)
	weights  <- rep(1, p)

	D <- EDI:::compute_pair_distance_matrix_cpp(pair_avg, weights)

	upper_sum <- sum(D[upper.tri(D)])
	r_ref     <- sum(pair_dist_r(pair_avg, weights)[upper.tri(matrix(0, m, m))])

	expect_equal(upper_sum, r_ref, tolerance = 1e-6,
	             label = "upper-triangle sum matches pure-R reference (fingerprint)")
})
