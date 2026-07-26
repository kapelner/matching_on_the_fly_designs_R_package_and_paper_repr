library(testthat)
library(EDI)

test_that("qr_reduce_preserve_cols_cpp preserves required columns when possible", {
	X = cbind(
		1,
		c(0, 1, 0, 1, 0, 1),
		c(10, 11, 10, 11, 10, 11),
		c(1, 1, 1, 1, 1, 1),
		c(0, 1, 0, 1, 0, 1)
	)

	reduced = EDI:::qr_reduce_preserve_cols_cpp(X, c(1L, 2L))

	expect_true(all(c(1L, 2L) %in% reduced$keep))
	expect_equal(qr(reduced$X_reduced)$rank, ncol(reduced$X_reduced))
	expect_equal(qr(X)$rank, qr(reduced$X_reduced)$rank)
})

test_that("qr_reduce_preserve_cols_cpp drops treatment when it is linearly dependent", {
	X = cbind(
		1,
		rep(1, 6),
		c(1, 2, 3, 4, 5, 6)
	)

	reduced = EDI:::qr_reduce_preserve_cols_cpp(X, c(1L, 2L))

	expect_true(1L %in% reduced$keep)
	expect_false(2L %in% reduced$keep)
	expect_equal(qr(reduced$X_reduced)$rank, ncol(reduced$X_reduced))
})

test_that("rank reduction handles empty and non-finite matrices without native crashes", {
	X_empty_rows = matrix(numeric(0), nrow = 0, ncol = 3)
	colnames(X_empty_rows) = paste0("x", seq_len(ncol(X_empty_rows)))
	reduced_empty_rows = EDI:::drop_linearly_dependent_cols(X_empty_rows)
	expect_equal(dim(reduced_empty_rows$M), c(0L, 3L))
	expect_equal(reduced_empty_rows$js, seq_len(3L))
	expect_equal(EDI:::matrix_rank_cpp(X_empty_rows), 0L)

	X_empty_cols = matrix(numeric(0), nrow = 4, ncol = 0)
	reduced_empty_cols = EDI:::drop_linearly_dependent_cols(X_empty_cols)
	expect_equal(dim(reduced_empty_cols$M), c(4L, 0L))
	expect_equal(reduced_empty_cols$js, integer(0))
	expect_equal(EDI:::matrix_rank_cpp(X_empty_cols), 0L)

	X_nonfinite = cbind(x1 = c(1, 2, NA), x2 = c(2, 4, 6))
	reduced_nonfinite = EDI:::drop_linearly_dependent_cols(X_nonfinite)
	expect_equal(reduced_nonfinite$M, X_nonfinite)
	expect_equal(reduced_nonfinite$js, seq_len(ncol(X_nonfinite)))
})
