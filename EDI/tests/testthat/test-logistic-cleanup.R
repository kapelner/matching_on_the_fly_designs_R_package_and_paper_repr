library(testthat)
library(EDI)

context("Logistic regression cleanup")

test_that("fast_logistic_regression_cpp works without smart_cold_start", {
  set.seed(123)
  n <- 100
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- rbinom(n, 1, 0.5)
  
  # This should now work without smart_cold_start argument
  res <- fast_logistic_regression_cpp(X, y)
  expect_equal(length(res$b), p)
  expect_true(res$iterations > 0)
})

test_that("fast_logistic_regression_with_var_cpp works without smart_cold_start", {
  set.seed(123)
  n <- 100
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- rbinom(n, 1, 0.5)
  
  res <- fast_logistic_regression_with_var_cpp(X, y)
  expect_equal(length(res$b), p)
  expect_true(is.matrix(res$fisher_information))
})
