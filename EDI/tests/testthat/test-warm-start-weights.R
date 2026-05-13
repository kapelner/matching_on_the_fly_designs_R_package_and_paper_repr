library(testthat)
library(EDI)

context("Warm-start weights")

test_that("fast_logistic_regression_cpp supports warm_start_weights", {
  set.seed(123)
  n <- 200
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rbinom(n, 1, 0.5)
  
  # 1. Initial fit
  fit1 <- fast_logistic_regression_cpp(X, y)
  expect_true(fit1$iterations > 1)
  
  # 2. Warm-start with converged weights and beta
  fit2 <- fast_logistic_regression_cpp(X, y, start_beta = fit1$b, warm_start_weights = fit1$w)
  
  # Should converge very quickly
  expect_true(fit2$iterations <= 2)
  expect_equal(fit1$b, fit2$b, tolerance = 1e-7)
})

test_that("fast_poisson_regression_cpp supports warm_start_weights", {
  set.seed(123)
  n <- 200
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rpois(n, 2)
  
  # 1. Initial fit
  fit1 <- fast_poisson_regression_cpp(X, y, optimization_alg = "irls")
  expect_true(fit1$iterations > 1)
  
  # In Poisson IRLS, w = mu (or mu * weights).
  fit2 <- fast_poisson_regression_cpp(X, y, start_beta = fit1$b, warm_start_weights = fit1$mu, optimization_alg = "irls")
  
  expect_true(fit2$iterations <= 2)
  expect_equal(fit1$b, fit2$b, tolerance = 1e-7)
})

test_that("fast_log_binomial_regression_cpp supports warm_start_weights", {
  set.seed(123)
  n <- 200
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  X[,1] <- 1 # Intercept
  # Ensure y=1 is not too common to help log-link convergence
  y <- rbinom(n, 1, 0.2)
  
  # 1. Initial fit
  fit1 <- fast_log_binomial_regression_cpp(X, y)
  
  fit2 <- fast_log_binomial_regression_cpp(X, y, start_beta = fit1$b, warm_start_weights = fit1$working_weights)
  
  expect_true(fit2$iterations <= 2)
  expect_equal(fit1$b, fit2$b, tolerance = 1e-7)
})

test_that("fast_robust_regression_cpp supports warm_start_weights", {
  set.seed(123)
  n <- 100
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- X %*% c(1, 2, 3) + rnorm(n)
  y[1:10] <- y[1:10] + 20 # Outliers
  
  # 1. Initial fit
  fit1 <- fast_robust_regression_cpp(X, y)
  
  fit2 <- fast_robust_regression_cpp(X, y, start_beta = fit1$b, warm_start_weights = fit1$w)
  
  # Robust might take more iterations because the scale is recomputed from the new start_beta residuals.
  # But it should still converge to the same result.
  expect_equal(fit1$b, fit2$b, tolerance = 1e-7)
})
