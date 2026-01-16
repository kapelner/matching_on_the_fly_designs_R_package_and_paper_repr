library(testthat)
library(SeqExpMatch)
library(survival)

context("Additional tests for optimized functions and argument permutations")

test_that("fast_beta_regression_mle argument permutations", {
  set.seed(1)
  n <- 100
  X <- cbind(1, rnorm(n))
  y <- rbeta(n, 2, 2)
  
  # Test with and without starting values
  res1 <- SeqExpMatch:::fast_beta_regression_mle(y, X)
  res2 <- SeqExpMatch:::fast_beta_regression_mle(y, X, start_beta = c(0, 0), start_phi = 5)
  expect_equal(as.numeric(res1$coefficients), as.numeric(res2$coefficients), tolerance = 1e-4)
  
  # Test with and without standard error computation
  res3 <- SeqExpMatch:::fast_beta_regression_mle(y, X, compute_std_errs = FALSE)
  res4 <- SeqExpMatch:::fast_beta_regression_mle(y, X, compute_std_errs = TRUE)
  expect_null(res3$std_errs)
  expect_false(is.null(res4$std_errs))
  expect_equal(length(res4$std_errs), 3) # beta0, beta1, phi
})

test_that("fast_weibull_regression robustness", {
  set.seed(1)
  n <- 50
  X <- matrix(rnorm(n), n, 1)
  y <- rexp(n)
  dead <- rep(1L, n)
  
  # Check it runs without error on small sample
  expect_error(SeqExpMatch:::fast_weibull_regression(y, dead, X), NA)
})

test_that("robust_survreg_with_surv_object fallback and optimization", {
  set.seed(1)
  n <- 100
  X <- matrix(rnorm(n), n, 1)
  y <- rexp(n)
  dead <- rbinom(n, 1, 0.9)
  surv_obj <- Surv(y, dead)
  
  # Test normal case (should use optimization)
  mod1 <- robust_survreg_with_surv_object(surv_obj, X, dist = "weibull")
  expect_s3_class(mod1, "survreg")
  
  # Test with non-weibull (should use original logic)
  mod2 <- robust_survreg_with_surv_object(surv_obj, X, dist = "lognormal")
  expect_s3_class(mod2, "survreg")
  
  # Test with collinearity (should trigger feature elimination)
  X_coll <- cbind(X, X) # Perfectly collinear
  # robust_survreg_with_surv_object uses cor() and eliminates columns
  # If it works, it should drop one X column.
  # Intercept + 1 feature = 2 coefficients.
  mod3 <- robust_survreg_with_surv_object(surv_obj, X_coll, dist = "weibull")
  expect_s3_class(mod3, "survreg")
  expect_equal(length(coef(mod3)), 2)
})

test_that("create_design_matrix consistency", {
  # We test this by initializing an inference object and checking the private method
  set.seed(1)
  n <- 20
  seq_des <- SeqDesignCRD$new(n = n, response_type = "continuous")
  X <- data.frame(x1 = rnorm(n))
  for(i in 1:n) seq_des$add_subject_to_experiment_and_assign(X[i, , drop=FALSE])
  seq_des$add_all_subject_responses(rnorm(n))
  
  inf_obj <- SeqDesignInferenceAllSimpleMeanDiff$new(seq_des)
  
  # Access private method
  dm <- inf_obj$.__enclos_env__$private$create_design_matrix()
  
  expect_equal(nrow(dm), n)
  expect_equal(ncol(dm), 1 + 1 + 1) # Intercept + w + x1
  expect_equal(as.numeric(dm[, 1]), rep(1, n))
  expect_equal(as.numeric(dm[, 2]), as.numeric(seq_des$get_w()))
})

test_that("fast_beta_regression_with_var consistency", {
  set.seed(1)
  n <- 100
  X <- matrix(rnorm(n), n, 1)
  y <- rbeta(n, 5, 5) # more central
  Xmm <- cbind(1, X)
  
  res <- SeqExpMatch:::fast_beta_regression_with_var(Xmm, y)
  expect_named(res, c("b", "phi", "ssq_b_2"))
  expect_true(res$ssq_b_2 > 0)
})
