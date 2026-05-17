library(testthat)
library(EDI)

context("Smart Default FALSE (Legacy) Paths")

test_that("Logistic regression works with smart_default = FALSE", {
    n = 100
    p = 5
    set.seed(42)
    X = matrix(rnorm(n * p), n, p)
    X[, 1] = 1
    w = rbinom(n, 1, 0.5)
    beta = rnorm(p) * 0.5
    y = rbinom(n, 1, plogis(X %*% beta))
    
    des = DesignFixedBernoulli$new(n = n, response_type = "incidence")
    des$add_all_subjects_to_experiment(as.data.frame(X[, -1]))
    des$overwrite_all_subject_assignments(w)
    des$add_all_subject_responses(y)
    
    # Cold start with legacy (all-zero) start
    inf = InferenceIncidLogRegr$new(des, smart_cold_start_default = FALSE)
    est = inf$compute_estimate()
    expect_true(is.finite(est))
    
    # Verify it converges to same solution as smart
    inf_smart = InferenceIncidLogRegr$new(des, smart_cold_start_default = TRUE)
    est_smart = inf_smart$compute_estimate()
    expect_equal(est, est_smart, tolerance = 1e-6)
})

test_that("Poisson regression works with smart_default = FALSE", {
    n = 100
    p = 5
    set.seed(42)
    X = matrix(rnorm(n * p), n, p)
    X[, 1] = 1
    w = rbinom(n, 1, 0.5)
    beta = rnorm(p) * 0.2
    y = rpois(n, exp(X %*% beta))
    
    des = DesignFixedBernoulli$new(n = n, response_type = "count")
    des$add_all_subjects_to_experiment(as.data.frame(X[, -1]))
    des$overwrite_all_subject_assignments(w)
    des$add_all_subject_responses(y)
    
    inf = InferenceCountPoisson$new(des, smart_cold_start_default = FALSE)
    est = inf$compute_estimate()
    expect_true(is.finite(est))
    
    inf_smart = InferenceCountPoisson$new(des, smart_cold_start_default = TRUE)
    est_smart = inf_smart$compute_estimate()
    expect_equal(est, est_smart, tolerance = 1e-6)
})

test_that("Negative Binomial works with smart_default = FALSE", {
    n = 100
    p = 5
    set.seed(42)
    X = matrix(rnorm(n * p), n, p)
    X[, 1] = 1
    w = rbinom(n, 1, 0.5)
    beta = rnorm(p) * 0.2
    y = rnbinom(n, size = 2, mu = exp(X %*% beta))
    
    des = DesignFixedBernoulli$new(n = n, response_type = "count")
    des$add_all_subjects_to_experiment(as.data.frame(X[, -1]))
    des$overwrite_all_subject_assignments(w)
    des$add_all_subject_responses(y)
    
    inf = InferenceCountNegBin$new(des, smart_cold_start_default = FALSE)
    est = inf$compute_estimate()
    expect_true(is.finite(est))
    
    inf_smart = InferenceCountNegBin$new(des, smart_cold_start_default = TRUE)
    est_smart = inf_smart$compute_estimate()
    expect_equal(est, est_smart, tolerance = 1e-6)
})

test_that("Weibull AFT works with smart_default = FALSE", {
    n = 100
    p = 5
    set.seed(42)
    X = matrix(rnorm(n * p), n, p)
    X[, 1] = 1
    w = rbinom(n, 1, 0.5)
    beta = rnorm(p) * 0.2
    y = rexp(n, exp(X %*% beta))
    dead = rbinom(n, 1, 0.8)
    
    des = DesignFixedBernoulli$new(n = n, response_type = "survival")
    des$add_all_subjects_to_experiment(as.data.frame(X[, -1]))
    des$overwrite_all_subject_assignments(w)
    des$add_all_subject_responses(y, dead = dead)
    
    inf = InferenceSurvivalWeibullRegr$new(des, smart_cold_start_default = FALSE)
    est = inf$compute_estimate()
    expect_true(is.finite(est))
    
    inf_smart = InferenceSurvivalWeibullRegr$new(des, smart_cold_start_default = TRUE)
    est_smart = inf_smart$compute_estimate()
    expect_equal(est, est_smart, tolerance = 1e-6)
})
