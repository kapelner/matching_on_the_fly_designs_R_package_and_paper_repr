test_that("Inference works for continuous", {
  n <- 10
  des <- SeqDesignCRD$new(n = n, response_type = "continuous", verbose = FALSE)
  set.seed(1)
  for (i in 1:n) {
    des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
  }
  des$add_all_subject_responses(rnorm(n))
  
  # Simple Mean Diff
  inf <- SeqDesignInferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
  est <- inf$compute_treatment_estimate()
  expect_true(is.numeric(est))
  
  # OLS
  inf_ols <- SeqDesignInferenceContinMultOLS$new(des, verbose = FALSE)
  est_ols <- inf_ols$compute_treatment_estimate()
  expect_true(is.numeric(est_ols))
})

test_that("Inference works for incidence", {
  n <- 30
  des <- SeqDesignCRD$new(n = n, response_type = "incidence", verbose = FALSE)
  set.seed(1)
  for (i in 1:n) {
    des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
  }
  des$add_all_subject_responses(rbinom(n, 1, 0.5))
  
  inf <- SeqDesignInferenceIncidUnivLogRegr$new(des, verbose = FALSE)
  est <- inf$compute_treatment_estimate()
  expect_true(is.numeric(est))
})

test_that("Inference works for count", {
  n <- 20
  des <- SeqDesignCRD$new(n = n, response_type = "count", verbose = FALSE)
  set.seed(1)
  for (i in 1:n) {
    des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
  }
  des$add_all_subject_responses(rpois(n, 5))
  
  inf <- SeqDesignInferenceCountUnivNegBinRegr$new(des, verbose = FALSE)
  est <- inf$compute_treatment_estimate()
  expect_true(is.numeric(est))
})

test_that("Inference works for proportion", {
  n <- 10
  des <- SeqDesignCRD$new(n = n, response_type = "proportion", verbose = FALSE)
  set.seed(1)
  for (i in 1:n) {
    des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
  }
  des$add_all_subject_responses(runif(n))
  
  inf <- SeqDesignInferencePropUniBetaRegr$new(des, verbose = FALSE)
  est <- inf$compute_treatment_estimate()
  expect_true(is.numeric(est))
})

test_that("Inference works for survival", {
  n <- 20
  des <- SeqDesignCRD$new(n = n, response_type = "survival")
  set.seed(1)
  for (i in 1:n) {
    des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
  }
  des$add_all_subject_responses(rexp(n), dead = rbinom(n, 1, 0.8))
  
  # KM Diff
  inf <- SeqDesignInferenceSurvivalKMDiff$new(des, verbose = FALSE)
  est <- inf$compute_treatment_estimate()
  expect_true(is.numeric(est))
  
  # Cox PH
  inf_cox <- SeqDesignInferenceSurvivalUniCoxPHRegr$new(des, verbose = FALSE)
  est_cox <- inf_cox$compute_treatment_estimate()
  expect_true(is.numeric(est_cox))
})
