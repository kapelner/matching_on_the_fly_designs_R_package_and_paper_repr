test_that("Different designs work with continuous response", {
  n <- 20
  designs <- list(
    CRD = SeqDesignCRD$new(n = n),
    Efron = SeqDesignEfron$new(n = n),
    Atkinson = SeqDesignAtkinson$new(n = n),
    iBCRD = SeqDesigniBCRD$new(n = n),
    KK14 = SeqDesignKK14$new(n = n),
    KK21 = SeqDesignKK21$new(n = n),
    KK21stepwise = SeqDesignKK21stepwise$new(n = n)
  )
  
  for (name in names(designs)) {
    des <- designs[[name]]
    set.seed(1)
    for (i in 1:n) {
      des$add_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
    }
    des$add_all_subject_responses(rnorm(n))
    
    # Check if we can run a simple inference
    inf <- SeqDesignInferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
    expect_true(is.numeric(inf$compute_treatment_estimate()), info = paste("Design:", name))
  }
})

test_that("KK14 with Morrison works", {
  n <- 20
  p <- 2
  des <- SeqDesignKK14$new(n = n, morrison = TRUE, p = p)
  set.seed(1)
  for (i in 1:n) {
    des$add_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
  }
  des$add_all_subject_responses(rnorm(n))
  
  inf <- SeqDesignInferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
  expect_true(is.numeric(inf$compute_treatment_estimate()))
})

test_that("Multivariate inference works", {
  n <- 40
  des <- SeqDesignCRD$new(n = n, response_type = "continuous")
  set.seed(1)
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  for (i in 1:n) {
    des$add_subject_to_experiment_and_assign(X[i, , drop = FALSE])
  }
  des$add_all_subject_responses(rnorm(n))
  
  # Continuous MultOLS
  inf_cont <- SeqDesignInferenceContinMultOLS$new(des, verbose = FALSE)
  expect_true(is.numeric(inf_cont$compute_treatment_estimate()))
  
  # Incidence MultiLogRegr
  des_inc <- SeqDesignCRD$new(n = n, response_type = "incidence")
  for (i in 1:n) {
    des_inc$add_subject_to_experiment_and_assign(X[i, , drop = FALSE])
  }
  des_inc$add_all_subject_responses(rbinom(n, 1, 0.5))
  inf_inc <- SeqDesignInferenceIncidMultiLogRegr$new(des_inc, verbose = FALSE)
  expect_true(is.numeric(inf_inc$compute_treatment_estimate()))
  
  # Count MultiNegBinRegr
  des_count <- SeqDesignCRD$new(n = n, response_type = "count")
  for (i in 1:n) {
    des_count$add_subject_to_experiment_and_assign(X[i, , drop = FALSE])
  }
  des_count$add_all_subject_responses(rpois(n, 5))
  inf_count <- SeqDesignInferenceCountMultiNegBinRegr$new(des_count, verbose = FALSE)
  expect_true(is.numeric(inf_count$compute_treatment_estimate()))
})

test_that("KK designs with KK inference work", {
  n <- 30
  designs <- list(
    KK14 = SeqDesignKK14$new(n = n),
    KK21 = SeqDesignKK21$new(n = n),
    KK21stepwise = SeqDesignKK21stepwise$new(n = n)
  )
  
  inf_classes <- list(
    KK14 = SeqDesignInferenceBaiAdjustedTKK14,
    KK21 = SeqDesignInferenceBaiAdjustedTKK21,
    KK21stepwise = SeqDesignInferenceBaiAdjustedTKK21
  )
  
  for (name in names(designs)) {
    des <- designs[[name]]
    set.seed(1)
    for (i in 1:n) {
      des$add_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
    }
    des$add_all_subject_responses(rnorm(n))
    
    inf <- inf_classes[[name]]$new(des, verbose = FALSE)
    expect_true(is.numeric(inf$compute_treatment_estimate()), info = paste("KK Design/Inf:", name))
  }
})
