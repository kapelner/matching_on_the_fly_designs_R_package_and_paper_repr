test_that("SeqDesignCRD works", {
  # Fixed sample
  des <- SeqDesignCRD$new(n = 10, verbose = FALSE)
  expect_true(des$is_fixed_sample_size())
  expect_equal(des$get_n(), 10)
  
  x_new <- data.frame(x1 = rnorm(1))
  w <- des$add_subject_to_experiment_and_assign(x_new)
  expect_true(w %in% c(0, 1))
  expect_equal(des$get_t(), 1)
  
  des$add_subject_response(1, 0.5)
  expect_equal(des$get_y()[1], 0.5)
  
  # Not fixed sample
  des <- SeqDesignCRD$new(n = NULL, verbose = FALSE)
  expect_false(des$is_fixed_sample_size())
  des$add_subject_to_experiment_and_assign(x_new)
  expect_equal(des$get_t(), 1)
})

test_that("SeqDesignEfron works", {
  des <- SeqDesignEfron$new(n = 10, prob_T = 0.5, verbose = FALSE)
  expect_true(des$is_fixed_sample_size())
  
  x_new <- data.frame(x1 = rnorm(1))
  w <- des$add_subject_to_experiment_and_assign(x_new)
  expect_true(w %in% c(0, 1))
})

test_that("SeqDesignAtkinson works", {
  des <- SeqDesignAtkinson$new(n = 10, verbose = FALSE)
  x_new <- data.frame(x1 = rnorm(1), x2 = rnorm(1))
  w <- des$add_subject_to_experiment_and_assign(x_new)
  expect_true(w %in% c(0, 1))
})

test_that("SeqDesignKK14 works", {
  des <- SeqDesignKK14$new(n = 10, verbose = FALSE)
  x_new <- data.frame(x1 = rnorm(1), x2 = rnorm(1))
  w <- des$add_subject_to_experiment_and_assign(x_new)
  expect_true(w %in% c(0, 1))
})

test_that("SeqDesignKK21 works", {
  des <- SeqDesignKK21$new(n = 10, verbose = FALSE)
  x_new <- data.frame(x1 = rnorm(1), x2 = rnorm(1))
  w <- des$add_subject_to_experiment_and_assign(x_new)
  expect_true(w %in% c(0, 1))
})

test_that("SeqDesignKK21stepwise works", {
  des <- SeqDesignKK21stepwise$new(n = 10, verbose = FALSE)
  x_new <- data.frame(x1 = rnorm(1), x2 = rnorm(1))
  w <- des$add_subject_to_experiment_and_assign(x_new)
  expect_true(w %in% c(0, 1))
})

test_that("SeqDesigniBCRD works", {
  des <- SeqDesigniBCRD$new(n = 10, verbose = FALSE)
  x_new <- data.frame(x1 = rnorm(1))
  w <- des$add_subject_to_experiment_and_assign(x_new)
  expect_true(w %in% c(0, 1))
})

test_that("Response types work", {
  types <- c("continuous", "incidence", "proportion", "count", "survival")
  for (rt in types) {
    des <- SeqDesignCRD$new(n = 10, response_type = rt, verbose = FALSE)
    expect_equal(des$get_response_type(), rt)
    
    x_new <- data.frame(x1 = rnorm(1))
    des$add_subject_to_experiment_and_assign(x_new)
    
    val <- switch(rt,
                  continuous = 1.5,
                  incidence = 1,
                  proportion = 0.5,
                  count = 5,
                  survival = 10)
    
    if (rt == "survival") {
      des$add_subject_response(1, val, dead = 1)
      expect_equal(des$get_dead()[1], 1)
    } else {
      des$add_subject_response(1, val)
    }
    expect_equal(des$get_y()[1], val)
  }
})
