.libPaths(c("local_R_lib", .libPaths()))
library(SeqExpMatch)
library(data.table)
library(testthat)

test_that("SeqDesignSPBR Stratification and Randomization", {
  n = 20
  strata_cols = c("gender", "age_cat")
  
  # Create some test data
  X = data.table(
    gender = rep(c("M", "F"), each = 10),
    age_cat = rep(rep(c("Young", "Old"), each = 5), 2)
  )
  
  seq_des = SeqDesignSPBR$new(strata_cols = strata_cols, block_size = 2, n = n)
  
  cat("Starting subjects addition...\n")
  assignments = numeric(n)
  for (i in 1:n) {
    assignments[i] = seq_des$add_subject_to_experiment_and_assign(X[i, ])
  }
  cat("Subjects addition finished.\n")
  
  check_stratum = function(idxs) {
    strat_w = assignments[idxs]
    # Blocks of size 2
    for (j in seq(1, length(idxs)-1, by=2)) {
      expect_equal(sum(strat_w[j:(j+1)]), 1)
    }
  }
  
  check_stratum(1:5) # Stratum 1
  check_stratum(6:10) # Stratum 2
  check_stratum(11:15) # Stratum 3
  check_stratum(16:20) # Stratum 4
  
  cat("Testing redraw...\n")
  old_w = copy(seq_des$get_w())
  seq_des$.__enclos_env__$private$redraw_w_according_to_design()
  new_w = seq_des$get_w()
  cat("Redraw finished.\n")
  
  # Verify redraw_w_according_to_design still maintains SPBR properties
  strat_keys = vapply(1:n, function(i) {
    seq_des$.__enclos_env__$private$get_strata_key(X[i, ])
  }, character(1))
  
  # Check balance in new_w per stratum
  unique_keys = unique(strat_keys)
  for (key in unique_keys) {
    idxs = which(strat_keys == key)
    strat_w = new_w[idxs]
    for (j in seq(1, length(idxs)-1, by=2)) {
      expect_equal(sum(strat_w[j:(j+1)]), 1)
    }
  }
  
  # Test Stratified Bootstrap Indices
  cat("Testing stratified bootstrap...\n")
  # Use an inference class to trigger bootstrap
  seq_des$add_all_subject_responses(rnorm(n))
  inf_obj = SeqDesignInferenceAllSimpleMeanDiff$new(seq_des)
  
  # Capture the indices by overriding or just checking the property
  # We'll just run it to make sure it doesn't crash
  bs_dist = inf_obj$approximate_bootstrap_distribution_beta_hat_T(B = 10)
  expect_equal(length(bs_dist), 10)
  expect_true(all(!is.na(bs_dist)))
  
  cat("SPBR tests passed!\n")
})
