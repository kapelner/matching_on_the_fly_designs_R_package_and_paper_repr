library(testthat)
library(EDI)

test_that("SimulationFramework handles multiple cells and summarizes correctly", {
	# Test grid of parameters
	results_file <- tempfile(fileext = ".csv")
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(FixedDesignBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		n = c(10, 20),
		p = c(1, 2),
		betaT = c(0, 0.5),
		Nrep = 2,
		inference_types_and_params = list(asymp_pval = list(delta = 0)),
		results_filename = results_file,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	sim$run()
	sm <- sim$summarize()
	
	# Grid: 2(n) * 2(p) * 2(betaT) = 8 cells
	expect_equal(nrow(sm), 8)
	expect_true(all(c("n", "p", "betaT", "MSE", "power") %in% names(sm)))
	
	# Power should be close to alpha (0.05) for betaT = 0, but Nrep=2 is too small for statistical check.
	# Just check it ran and produced finite values.
	expect_true(all(is.finite(sm$MSE)))
})

test_that("SimulationFramework continue logic works for both csv and csv.bz2", {
	for (ext in c(".csv", ".csv.bz2")) {
		results_file <- tempfile(fileext = ext)
		
		# First run: 1 rep
		sim1 <- SimulationFramework$new(
			response_type = "continuous",
			n = 10L, Nrep = 1L,
			results_filename = results_file,
			verbose = FALSE,
			continue_from_last_result_row = FALSE
		)
		sim1$run()
		expect_equal(nrow(sim1$get_results()), 1)
		
		# Second run: continue, total Nrep = 3
		sim2 <- SimulationFramework$new(
			response_type = "continuous",
			n = 10L, Nrep = 3L,
			results_filename = results_file,
			verbose = FALSE,
			continue_from_last_result_row = TRUE
		)
		sim2$run()
		res <- sim2$get_results()
		expect_equal(nrow(res), 3)
		expect_equal(sort(res$rep), 1:3)
	}
})

test_that("SimulationFramework handles seed for reproducibility", {
	results_file1 <- tempfile(fileext = ".csv")
	sim1 <- SimulationFramework$new(
		response_type = "continuous",
		n = 10L, Nrep = 2L,
		seed = 12345,
		results_filename = results_file1,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	sim1$run()
	res1 <- sim1$get_results()
	
	results_file2 <- tempfile(fileext = ".csv")
	sim2 <- SimulationFramework$new(
		response_type = "continuous",
		n = 10L, Nrep = 2L,
		seed = 12345,
		results_filename = results_file2,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	sim2$run()
	res2 <- sim2$get_results()
	
	expect_equal(res1$estimate, res2$estimate)
})

test_that("SimulationFramework validates response_type", {
	expect_error(
		SimulationFramework$new(response_type = "invalid"),
		"must be one of"
	)
})

test_that("SimulationFramework handles factor covariate in design initialization", {
	# Some designs might expect factors. SimulationFramework auto-injects some defaults.
	# Let's test if it handles a custom categorical covariate.
	
	set.seed(55)
	n <- 20
	X <- data.frame(
		x1 = rnorm(n),
		cat1 = factor(sample(letters[1:3], n, replace = TRUE))
	)
	
	sim <- SimulationFramework$new(
		response_type = "continuous",
		X_mat = as.matrix(model.matrix(~ x1 + cat1, data = X)[, -1]),
		Nrep = 1,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	# If X_mat is supplied, it must be numeric. model.matrix handled it.
	sim$run()
	expect_equal(nrow(sim$get_results()), 1)
})

test_that("SimulationFramework summarize handles all-NA results", {
	sim <- SimulationFramework$new(
		response_type = "continuous",
		Nrep = 1,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	priv <- sim$.__enclos_env__$private
	priv$has_run <- TRUE
	# Inject a result with NA estimate
	priv$raw_results <- list(data.table::data.table(
		response_type = "continuous",
		cond_exp_func_model = "linear",
		n = 100, p = 5, betaT = 1, rep = 1,
		design = "FixedDesignBernoulli",
		inference = "InferenceAllSimpleMeanDiff",
		inference_type = "asymp_pval",
		estimate = NA_real_,
		pval = NA_real_,
		true_estimand = 1.0
	))
	priv$results_idx <- 1L
	
	sm <- sim$summarize()
	expect_true(is.na(sm$MSE))
	expect_equal(sm$n_est, 0L)
})

test_that("SimulationFramework respects clamping parameters", {
	# Force an extreme signal that would normally result in 0/1 probabilities
	sim <- SimulationFramework$new(
		response_type = "incidence",
		n = 10,
		betaT = 100, # huge effect
		incidence_clamp = 0.01,
		Nrep = 1,
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)
	
	priv <- sim$.__enclos_env__$private
	# We want to check if the true_estimand calculation uses the clamp.
	# It's calculated in .run_single_replication_in_worker and returned in raw results.
	sim$run()
	res <- sim$get_results()
	
	# With betaT = 100, p_t should be clamped to 0.99.
	# If eta is large, p_c might be clamped to 0.99 too, making diff 0.
	# If eta is 0, p_c is 0.5. diff = 0.99 - 0.5 = 0.49.
	# If we didn't clamp, p_t would be 1.0.
	
	# The ATE calculation in .run_single_replication_in_worker uses the clamp.
	expect_true(all(res$true_estimand <= 0.99))
})
