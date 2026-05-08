test_that("SimulationFramework accepts merged design classes and params", {
	set.seed(20260429)
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(
			FixedDesignBernoulli,
			FixedDesignBlocking = list(B_target = 2L)
		),
		inference_classes_and_params = list(
			InferenceAllSimpleMeanDiff = list(max_resample_attempts = 10L)
		),
		n = 8L,
		p = 2L,
		Nrep = 1L,
		betaT = 0,
		results_filename = tempfile(fileext = ".csv"),
		continue_from_last_result_row = FALSE,
		verbose = FALSE,
		inference_types_and_params = list(asymp_pval = list(delta = 0)),
		turn_off_asserts_for_speed = FALSE
	)

	sim$run()
	raw <- sim$get_results()
	sm <- sim$summarize()

	expect_true("ci_lo" %in% names(raw))
	expect_true("ci_hi" %in% names(raw))
	expect_false("ci_a" %in% names(raw))
	expect_false("ci_b" %in% names(raw))
	expect_true("asymp_pval" %in% raw$inference_type)
	expect_true("asymp_pval" %in% sm$inference_type)
	expect_equal(sort(sm$design), sort(c(
		"FixedDesignBernoulli",
		"FixedDesignBlocking (B_target=2L)"
	)))
	expect_true("design_params" %in% names(sm))
	expect_true("B_target=2L" %in% sm$design_params)
	expect_true("max_resample_attempts=10L" %in% sm$inference_params)
	expect_true("asymp_pval(delta=0)" %in% sm$inference_type_params)
})

test_that("SimulationFramework preserves FixedDesignBlocking NULL strata_cols semantics", {
	set.seed(20260506)
	sim <- SimulationFramework$new(
		response_type = "incidence",
		design_classes_and_params = list(
			FixedDesignBlocking = list(B_target = 4L, exact_num_blocks = TRUE)
		),
		inference_classes_and_params = list(InferenceIncidCMH),
		inference_types_and_params = list(asymp_pval = list()),
		n = 8L,
		p = 1L,
		Nrep = 1L,
		betaT = 0,
		results_filename = tempfile(fileext = ".csv"),
		continue_from_last_result_row = FALSE,
		verbose = FALSE,
		turn_off_asserts_for_speed = FALSE
	)

	priv <- sim$.__enclos_env__$private
	priv$current_n <- 8L
	priv$current_p <- 1L
	priv$current_betaT <- 0
	priv$current_cond_exp_func_model <- "linear"
	priv$current_response_type <- "incidence"
	dat <- list(
		X = data.frame(x = 1:8),
		y_linear_model = as.numeric(scale(1:8))
	)
	combos <- priv$.build_valid_combos_for_current_cell(dat)
	designs <- unique(vapply(combos, `[[`, "", "design"))
	expect_true("FixedDesignBlocking (B_target=4L, exact_num_blocks=TRUE)" %in% designs)
})

test_that("SimulationFramework defaults to csv.bz2 results and rejects unsupported extensions", {
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(FixedDesignBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		verbose = FALSE,
		continue_from_last_result_row = FALSE
	)

	expect_equal(sim$.__enclos_env__$private$results_filename, "simulation_framework_results.csv.bz2")

	expect_error(
		SimulationFramework$new(
			response_type = "continuous",
			design_classes_and_params = list(FixedDesignBernoulli),
			inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
			results_filename = tempfile(fileext = ".txt"),
			verbose = FALSE,
			continue_from_last_result_row = FALSE
		),
		"results_filename must end in either '.csv' or '.csv.bz2'"
	)
})

test_that("SimulationFramework validates merged inference constructor params", {
	expect_error(
		SimulationFramework$new(
			response_type = "continuous",
			design_classes_and_params = list(FixedDesignBernoulli),
			inference_classes_and_params = list(
				InferenceAllSimpleMeanDiff = list(not_a_constructor_arg = TRUE)
			),
			results_filename = tempfile(fileext = ".csv"),
			continue_from_last_result_row = FALSE,
			verbose = FALSE
		),
		"not accepted by initialize"
	)
})

test_that("SimulationFramework validates inference type params against function arguments", {
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(FixedDesignBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(asymp_pval = list(not_an_arg = TRUE)),
		n = 8L,
		p = 2L,
		Nrep = 1L,
		betaT = 0,
		results_filename = tempfile(fileext = ".csv"),
		continue_from_last_result_row = FALSE,
		verbose = FALSE
	)

	expect_error(sim$run(), "not accepted by compute_asymp_two_sided_pval")
})

test_that("SimulationFramework true mean-difference estimands match DGP scale", {
	get_true_estimand <- function(response_type, y_linear_model, betaT, ...) {
		sim <- SimulationFramework$new(
			response_type = response_type,
			design_classes_and_params = list(FixedDesignBernoulli),
			inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
			inference_types_and_params = list(asymp_pval = list(delta = 0)),
			n = length(y_linear_model),
			p = 1L,
			Nrep = 1L,
			betaT = betaT,
			X_mat = matrix(y_linear_model, ncol = 1L),
			results_filename = tempfile(fileext = ".csv"),
			continue_from_last_result_row = FALSE,
			verbose = FALSE,
			turn_off_asserts_for_speed = FALSE,
			...
		)
		priv <- sim$.__enclos_env__$private
		priv$compute_true_mean_diff_ate(y_linear_model)
	}

	eta <- c(-1, -0.25, 0.5, 1.25)
	betaT <- 0.7

	expect_equal(
		get_true_estimand("incidence", eta, betaT, incidence_clamp = 1e-6),
		mean(pmin(1 - 1e-6, pmax(1e-6, plogis(eta + betaT))) -
			pmin(1 - 1e-6, pmax(1e-6, plogis(eta)))),
		tolerance = 1e-12
	)
	expect_equal(
		get_true_estimand("proportion", eta, betaT, proportion_clamp = 1e-6),
		mean(pmin(1 - 1e-6, pmax(1e-6, plogis(eta + betaT))) -
			pmin(1 - 1e-6, pmax(1e-6, plogis(eta)))),
		tolerance = 1e-12
	)
	expect_equal(
		get_true_estimand("count", eta, betaT, count_clamp = 1e-6),
		mean(pmax(1e-6, exp(eta + betaT)) - pmax(1e-6, exp(eta))),
		tolerance = 1e-12
	)
	expect_equal(
		get_true_estimand("ordinal", eta, betaT, sd_noise = 0, n_ordinal_levels = 4L),
		mean(pmin(4L, pmax(1L, round(eta + betaT))) - pmin(4L, pmax(1L, round(eta)))),
		tolerance = 1e-12
	)
})

test_that("SimulationFramework accepts vector grids for n, p, betaT, and cond_exp_func_model", {
	sim <- suppressWarnings(SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(FixedDesignBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(asymp_pval = list(delta = 0)),
		n = c(6L, 8L),
		p = c(2L, 5L),
		Nrep = 1L,
		betaT = c(0, 1),
		cond_exp_func_model = c("linear", "nonlinear"),
		results_filename = tempfile(fileext = ".csv"),
		continue_from_last_result_row = FALSE,
		verbose = FALSE,
		turn_off_asserts_for_speed = FALSE
	))

	sim$run()
	raw <- sim$get_results()
	sm <- sim$summarize()

	expect_true(all(c("cond_exp_func_model", "n", "p", "betaT") %in% names(raw)))
	expect_true(all(c("cond_exp_func_model", "n", "p", "betaT") %in% names(sm)))
	expect_equal(data.table::uniqueN(sm[, .(cond_exp_func_model, n, p, betaT)]), 12L)
	expect_false(any(sm$cond_exp_func_model == "nonlinear" & sm$p < 5L))
})

test_that("SimulationFramework can write and reload csv.bz2 results", {
	results_file <- tempfile(fileext = ".csv.bz2")
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(FixedDesignBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(asymp_pval = list(delta = 0)),
		n = 8L,
		p = 2L,
		Nrep = 1L,
		betaT = 0,
		results_filename = results_file,
		continue_from_last_result_row = TRUE,
		verbose = FALSE,
		turn_off_asserts_for_speed = FALSE
	)
	sim$run()
	raw_first <- sim$get_results()

	expect_true(file.exists(results_file))
	expect_false(file.exists(file.path(
		dirname(results_file),
		paste0(".", sub("\\.csv\\.bz2$", "", basename(results_file)), "__staging.csv")
	)))

	sim_resume <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(FixedDesignBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(asymp_pval = list(delta = 0)),
		n = 8L,
		p = 2L,
		Nrep = 1L,
		betaT = 0,
		results_filename = results_file,
		continue_from_last_result_row = TRUE,
		verbose = FALSE,
		turn_off_asserts_for_speed = FALSE
	)
	sim_resume$run()
	raw_second <- sim_resume$get_results()

	expect_equal(nrow(raw_second), nrow(raw_first))
	expect_equal(raw_second$estimate, raw_first$estimate)
})

test_that("SimulationFramework summarize preserves raw metric precision", {
	sim <- SimulationFramework$new(
		response_type = "continuous",
		design_classes_and_params = list(FixedDesignBernoulli),
		inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
		inference_types_and_params = list(
			asymp_ci = list(),
			asymp_pval = list(delta = 0)
		),
		n = 4L,
		p = 1L,
		Nrep = 1L,
		betaT = 0,
		results_filename = tempfile(fileext = ".csv"),
		continue_from_last_result_row = FALSE,
		verbose = FALSE
	)
	priv <- sim$.__enclos_env__$private
	priv$has_run <- TRUE
	priv$valid_combos <- list(list(
		response_type = "continuous",
		cond_exp_func_model = "linear",
		n = 4L,
		p = 1L,
		betaT = 0,
		design = "FixedDesignBernoulli",
		inference = "InferenceAllSimpleMeanDiff",
		inference_type = "asymp_pval"
	))
	priv$raw_results <- list(
		list(
			response_type = "continuous",
			rep = 1L,
			cond_exp_func_model = "linear",
			n = 4L,
			p = 1L,
			betaT = 0,
			design = "FixedDesignBernoulli",
			inference = "InferenceAllSimpleMeanDiff",
			inference_type = "asymp_pval",
			estimate = 1.234567,
			ci_lo = -0.5,
			ci_hi = 1.5,
			pval = 0.5,
			true_estimand = 0.5
		),
		list(
			response_type = "continuous",
			rep = 2L,
			cond_exp_func_model = "linear",
			n = 4L,
			p = 1L,
			betaT = 0,
			design = "FixedDesignBernoulli",
			inference = "InferenceAllSimpleMeanDiff",
			inference_type = "asymp_pval",
			estimate = 0.111111,
			ci_lo = 0.6,
			ci_hi = 1.0,
			pval = 0.987654,
			true_estimand = 0.5
		),
		list(
			response_type = "continuous",
			rep = 3L,
			cond_exp_func_model = "linear",
			n = 4L,
			p = 1L,
			betaT = 0,
			design = "FixedDesignBernoulli",
			inference = "InferenceAllSimpleMeanDiff",
			inference_type = "asymp_pval",
			estimate = 0.222222,
			ci_lo = 0.0,
			ci_hi = 1.2,
			pval = 0.012345,
			true_estimand = 0.5
		)

	)
	priv$results_idx <- 3L

	sm <- sim$summarize()
	expect_equal(sm$MSE, mean(c((1.234567 - 0.5)^2, (0.111111 - 0.5)^2, (0.222222 - 0.5)^2)), tolerance = 1e-15)
	expect_equal(sm$power, 1 / 3, tolerance = 1e-15)
	expect_equal(sm$ci_length, mean(c(2.0, 0.4, 1.2)), tolerance = 1e-15)
	expect_equal(sm$coverage, 2 / 3, tolerance = 1e-15)
	expect_false(identical(sm$coverage, round(sm$coverage, 3)))
	expect_false(identical(sm$MSE, round(sm$MSE, 5)))
})
