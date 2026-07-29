library(EDI)

build_resampling_smoke_inference = function(n = 20L, seed = 123L){
	set.seed(seed)
	des = DesignFixedBernoulli$new(n = n, response_type = "continuous", seed = seed)
	des$add_all_subjects_to_experiment(data.frame(x = rnorm(n)))
	des$assign_w_to_all_subjects()
	w = des$get_w()
	for (t in seq_len(n)) {
		des$add_one_subject_response(t, w[t] + stats::rnorm(1, sd = 0.05))
	}
	inf = InferenceAllSimpleMeanDiff$new(des)
	inf$num_cores = 1L
	inf$.__enclos_env__$private$seed = seed
	inf
}

test_that("m-out-of-n bootstrap methods are available and return scalar inference", {
	inf = build_resampling_smoke_inference()
	expect_true(is.function(inf$approximate_m_out_of_n_bootstrap_distribution_beta_hat_T))
	expect_true(is.function(inf$compute_m_out_of_n_bootstrap_two_sided_pval))
	expect_true(is.function(inf$compute_m_out_of_n_bootstrap_confidence_interval))
	expect_true(is.function(inf$select_optimal_m_out_of_n_bootstrap))

	d = inf$approximate_m_out_of_n_bootstrap_distribution_beta_hat_T(B = 11, m = 8, show_progress = FALSE)
	expect_length(d, 11)
	expect_true(any(is.finite(d)))

	ci = inf$compute_m_out_of_n_bootstrap_confidence_interval(B = 21, m = 8, show_progress = FALSE)
	expect_length(ci, 2)
	expect_true(all(is.finite(ci)))
	expect_lte(ci[1], ci[2])

	p = inf$compute_m_out_of_n_bootstrap_two_sided_pval(B = 21, m = 8, show_progress = FALSE)
	expect_true(is.finite(p))
	expect_gte(p, 0)
	expect_lte(p, 1)
})

test_that("NULL is the only automatic m/b size sentinel", {
	inf = build_resampling_smoke_inference()
	expected_size = as.integer(floor(20 ^ 0.7))

	expect_null(formals(inf$approximate_m_out_of_n_bootstrap_distribution_beta_hat_T)$m)
	expect_null(formals(inf$compute_m_out_of_n_bootstrap_two_sided_pval)$m)
	expect_null(formals(inf$compute_m_out_of_n_bootstrap_confidence_interval)$m)
	expect_null(formals(inf$approximate_subsampling_distribution_beta_hat_T)$b)
	expect_null(formals(inf$compute_subsampling_two_sided_pval)$b)
	expect_null(formals(inf$compute_subsampling_confidence_interval)$b)

	m_dbg = inf$approximate_m_out_of_n_bootstrap_distribution_beta_hat_T(B = 7, show_progress = FALSE, debug = TRUE)
	expect_equal(m_dbg$m, expected_size)

	b_dbg = inf$approximate_subsampling_distribution_beta_hat_T(B = 7, show_progress = FALSE, debug = TRUE)
	expect_equal(b_dbg$b, expected_size)

	expect_error(
		suppressWarnings(inf$approximate_m_out_of_n_bootstrap_distribution_beta_hat_T(B = 7, m = "auto", show_progress = FALSE)),
		"m must satisfy",
		fixed = TRUE
	)
	expect_error(
		suppressWarnings(inf$approximate_subsampling_distribution_beta_hat_T(B = 7, b = "auto", show_progress = FALSE)),
		"b must satisfy",
		fixed = TRUE
	)
})

test_that("PRW subsampling methods are available and return scalar inference", {
	inf = build_resampling_smoke_inference()
	expect_true(is.function(inf$approximate_subsampling_distribution_beta_hat_T))
	expect_true(is.function(inf$compute_subsampling_two_sided_pval))
	expect_true(is.function(inf$compute_subsampling_confidence_interval))
	expect_true(is.function(inf$select_optimal_b_subsampling))
	expect_true(is.function(inf$compute_subsampling_sensitivity))

	d = inf$approximate_subsampling_distribution_beta_hat_T(B = 11, b = 8, show_progress = FALSE)
	expect_length(d, 11)
	expect_true(any(is.finite(d)))

	ci = inf$compute_subsampling_confidence_interval(B = 21, b = 8, show_progress = FALSE)
	expect_length(ci, 2)
	expect_true(all(is.finite(ci)))
	expect_lte(ci[1], ci[2])

	p = inf$compute_subsampling_two_sided_pval(B = 21, b = 8, show_progress = FALSE)
	expect_true(is.finite(p))
	expect_gte(p, 0)
	expect_lte(p, 1)
})

test_that("m-out-of-n and PRW CIs reject extreme finite intervals under hardening", {
	inf = build_resampling_smoke_inference()
	priv = inf$.__enclos_env__$private
	priv$harden = TRUE
	old_width_threshold = priv$bootstrap_extreme_ci_width_threshold
	on.exit({
		priv$bootstrap_extreme_ci_width_threshold = old_width_threshold
	})
	priv$bootstrap_extreme_ci_width_threshold = 1e-12

	m_ci = inf$compute_m_out_of_n_bootstrap_confidence_interval(B = 21, m = 8, show_progress = FALSE)
	expect_true(all(is.na(m_ci)))
	expect_true(inf$is_nonestimable("estimate"))
	expect_match(inf$get_nonestimable_reason(), "m_out_of_n_extreme_confidence_interval", fixed = TRUE)

	inf_sub = build_resampling_smoke_inference()
	priv_sub = inf_sub$.__enclos_env__$private
	priv_sub$harden = TRUE
	old_sub_width_threshold = priv_sub$bootstrap_extreme_ci_width_threshold
	on.exit({
		priv_sub$bootstrap_extreme_ci_width_threshold = old_sub_width_threshold
	}, add = TRUE)
	priv_sub$bootstrap_extreme_ci_width_threshold = 1e-12
	sub_ci = inf_sub$compute_subsampling_confidence_interval(B = 21, b = 8, show_progress = FALSE)
	expect_true(all(is.na(sub_ci)))
	expect_true(inf_sub$is_nonestimable("estimate"))
	expect_match(inf_sub$get_nonestimable_reason(), "subsampling_extreme_confidence_interval", fixed = TRUE)
})

test_that("minimum-volatility selector is shared by m and b methods", {
	inf = build_resampling_smoke_inference()

	m_sel = inf$select_optimal_m_out_of_n_bootstrap(
		B = 9,
		m_grid = c(6L, 7L, 8L),
		show_progress = FALSE,
		min_finite_fraction = 0
	)
	expect_s3_class(m_sel, "EDIMOutOfNBootstrapMSelection")
	expect_true(is.finite(m_sel$m_optimal))
	expect_true(m_sel$m_optimal %in% c(6L, 7L, 8L))
	expect_true(all(c("m", "ci_width", "objective_value", "volatility", "eligible") %in% names(m_sel$grid_table)))

	b_sel = inf$select_optimal_b_subsampling(
		B = 9,
		b_grid = c(6L, 7L, 8L),
		show_progress = FALSE,
		min_finite_fraction = 0
	)
	expect_s3_class(b_sel, "EDISubsamplingBSelection")
	expect_true(is.finite(b_sel$b_optimal))
	expect_true(b_sel$b_optimal %in% c(6L, 7L, 8L))
	expect_true(all(c("b", "ci_width", "objective_value", "volatility", "eligible") %in% names(b_sel$grid_table)))

	sens = inf$compute_subsampling_sensitivity(
		B = 9,
		b_grid = c(6L, 7L, 8L),
		show_progress = FALSE
	)
	expect_s3_class(sens, "EDISubsamplingSensitivity")
	expect_true("grid_table" %in% names(sens))
})

test_that("debug distributions include centered scaled diagnostics", {
	inf = build_resampling_smoke_inference()

	m_dbg = inf$approximate_m_out_of_n_bootstrap_distribution_beta_hat_T(B = 7, m = 8, show_progress = FALSE, debug = TRUE)
	expect_true(all(c("values", "centered_scaled_values", "finite_fraction", "m", "n_units") %in% names(m_dbg)))
	expect_length(m_dbg$values, 7)
	expect_length(m_dbg$centered_scaled_values, 7)

	b_dbg = inf$approximate_subsampling_distribution_beta_hat_T(B = 7, b = 8, show_progress = FALSE, debug = TRUE)
	expect_true(all(c("values", "centered_scaled_values", "finite_fraction", "b", "n_units") %in% names(b_dbg)))
	expect_length(b_dbg$values, 7)
	expect_length(b_dbg$centered_scaled_values, 7)
})
