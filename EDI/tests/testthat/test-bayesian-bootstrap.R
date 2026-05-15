context("Bayesian bootstrap")

make_seq_design_for_bayes_boot = function(response_type, y){
	des = DesignSeqOneByOneBernoulli$new(n = length(y), response_type = response_type)
	for (i in seq_along(y)) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = i / 10))
	}
	des$add_all_subject_responses(y)
	des
}

make_kk_design_for_weighted_bayes_boot = function(response_type, y, n_pairs = 3L, n_single = 2L){
	n = 2L * n_pairs + n_single
	des = DesignSeqOneByOneKK14$new(n = n, response_type = response_type, verbose = FALSE)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = i / 10, x2 = (i %% 3) / 10))
	}
	des$.__enclos_env__$private$m <- c(rep(seq_len(n_pairs), each = 2L), rep(0L, n_single))
	des$add_all_subject_responses(y)
	des
}

test_that("jackknife descendants expose Bayesian bootstrap methods", {
	des = make_seq_design_for_bayes_boot("count", c(0L, 1L, 1L, 2L, 3L, 1L, 0L, 2L))
	inf = InferenceCountPoisson$new(des)
	expect_true(inherits(inf, "InferenceJackknife"))
	expect_true(is.function(inf$approximate_bayesian_bootstrap_distribution_beta_hat_T))
	expect_true(is.function(inf$compute_bayesian_bootstrap_two_sided_pval))
	expect_true(is.function(inf$compute_bayesian_bootstrap_confidence_interval))
})

test_that("equal subject weights recover the original point estimate", {
	des = make_seq_design_for_bayes_boot("count", c(0L, 1L, 1L, 2L, 3L, 1L, 0L, 2L))
	inf = InferenceCountPoisson$new(des)
	est = inf$compute_estimate()
	n = des$get_n()
	inf$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n),
		unit_group_id = rep(1L, n),
		n_units = n
	)
	weighted_est = inf$compute_estimate_with_bootstrap_weights(rep(1, n))
	expect_equal(as.numeric(weighted_est), as.numeric(est), tolerance = 1e-5)
})

test_that("weighted logistic hook returns a finite estimate", {
	des = make_seq_design_for_bayes_boot("incidence", c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L))
	inf = InferenceIncidLogRegr$new(des)
	n = des$get_n()
	inf$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n),
		unit_group_id = rep(1L, n),
		n_units = n
	)
	weighted_est = inf$compute_estimate_with_bootstrap_weights(rep(1, n))
	expect_true(is.finite(as.numeric(weighted_est)))
})

test_that("Bayesian bootstrap smoke test runs on a stable first-wave family", {
	des = make_seq_design_for_bayes_boot("count", c(0L, 1L, 1L, 2L, 3L, 1L, 0L, 2L))
	inf = InferenceCountPoisson$new(des)
	inf$num_cores = 1L
	set.seed(20260515)
	boot = inf$approximate_bayesian_bootstrap_distribution_beta_hat_T(B = 8L, show_progress = FALSE)
	expect_length(boot, 8L)
	expect_true(all(is.finite(boot)))
	set.seed(20260515)
	ci = inf$compute_bayesian_bootstrap_confidence_interval(B = 8L, show_progress = FALSE, type = "percentile")
	expect_length(ci, 2L)
	expect_true(all(is.finite(ci)))
})

test_that("next-wave weighted hooks return finite estimates on simple and g-computation paths", {
	des_cont = make_seq_design_for_bayes_boot("continuous", c(0, 1, 2, 3, 4, 5, 6, 7))
	inf_cont = InferenceAllSimpleMeanDiff$new(des_cont)
	n_cont = des_cont$get_n()
	inf_cont$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_cont),
		unit_group_id = rep(1L, n_cont),
		n_units = n_cont
	)
	expect_true(is.finite(as.numeric(inf_cont$compute_estimate_with_bootstrap_weights(rep(1, n_cont)))))

	des_incid = make_seq_design_for_bayes_boot("incidence", c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L))
	inf_incid = InferenceIncidGCompRiskDiff$new(des_incid)
	n_incid = des_incid$get_n()
	inf_incid$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_incid),
		unit_group_id = rep(1L, n_incid),
		n_units = n_incid
	)
	expect_true(is.finite(as.numeric(inf_incid$compute_estimate_with_bootstrap_weights(rep(1, n_incid)))))

	des_prop = make_seq_design_for_bayes_boot("proportion", c(0.1, 0.8, 0.2, 0.7, 0.6, 0.3, 0.9, 0.4))
	inf_prop = InferencePropGCompMeanDiff$new(des_prop)
	n_prop = des_prop$get_n()
	inf_prop$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_prop),
		unit_group_id = rep(1L, n_prop),
		n_units = n_prop
	)
	expect_true(is.finite(as.numeric(inf_prop$compute_estimate_with_bootstrap_weights(rep(1, n_prop)))))

	des_ord = make_seq_design_for_bayes_boot("ordinal", c(1L, 2L, 1L, 3L, 2L, 1L, 3L, 2L))
	inf_ord = InferenceOrdinalGCompMeanDiff$new(des_ord)
	n_ord = des_ord$get_n()
	inf_ord$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_ord),
		unit_group_id = rep(1L, n_ord),
		n_units = n_ord
	)
	expect_true(is.finite(as.numeric(inf_ord$compute_estimate_with_bootstrap_weights(rep(1, n_ord)))))
})

test_that("next-wave weighted hooks return finite estimates on KK GEE paths", {
	skip_if_not_installed("geepack")

	des_incid = make_kk_design_for_weighted_bayes_boot("incidence", c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L))
	inf_incid = InferenceIncidKKGEE$new(des_incid, use_rcpp = TRUE, verbose = FALSE)
	n_incid = des_incid$get_n()
	inf_incid$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_incid),
		unit_group_id = rep(1L, n_incid),
		n_units = n_incid
	)
	expect_true(is.finite(as.numeric(inf_incid$compute_estimate_with_bootstrap_weights(rep(1, n_incid)))))

	des_count = make_kk_design_for_weighted_bayes_boot("count", c(0L, 1L, 1L, 2L, 3L, 1L, 0L, 2L))
	inf_count = InferenceCountPoissonKKGEE$new(des_count, use_rcpp = TRUE, verbose = FALSE)
	n_count = des_count$get_n()
	inf_count$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_count),
		unit_group_id = rep(1L, n_count),
		n_units = n_count
	)
	expect_true(is.finite(as.numeric(inf_count$compute_estimate_with_bootstrap_weights(rep(1, n_count)))))

	des_prop = make_kk_design_for_weighted_bayes_boot("proportion", c(0.1, 0.8, 0.2, 0.7, 0.6, 0.3, 0.9, 0.4))
	inf_prop = InferencePropKKGEE$new(des_prop, use_rcpp = TRUE, verbose = FALSE)
	n_prop = des_prop$get_n()
	inf_prop$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_prop),
		unit_group_id = rep(1L, n_prop),
		n_units = n_prop
	)
	expect_true(is.finite(as.numeric(inf_prop$compute_estimate_with_bootstrap_weights(rep(1, n_prop)))))

	des_ord = make_kk_design_for_weighted_bayes_boot("ordinal", c(1L, 2L, 1L, 3L, 2L, 1L, 3L, 2L))
	inf_ord = InferenceOrdinalKKGEE$new(des_ord, verbose = FALSE)
	n_ord = des_ord$get_n()
	inf_ord$.__enclos_env__$private$current_bayesian_bootstrap_context = list(
		row_to_unit = seq_len(n_ord),
		unit_group_id = rep(1L, n_ord),
		n_units = n_ord
	)
	expect_true(is.finite(as.numeric(inf_ord$compute_estimate_with_bootstrap_weights(rep(1, n_ord)))))
})

test_that("Bayesian bootstrap matches mirai-backed parallel execution", {
	skip_if_not_installed("mirai")

	on.exit(unset_num_cores(), add = TRUE)

	des = make_seq_design_for_bayes_boot("count", c(0L, 1L, 1L, 2L, 3L, 1L, 0L, 2L))
	serial_inf = InferenceCountPoisson$new(des)
	mirai_inf = InferenceCountPoisson$new(des)

	serial_inf$num_cores = 1L
	set.seed(20260515)
	serial_boot = serial_inf$approximate_bayesian_bootstrap_distribution_beta_hat_T(
		B = 11L,
		show_progress = FALSE
	)

	set_num_cores(2L, force_mirai = TRUE)
	mirai_inf$num_cores = 2L
	set.seed(20260515)
	mirai_boot = mirai_inf$approximate_bayesian_bootstrap_distribution_beta_hat_T(
		B = 11L,
		show_progress = FALSE
	)

	expect_equal(unname(mirai_boot), unname(serial_boot), tolerance = 1e-10)
})
