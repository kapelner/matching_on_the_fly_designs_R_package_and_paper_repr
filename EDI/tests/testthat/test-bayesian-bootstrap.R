context("Bayesian bootstrap")

make_seq_design_for_bayes_boot = function(response_type, y){
	des = DesignSeqOneByOneBernoulli$new(n = length(y), response_type = response_type)
	for (i in seq_along(y)) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = i / 10))
	}
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
	expect_equal(as.numeric(weighted_est), as.numeric(est), tolerance = 1e-8)
})

test_that("Bayesian bootstrap smoke test runs on supported GLM paths", {
	des = make_seq_design_for_bayes_boot("incidence", c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L))
	inf = InferenceIncidLogRegr$new(des)
	boot = inf$approximate_bayesian_bootstrap_distribution_beta_hat_T(B = 8, show_progress = FALSE)
	expect_length(boot, 8L)
	expect_true(all(is.finite(boot)))
	ci = inf$compute_bayesian_bootstrap_confidence_interval(B = 8, show_progress = FALSE, type = "percentile")
	expect_length(ci, 2L)
	expect_true(all(is.finite(ci)))
})
