test_that("jackknife aliases and public distribution API are coherent", {
	des <- DesignSeqOneByOneBernoulli$new(n = 8, response_type = "continuous")
	for (i in seq_len(8)) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x1 = i / 10))
	}
	des$add_all_subject_responses(c(0, 1, 2, 3, 4, 5, 6, 7))

	inf <- InferenceAllSimpleMeanDiff$new(des, verbose = FALSE)

	jack <- inf$approximate_jackknife_distribution_beta_hat_T()
	expect_length(jack, des$get_n())

	full_est <- inf$compute_estimate()
	jack_bar <- mean(jack)
	manual_bias <- (length(jack) - 1) * (jack_bar - full_est)
	manual_est <- full_est - manual_bias
	manual_se <- sqrt(((length(jack) - 1) / length(jack)) * sum((jack - jack_bar)^2))

	expect_equal(inf$compute_jackknife_bias_estimate(), manual_bias, tolerance = 1e-12)
	expect_equal(inf$compute_jackknife_estimate(), manual_est, tolerance = 1e-12)
	expect_equal(inf$compute_jackknife_corrected_estimate(), manual_est, tolerance = 1e-12)
	expect_equal(inf$compute_jackknife_std_error(), manual_se, tolerance = 1e-12)
	expect_equal(inf$compute_jackknife_standard_error(), manual_se, tolerance = 1e-12)
})

test_that("clustered designs support cluster jackknife deletion", {
	des <- DesignFixedCluster$new(
		n = 8,
		response_type = "continuous",
		cluster_col = "cl",
		verbose = FALSE
	)
	des$add_all_subjects_to_experiment(
		data.frame(
			x1 = seq_len(8),
			cl = factor(rep(1:4, each = 2), levels = 1:4)
		)
	)
	des$overwrite_all_subject_assignments(rep(c(1, 1, 0, 0), each = 2))
	des$add_all_subject_responses(c(8, 10, 1, 3, 9, 11, 2, 4))

	inf <- InferenceAllSimpleMeanDiff$new(des, verbose = FALSE)
	cluster_jack <- inf$approximate_jackknife_distribution_beta_hat_T(unit = "cluster")
	expect_length(cluster_jack, 4L)

	cluster_ids <- as.character(des$get_X_raw()$cl)
	full_est <- inf$compute_estimate()
	manual_jack <- vapply(unique(cluster_ids), function(cluster_id) {
		keep <- cluster_ids != cluster_id
		y <- des$get_y()[keep]
		w <- des$get_w()[keep]
		mean(y[w == 1]) - mean(y[w == 0])
	}, numeric(1))
	jack_bar <- mean(manual_jack)
	manual_bias <- (length(manual_jack) - 1) * (jack_bar - full_est)
	manual_est <- full_est - manual_bias
	manual_se <- sqrt(((length(manual_jack) - 1) / length(manual_jack)) * sum((manual_jack - jack_bar)^2))

	expect_equal(cluster_jack, manual_jack, tolerance = 1e-12)
	expect_equal(inf$compute_jackknife_bias_estimate(unit = "cluster"), manual_bias, tolerance = 1e-12)
	expect_equal(inf$compute_jackknife_estimate(unit = "cluster"), manual_est, tolerance = 1e-12)
	expect_equal(inf$compute_jackknife_standard_error(unit = "cluster"), manual_se, tolerance = 1e-12)
})
