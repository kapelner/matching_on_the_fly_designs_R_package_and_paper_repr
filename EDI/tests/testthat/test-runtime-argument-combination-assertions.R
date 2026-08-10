test_that("bootstrap argument helper enforces min_number_usable_samples within B", {
	expect_error(
		EDI:::assertBootstrapArgs(B = 2L, min_number_usable_samples = 3L, context = "unit test bootstrap"),
		"min_number_usable_samples.*B"
	)
	expect_silent(EDI:::assertBootstrapArgs(B = 3L, min_number_usable_samples = 3L, context = "unit test bootstrap"))

	old = EDI:::should_run_asserts()
	on.exit(EDI::toggle_asserts(old), add = TRUE)
	EDI::toggle_asserts(FALSE)
	expect_silent(EDI:::assertBootstrapArgs(B = 2L, min_number_usable_samples = 3L, context = "unit test bootstrap"))
})

test_that("public bootstrap methods use promoted runtime argument checks", {
	des = DesignFixedBernoulli$new(n = 4, response_type = "continuous", seed = 20260810L)
	des$add_all_subjects_to_experiment(data.frame(x = 1:4))
	des$assign_w_to_all_subjects(w_precomputed = c(0, 1, 0, 1))
	des$add_all_subject_responses(c(1, 2, 3, 4))
	inf = InferenceAllSimpleMeanDiff$new(des)

	expect_error(
		inf$compute_bootstrap_confidence_interval(B = 2L, min_number_usable_samples = 3L, show_progress = FALSE),
		"min_number_usable_samples.*B"
	)
	expect_error(
		inf$compute_bootstrap_two_sided_pval(B = 2L, min_number_usable_samples = 3L, show_progress = FALSE),
		"min_number_usable_samples.*B"
	)
})

test_that("formula context helper gives direct missing-variable errors", {
	expect_silent(EDI:::assertFormulaContext(~ x + z, data.frame(x = 1, z = 2), context = "model_formula"))
	expect_error(
		EDI:::assertFormulaContext(~ x + missing_col, data.frame(x = 1), context = "model_formula"),
		"missing_col"
	)

	des = DesignFixedBernoulli$new(n = 4, response_type = "continuous", seed = 20260810L)
	des$add_all_subjects_to_experiment(data.frame(x = 1:4))
	des$assign_w_to_all_subjects(w_precomputed = c(0, 1, 0, 1))
	des$add_all_subject_responses(c(1, 2, 3, 4))
	expect_error(
		InferenceAllSimpleMeanDiff$new(des, model_formula = ~ missing_col),
		"missing_col"
	)
})

test_that("strata and cluster helper catches promoted runtime relationships", {
	expect_error(
		EDI:::assertStrataClusterArgs(strata_cols = "block", cluster_col = "block", context = "unit test design"),
		"cluster_col.*strata_cols"
	)
	expect_error(
		EDI:::assertStrataClusterArgs(strata_cols = "block", cluster_col = "cluster", data = data.frame(block = factor("a")), context = "unit test design"),
		"cluster"
	)
	expect_error(
		EDI:::assertStrataClusterArgs(strata_cols = "block", data = data.frame(block = 1), strata_cols_must_be_factor = TRUE, context = "unit test design"),
		"non-categorical"
	)
	expect_silent(
		EDI:::assertStrataClusterArgs(strata_cols = "block", cluster_col = "cluster", data = data.frame(block = "a", cluster = "c"), strata_cols_must_be_factor = TRUE, context = "unit test design")
	)
})

test_that("public design paths use promoted strata and cluster runtime checks", {
	expect_error(
		DesignFixedBlockedCluster$new(strata_cols = "group_id", cluster_col = "group_id", response_type = "continuous", n = 4),
		"cluster_col.*strata_cols"
	)

	des_cluster = DesignFixedCluster$new(cluster_col = "cluster_id", response_type = "continuous", n = 4)
	expect_error(
		des_cluster$add_all_subjects_to_experiment(data.frame(x = 1:4)),
		"cluster_id"
	)

	des_seq = DesignSeqOneByOneSPBR$new(strata_cols = "stratum", response_type = "continuous", n = 4)
	expect_error(
		des_seq$add_one_subject_to_experiment_and_assign(data.frame(stratum = 1)),
		"non-categorical|Continuous covariates are not allowed"
	)
})
