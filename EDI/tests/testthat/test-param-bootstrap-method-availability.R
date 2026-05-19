test_that("unsupported families do not expose parametric LR bootstrap methods", {
	make_bernoulli_design = function(response_type, responses, n = length(responses)) {
		des = DesignSeqOneByOneBernoulli$new(n = n, response_type = response_type)
		for (i in seq_len(n)) {
			des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
		}
		des$add_all_subject_responses(responses)
		des
	}

	make_kk_design = function(responses) {
		des = DesignSeqOneByOneKK14$new(n = length(responses), response_type = "continuous")
		for (i in seq_along(responses)) {
			des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
		}
		des$add_all_subject_responses(responses)
		des
	}

	expect_no_param_bootstrap_methods = function(obj) {
		expect_false("compute_lik_ratio_bootstrap_two_sided_pval" %in% names(obj))
		expect_false("compute_lik_ratio_bootstrap_confidence_interval" %in% names(obj))
		expect_null(obj$compute_lik_ratio_bootstrap_two_sided_pval)
		expect_null(obj$compute_lik_ratio_bootstrap_confidence_interval)
	}

	unsupported_objects = list(
		InferenceIncidRiskDiff$new(
			make_bernoulli_design("incidence", c(0, 1, 0, 1, 1, 0))
		),
		InferencePropFractionalLogit$new(
			make_bernoulli_design("proportion", c(0.2, 0.4, 0.6, 0.8, 0.3, 0.7))
		),
		InferenceCountHurdleNegBin$new(
			make_bernoulli_design("count", c(0, 1, 2, 0, 3, 1))
		),
		InferenceAllKKMeanDiffIVWC$new(
			make_kk_design(c(2.1, 1.4, 3.2, 0.8, 2.7, 1.9))
		)
	)

	for (obj in unsupported_objects) {
		expect_no_param_bootstrap_methods(obj)
	}
})

test_that("supported families still expose parametric LR bootstrap methods", {
	des = DesignSeqOneByOneBernoulli$new(n = 6, response_type = "incidence")
	for (i in 1:6) {
		des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
	}
	des$add_all_subject_responses(c(0, 1, 0, 1, 1, 0))
	obj = InferenceIncidLogRegr$new(des)

	expect_true("compute_lik_ratio_bootstrap_two_sided_pval" %in% names(obj))
	expect_true("compute_lik_ratio_bootstrap_confidence_interval" %in% names(obj))
	expect_true(is.function(obj$compute_lik_ratio_bootstrap_two_sided_pval))
	expect_true(is.function(obj$compute_lik_ratio_bootstrap_confidence_interval))
})
