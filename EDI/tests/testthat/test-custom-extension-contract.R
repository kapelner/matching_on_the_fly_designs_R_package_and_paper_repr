test_that("custom extension base classes are internal but usable for subclassing", {
	exports = getNamespaceExports("EDI")
	expect_false("InferenceCustomAsymp" %in% exports)
	expect_false("InferenceCustomRand" %in% exports)
	expect_false("InferenceCustomBoot" %in% exports)
	expect_false("DesignCustomSequential" %in% exports)
	expect_false("DesignCustomFixed" %in% exports)

	expect_true(is.environment(getNamespace("EDI")))
	expect_s3_class(getFromNamespace("InferenceCustomAsymp", "EDI"), "R6ClassGenerator")
	expect_s3_class(getFromNamespace("InferenceCustomRand", "EDI"), "R6ClassGenerator")
	expect_s3_class(getFromNamespace("InferenceCustomBoot", "EDI"), "R6ClassGenerator")
	expect_s3_class(getFromNamespace("DesignCustomSequential", "EDI"), "R6ClassGenerator")
	expect_s3_class(getFromNamespace("DesignCustomFixed", "EDI"), "R6ClassGenerator")
})

test_that("custom asymptotic inference works from an external-package-like environment", {
	ext_env = new.env(parent = globalenv())
	ext_env$R6Class = R6::R6Class
	ext_env$InferenceCustomAsymp = getFromNamespace("InferenceCustomAsymp", "EDI")

	evalq({
		ExternalMeanDiff = R6Class(
			"ExternalMeanDiff",
			inherit = InferenceCustomAsymp,
			public = list(
				fit = function(estimate_only = FALSE) {
					dat = self$get_analysis_data()
					y_t = dat$y[dat$w == 1]
					y_c = dat$y[dat$w == 0]
					est = mean(y_t) - mean(y_c)
					if (estimate_only) return(list(estimate = est))
					se = sqrt(stats::var(y_t) / length(y_t) + stats::var(y_c) / length(y_c))
					df = length(y_t) + length(y_c) - 2
					list(
						estimate = est,
						se = se,
						df = df,
						is_z = FALSE,
						model = list(n_t = length(y_t), n_c = length(y_c))
					)
				}
			)
		)
	}, envir = ext_env)

	des = FixedDesignBernoulli$new(n = 20, response_type = "continuous", verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = seq_len(20)))
	des$overwrite_all_subject_assignments(rep(c(0, 1), each = 10))
	des$add_all_subject_responses(c(1:10, 12:21))

	inf = ext_env$ExternalMeanDiff$new(des, verbose = FALSE)
	expect_s3_class(inf, "ExternalMeanDiff")
	expect_equal(inf$get_response(), c(1:10, 12:21))
	expect_equal(inf$get_treatment(), rep(c(0, 1), each = 10))
	expect_equal(inf$get_response_type(), "continuous")
	expect_identical(inf$get_design_object(), des)
	expect_equal(nrow(inf$get_analysis_data()), 20)
	expect_true("x" %in% names(inf$get_analysis_data()))

	expect_equal(inf$compute_estimate(), 11)
	expect_true(is.finite(inf$compute_asymp_two_sided_pval()))
	expect_length(inf$compute_asymp_confidence_interval(), 2)
	expect_equal(inf$get_mod()$n_t, 10)

	set.seed(20260420)
	boot_p = inf$compute_bootstrap_two_sided_pval(
		B = 21,
		na.rm = TRUE,
		min_number_usable_samples = 5L
	)
	expect_true(is.finite(boot_p))
	expect_true(boot_p >= 0 && boot_p <= 1)
})

test_that("custom design extension bases delegate user assignment rules", {
	CustomFixedBase = getFromNamespace("DesignCustomFixed", "EDI")
	CustomSequentialBase = getFromNamespace("DesignCustomSequential", "EDI")

	AlternatingFixed = R6::R6Class(
		"AlternatingFixed",
		inherit = CustomFixedBase,
		public = list(
			draw_assignments = function(r = 1) {
				matrix(rep(rep(c(0, 1), length.out = self$get_n()), r), nrow = self$get_n(), ncol = r)
			}
		)
	)
	fixed = AlternatingFixed$new(n = 6, response_type = "continuous", verbose = FALSE)
	fixed$add_all_subjects_to_experiment(data.frame(x = 1:6))
	fixed$assign_w_to_all_subjects()
	expect_equal(fixed$get_w(), c(0, 1, 0, 1, 0, 1))
	expect_equal(dim(fixed$draw_ws_according_to_design(3)), c(6, 3))

	ExternalSequential = R6::R6Class(
		"ExternalSequential",
		inherit = CustomSequentialBase,
		public = list(
			assignment_rule = function() {
				as.numeric(self$get_t() %% 2L == 0L)
			}
		)
	)
	seq_des = ExternalSequential$new(n = 4, response_type = "continuous", verbose = FALSE)
	for (i in 1:4) {
		seq_des$add_one_subject_to_experiment_and_assign(data.frame(x = i))
	}
	expect_equal(seq_des$get_w(), c(0, 1, 0, 1))
})
