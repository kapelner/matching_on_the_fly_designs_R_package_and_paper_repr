library(testthat)
library(EDI)

capability_test_metadata = function(likelihood_tier = "none", capabilities = character()) {
	list(
		abstract = FALSE,
		exported = FALSE,
		response_types = "continuous",
		design_families = "all",
		compatibility = EDI:::always_compatible_inference_metadata,
		likelihood_tier = likelihood_tier,
		required_packages = character(),
		capabilities = capabilities
	)
}

likelihood_ratio_public = function() {
	list(
		compute_lik_ratio_two_sided_pval = function(delta = 0) NA_real_,
		compute_lik_ratio_confidence_interval = function(alpha = 0.05) c(NA_real_, NA_real_)
	)
}

parametric_likelihood_bootstrap_public = function() {
	c(likelihood_ratio_public(), list(
		compute_lik_ratio_bootstrap_two_sided_pval = function(delta = 0, B = 199) NA_real_,
		compute_lik_ratio_bootstrap_confidence_interval = function(alpha = 0.05, B = 199) c(NA_real_, NA_real_),
		get_last_param_bootstrap_diagnostics = function() NULL
	))
}

likelihood_ratio_private = function(include_null_simulator = FALSE) {
	private = list(
		get_likelihood_test_spec = function() list()
	)
	if (isTRUE(include_null_simulator)) {
		private$simulate_under_lik_null = function(spec, delta, null_fit) list()
	}
	private
}

test_that("capability tables declare requirements and public API methods", {
	expect_true(is.list(EDI:::capability_requires))
	expect_true(is.list(EDI:::public_methods_for_capability))
	expect_true(all(c(
		"exact_test",
		"exact_binomial_incidence",
		"exact_fisher_incidence",
		"exact_zhang_incidence",
		"likelihood_ratio",
		"estimating_equation_likelihood_ratio",
		"parametric_likelihood_bootstrap"
	) %in% names(EDI:::capability_requires)))
	expect_true(all(c(
		"exact_test",
		"likelihood_ratio",
		"estimating_equation_likelihood_ratio",
		"parametric_likelihood_bootstrap"
	) %in% names(EDI:::public_methods_for_capability)))
	expect_identical(
		EDI:::capability_requires$parametric_likelihood_bootstrap$capabilities,
		"likelihood_ratio"
	)
	expect_true("simulate_under_lik_null" %in%
		EDI:::capability_requires$parametric_likelihood_bootstrap$private_methods)
	expect_true(all(c(
		"compute_exact_confidence_interval",
		"compute_exact_two_sided_pval_for_treatment_effect"
	) %in% EDI:::public_methods_for_capability$exact_test))
})

test_that("root supports and capabilities are metadata queries", {
	on.exit({
		EDI:::populate_inference_component_registry()
		EDI:::populate_inference_class_registry()
	}, add = TRUE)
	EDI:::register_inference_component(EDI:::InferenceComponent(
		name = "InferenceTemporaryQueryComponent",
		status = "active",
		file = "test",
		provides_capabilities = "temporary_component_query"
	))
	gen = EDI:::define_inference_class(
		classname = "InferenceTemporaryCapabilityQueryHost",
		inherit = Inference,
		components = "InferenceTemporaryQueryComponent",
		public = list(initialize = function() invisible(self)),
		metadata = list(likelihood_tier = "none")
	)
	EDI:::register_inference_class(
		name = "InferenceTemporaryCapabilityQueryHost",
		parent = "Inference",
		metadata = capability_test_metadata(
			likelihood_tier = "none",
			capabilities = "temporary_class_query"
		),
		direct_components = "InferenceTemporaryQueryComponent"
	)
	obj = gen$new()

	expect_setequal(obj$capabilities(), c("temporary_component_query", "temporary_class_query"))
	expect_identical(
		obj$supports(c("temporary_component_query", "temporary_class_query", "missing")),
		c(
			temporary_component_query = TRUE,
			temporary_class_query = TRUE,
			missing = FALSE
		)
	)
})

test_that("likelihood_tier none exposes no likelihood APIs on factory classes", {
	likelihood_api = unique(unlist(EDI:::public_methods_for_capability[c(
		"likelihood_ratio",
		"estimating_equation_likelihood_ratio",
		"parametric_likelihood_bootstrap"
	)], use.names = FALSE))
	gen = EDI:::define_inference_class(
		classname = "InferenceTemporaryNoLikelihoodApis",
		metadata = list(likelihood_tier = "none")
	)

	expect_false(any(likelihood_api %in% names(gen$public_methods)))
	expect_error(
		EDI:::define_inference_class(
			classname = "InferenceTemporaryNoneWithLikelihoodRatio",
			public = likelihood_ratio_public(),
			private = likelihood_ratio_private(),
			metadata = list(
				likelihood_tier = "none",
				capabilities = "likelihood_ratio"
			)
		),
		"disallowed likelihood tier"
	)
})

test_that("quasi likelihood tier requires estimating-equation likelihood-ratio capability", {
	expect_error(
		EDI:::define_inference_class(
			classname = "InferenceTemporaryQuasiWithLikelihoodRatio",
			public = likelihood_ratio_public(),
			private = likelihood_ratio_private(),
			metadata = list(
				likelihood_tier = "quasi",
				capabilities = "likelihood_ratio"
			)
		),
		"disallowed likelihood tier"
	)
	expect_silent(EDI:::define_inference_class(
		classname = "InferenceTemporaryQuasiEstimatingEquationLR",
		public = likelihood_ratio_public(),
		private = likelihood_ratio_private(),
		metadata = list(
			likelihood_tier = "quasi",
			capabilities = "estimating_equation_likelihood_ratio"
		)
	))
})

test_that("partial and full likelihood tiers do not imply parametric bootstrap", {
	bootstrap_api = EDI:::public_methods_for_capability$parametric_likelihood_bootstrap
	gen = EDI:::define_inference_class(
		classname = "InferenceTemporaryFullLikelihoodNoParamBootstrap",
		public = likelihood_ratio_public(),
		private = likelihood_ratio_private(),
		metadata = list(
			likelihood_tier = "full",
			capabilities = "likelihood_ratio"
		)
	)

	expect_false(any(bootstrap_api %in% names(gen$public_methods)))
	expect_error(
		EDI:::define_inference_class(
			classname = "InferenceTemporaryFullLikelihoodMissingNullSimulator",
			public = parametric_likelihood_bootstrap_public(),
			private = likelihood_ratio_private(include_null_simulator = FALSE),
			metadata = list(
				likelihood_tier = "full",
				capabilities = c("likelihood_ratio", "parametric_likelihood_bootstrap")
			)
		),
		"simulate_under_lik_null"
	)
	expect_silent(EDI:::define_inference_class(
		classname = "InferenceTemporaryFullLikelihoodParamBootstrap",
		public = parametric_likelihood_bootstrap_public(),
		private = likelihood_ratio_private(include_null_simulator = TRUE),
		metadata = list(
			likelihood_tier = "full",
			capabilities = c("likelihood_ratio", "parametric_likelihood_bootstrap")
		)
	))

	expect_silent(EDI:::define_inference_class(
		classname = "InferenceTemporaryPartialLikelihoodNoParamBootstrap",
		public = likelihood_ratio_public(),
		private = likelihood_ratio_private(),
		metadata = list(
			likelihood_tier = "partial",
			capabilities = "likelihood_ratio"
		)
	))
})
