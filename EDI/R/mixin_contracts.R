#' Internal Mixin Host Contracts
#'
#' Pattern-1 mixins are lists spliced into an R6 class's \code{public} and
#' \code{private} lists. \code{EDI_MIXIN_CONTRACTS} documents the host-private
#' methods and state which a mixin requires but does not define itself.
#' Empty vectors mean that the mixin is self-contained (or is an intentionally
#' empty future extension point). These contracts are deliberately narrow: a
#' method supplied by the mixin itself is not repeated as a host requirement.
#'
#' Single-class-use extensions (\code{InferenceExt*}) are file-splits, not
#' reusable mixins, and are deliberately not tracked here -- this registry
#' exists to guard against silent method-name collisions when two or more
#' mixins are spliced into the same host, which cannot happen for an
#' extension spliced into exactly one class.
#'
#' \code{EDI_MIXIN_COMPOSITIONS} lists every base class that combines two or
#' more mixins. \code{EDI_MIXIN_ALLOWED_COLLISIONS} records the sole deliberate
#' overwrite: the compound KK mixin replaces the pass-through implementation of
#' \code{compute_basic_match_data()}. Tests use these objects to guard against
#' silent method-name overwrites as new mixins are added.
#'
#' @keywords internal
#' @noRd
EDI_MIXIN_CONTRACTS = list(
	InferenceMixinCordeiroFerrariApprox = list(
		file = "inference_mixin_cordeiro_ferrari_approx.R",
		private_methods = character(),
		private_state = character()
	),
	InferenceMixinKKGEEShared = list(
		file = "inference_mixin_kk_gee_shared.R",
		private_methods = c(
			"cache_nonestimable_estimate", "cache_nonestimable_se", "clear_nonestimable_state",
			"compute_z_or_t_ci_from_s_and_df", "compute_z_or_t_two_sided_pval_from_s_and_df",
			"create_design_matrix", "expand_subject_or_block_weights_to_row_weights",
			"fit_with_hardened_qr_column_dropping", "gee_family", "gee_response_type",
			"get_fit_warm_start_fisher", "get_fit_warm_start_for_length",
			"get_fit_warm_start_weights", "set_fit_warm_start", "shared_gee_dispatch"
		),
		private_state = c("any_censoring", "cached_values", "harden", "n", "w", "y")
	),
	InferenceMixinKKGLMMShared = list(
		file = "inference_mixin_kk_glmm_shared.R",
		private_methods = c(
			"cache_nonestimable_estimate", "compute_standard_error_from_information_matrix",
			"compute_z_or_t_ci_from_s_and_df", "compute_z_or_t_two_sided_pval_from_s_and_df",
			"create_design_matrix", "expand_subject_or_block_weights_to_row_weights",
			"fit_with_hardened_qr_column_dropping", "glmm_family", "glmm_response_type",
			"shared"
		),
		private_state = c("any_censoring", "cached_values", "harden", "n", "w", "y")
	),
	InferenceMixinKKPassThrough = list(
		file = "inference_mixin_kk_passthrough.R",
		private_methods = c(
			"assert_valid_bootstrap_type", "cache_nonestimable_estimate",
			"effective_parallel_cores",
			"expand_subject_or_block_weights_to_row_weights", "get_X", "has_private_method",
			"object_has_private_method", "par_lapply"
		),
		private_state = c("cached_values", "des_obj_priv_int", "has_match_structure", "n")
	),
	InferenceMixinKKPassThroughCompound = list(
		file = "inference_mixin_kk_passthrough_compound.R",
		private_methods = c("cache_nonestimable_estimate", "compute_basic_kk_match_data_impl"),
		private_state = c("cached_values", "has_match_structure")
	),
	InferenceMixinLemonteGradientApprox = list(
		file = "inference_mixin_lemonte_gradient_approx.R",
		private_methods = character(),
		private_state = character()
	),
	InferenceMixinOffOptimumLikelihoodEval = list(
		file = "inference_mixin_off_optimum_likelihood_eval.R",
		private_methods = c("get_default_information_source", "get_likelihood_test_spec"),
		private_state = character()
	)
)

EDI_MIXIN_COMPOSITIONS = list(
	InferenceKKPassThroughCompound = c(
		"InferenceMixinKKPassThrough", "InferenceMixinKKPassThroughCompound"
	),
	InferenceKKPassThroughCompoundNoParamBootstrap = c(
		"InferenceMixinKKPassThrough", "InferenceMixinKKPassThroughCompound"
	)
)

EDI_MIXIN_ALLOWED_COLLISIONS = list(
	InferenceKKPassThroughCompound = list(
		private = "compute_basic_match_data",
		public = character()
	),
	InferenceKKPassThroughCompoundNoParamBootstrap = list(
		private = "compute_basic_match_data",
		public = character()
	)
)

EDI_MIXIN_ALLOWED_OVERRIDES = list(
	InferenceKKPassThroughCompound = list(
		private = character(),
		public = c("approximate_bootstrap_distribution_beta_hat_T", "compute_estimate_with_bootstrap_weights")
	),
	InferenceKKPassThroughCompoundNoParamBootstrap = list(
		private = character(),
		public = c("approximate_bootstrap_distribution_beta_hat_T", "compute_estimate_with_bootstrap_weights")
	)
)

EDI_MIXIN_DEPENDENCIES = list(
	InferenceMixinKKPassThroughCompound = "InferenceMixinKKPassThrough"
)

EDI_LEGACY_MIXIN_COMPONENT_NAMES = c(
	InferenceMixinKKGEEShared = "KKGEE",
	InferenceMixinKKGLMMShared = "KKGLMM",
	InferenceMixinKKPassThrough = "KKPassThrough",
	InferenceMixinKKPassThroughCompound = "KKCompound",
	InferenceMixinOffOptimumLikelihoodEval = "OffOptimumLikelihoodEval"
)

EDI_INFERENCE_COMPONENTS = new.env(parent = emptyenv())

EDI_COMPONENT_ALLOWED_STATUSES = c("active", "scaffold")

EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS = c("none", "quasi", "partial", "full")

capability_requires = list(
	exact_test = list(),
	exact_binomial_incidence = list(
		capabilities = "exact_test"
	),
	exact_fisher_incidence = list(
		capabilities = "exact_test"
	),
	exact_zhang_incidence = list(
		capabilities = "exact_test"
	),
	randomization_test = list(),
	randomization_ci = list(
		capabilities = "randomization_test"
	),
	randomization_bootstrap = list(
		capabilities = c("randomization_test", "nonparametric_bootstrap")
	),
	jackknife = list(),
	wald = list(),
	likelihood_tests = list(
		likelihood_tier = c("quasi", "partial", "full"),
		private_methods = "get_likelihood_test_spec"
	),
	likelihood_ratio = list(
		likelihood_tier = c("partial", "full"),
		private_methods = "get_likelihood_test_spec"
	),
	estimating_equation_likelihood_ratio = list(
		likelihood_tier = "quasi",
		private_methods = "get_likelihood_test_spec"
	),
	parametric_likelihood_bootstrap = list(
		likelihood_tier = c("partial", "full"),
		capabilities = "likelihood_ratio",
		private_methods = c("get_likelihood_test_spec", "simulate_under_lik_null")
	),
	bayesian_bootstrap = list(
		public_methods = "compute_estimate_with_bootstrap_weights"
	),
	nonparametric_bootstrap = list(
		public_methods = "compute_estimate"
	),
	standard_model_cache = list(
		capabilities = "likelihood_tests"
	),
	count_likelihood_plumbing = list(
		likelihood_tier = c("quasi", "full"),
		capabilities = "likelihood_tests"
	),
	bartlett_approximation = list(
		capabilities = "parametric_likelihood_bootstrap",
		private_methods = c("run_param_bootstrap_replicates", "param_bootstrap_lr_extreme")
	),
	kk_passthrough = list(),
	kk_compound = list(),
	kk_gee = list(likelihood_tier = "quasi"),
	kk_glmm = list(likelihood_tier = c("partial", "full")),
	off_optimum_likelihood_eval = list(
		likelihood_tier = c("partial", "full"),
		private_methods = c("get_default_information_source", "get_likelihood_test_spec")
	)
)

public_methods_for_capability = list(
	exact_test = c(
		"compute_exact_confidence_interval",
		"compute_exact_two_sided_pval_for_treatment_effect"
	),
	randomization_test = c(
		"approximate_randomization_distribution_beta_hat_T",
		"compute_rand_two_sided_pval"
	),
	randomization_ci = c(
		"compute_rand_confidence_interval"
	),
	randomization_bootstrap = c(
		"approximate_rand_bootstrap_distribution_beta_hat_T",
		"compute_rand_bootstrap_two_sided_pval"
	),
	jackknife = c(
		"approximate_jackknife_distribution_beta_hat_T",
		"compute_jackknife_estimate",
		"compute_jackknife_bias_estimate",
		"compute_jackknife_std_error",
		"compute_jackknife_wald_two_sided_pval",
		"compute_jackknife_wald_confidence_interval"
	),
	wald = c(
		"compute_asymp_confidence_interval",
		"compute_asymp_two_sided_pval",
		"compute_wald_two_sided_pval",
		"compute_wald_confidence_interval"
	),
	likelihood_tests = c(
		"set_testing_type",
		"set_information_preference",
		"get_testing_type",
		"get_information_preference",
		"get_information_source_used",
		"get_supported_testing_types",
		"get_supported_information_preferences",
		"compute_score_two_sided_pval",
		"compute_score_confidence_interval",
		"compute_lik_ratio_two_sided_pval",
		"compute_lik_ratio_confidence_interval",
		"compute_gradient_two_sided_pval",
		"compute_gradient_confidence_interval"
	),
	likelihood_ratio = c(
		"compute_lik_ratio_two_sided_pval",
		"compute_lik_ratio_confidence_interval"
	),
	estimating_equation_likelihood_ratio = c(
		"compute_lik_ratio_two_sided_pval",
		"compute_lik_ratio_confidence_interval"
	),
	parametric_likelihood_bootstrap = c(
		"compute_lik_ratio_bootstrap_two_sided_pval",
		"compute_lik_ratio_bootstrap_confidence_interval",
		"get_last_param_bootstrap_diagnostics"
	),
	bayesian_bootstrap = c(
		"approximate_bayesian_bootstrap_distribution_beta_hat_T",
		"compute_bayesian_bootstrap_two_sided_pval",
		"compute_bayesian_bootstrap_confidence_interval"
	),
	nonparametric_bootstrap = c(
		"approximate_bootstrap_distribution_beta_hat_T",
		"compute_bootstrap_two_sided_pval",
		"compute_bootstrap_confidence_interval"
	)
)

EDI_COMPONENT_SPECS = list(
	RandomizationTest = list(
		status = "active",
		source_name = "InferenceRand",
		file = "inference_all_abstract_rand.R",
		dependencies = character(),
		owns_state = c(
			"custom_randomization_statistic_function", "compiled_cpp_stat_fn",
			"compiled_cpp_stat_src", "randomization_mc_control"
		),
		provides_capabilities = "randomization_test",
		allowed_likelihood_tiers = EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS,
		declare_body_references_optional = TRUE
	),
	RandomizationCI = list(
		status = "active",
		source_name = "InferenceRandCI",
		file = "inference_all_abstract_rand_ci.R",
		dependencies = "RandomizationTest",
		provides_capabilities = "randomization_ci",
		allowed_likelihood_tiers = EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS,
		declare_body_references_optional = TRUE
	),
	NonparametricBootstrap = list(
		status = "active",
		source_name = "InferenceNonParamBootstrap",
		file = "inference_all_abstract_non_param_boot.R",
		dependencies = "RandomizationCI",
		owns_state = c(
			"boot_distr_cache", "jack_distr_cache",
			"bootstrap_extreme_estimate_threshold",
			"bootstrap_extreme_ci_width_threshold"
		),
		provides_capabilities = "nonparametric_bootstrap",
		allowed_likelihood_tiers = EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS,
		declare_body_references_optional = TRUE
	),
	RandomizationBootstrap = list(
		status = "active",
		source_name = "InferenceRandBootstrap",
		file = "inference_all_abstract_rand_bootstrap.R",
		dependencies = "NonparametricBootstrap",
		owns_state = c("rand_boot_draws_counter", "brt_mc_control"),
		provides_capabilities = "randomization_bootstrap",
		allowed_likelihood_tiers = EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS,
		declare_body_references_optional = TRUE
	),
	BayesianBootstrap = list(
		status = "active",
		source_name = "InferenceBayesianBootstrap",
		file = "inference_all_abstract_bayesian_bootstrap.R",
		dependencies = "RandomizationBootstrap",
		owns_state = c(
			"current_bayesian_bootstrap_context",
			"current_bayesian_bootstrap_subject_or_block_weights"
		),
		provides_capabilities = "bayesian_bootstrap",
		allowed_likelihood_tiers = EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS,
		declare_body_references_optional = TRUE
	),
	Jackknife = list(
		status = "active",
		source_name = "InferenceJackknife",
		file = "inference_all_abstract_jackknife.R",
		dependencies = character(),
		provides_capabilities = "jackknife",
		allowed_likelihood_tiers = EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS,
		declare_body_references_optional = TRUE
	),
	Wald = list(
		status = "active",
		source_name = "InferenceAsymp",
		file = "inference_all_abstract_asymp.R",
		dependencies = "Jackknife",
		owns_state = "cached_mod",
		provides_capabilities = "wald",
		allowed_likelihood_tiers = EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS,
		declare_body_references_optional = TRUE
	),
	ExactTest = list(
		status = "active",
		source_name = "InferenceExact",
		file = "inference_all_abstract_exact.R",
		dependencies = character(),
		owns_state = "default_exact_type",
		provides_capabilities = "exact_test",
		allowed_likelihood_tiers = "none",
		declare_body_references_optional = TRUE
	),
	ExactBinomialIncidence = list(
		status = "active",
		source_name = "InferenceIncidExactBinomial",
		file = "inference_incidence_exact_binomial.R",
		dependencies = "ExactTest",
		requires_capabilities = "exact_test",
		provides_capabilities = "exact_binomial_incidence",
		allowed_likelihood_tiers = "none",
		declare_body_references_optional = TRUE
	),
	ExactFisherIncidence = list(
		status = "active",
		source_name = "InferenceIncidExactFisher",
		file = "inference_indicidence_exact_fisher.R",
		dependencies = "ExactTest",
		requires_capabilities = "exact_test",
		provides_capabilities = "exact_fisher_incidence",
		allowed_likelihood_tiers = "none",
		declare_body_references_optional = TRUE
	),
	ExactZhangIncidence = list(
		status = "active",
		source_name = "InferenceIncidenceExactZhang",
		file = "inference_incidence_exact_zhang.R",
		dependencies = "ExactTest",
		requires_capabilities = "exact_test",
		provides_capabilities = "exact_zhang_incidence",
		allowed_likelihood_tiers = "none",
		declare_body_references_optional = TRUE
	),
	LikelihoodTests = list(
		status = "active",
		source_name = "InferenceAsympLik",
		file = "inference_all_abstract_asymp_lik.R",
		dependencies = "Wald",
		owns_state = c(
			"likelihood_ci_max_abs", "testing_type",
			"information_preference", "information_source_used"
		),
		provides_capabilities = "likelihood_tests",
		allowed_likelihood_tiers = c("quasi", "partial", "full"),
		declare_body_references_optional = TRUE
	),
	ParametricLikelihoodBootstrap = list(
		status = "active",
		source_name = "InferenceParamBootstrap",
		file = "inference_all_abstract_param_boot.R",
		dependencies = "LikelihoodTests",
		owns_state = c(
			"bartlett_factor_mc_min_usable_fraction",
			"bartlett_factor_mc_max_attempts_per_replicate",
			"param_bootstrap_extreme_lr_threshold"
		),
		provides_capabilities = "parametric_likelihood_bootstrap",
		allowed_likelihood_tiers = c("partial", "full"),
		declare_body_references_optional = TRUE
	),
	StandardModelCache = list(
		status = "active",
		source_name = "InferenceAsympLikStdModCache",
		file = "inference_all_abstract_asymp_lik_std_mod_cache.R",
		dependencies = "LikelihoodTests",
		provides_capabilities = "standard_model_cache",
		allowed_likelihood_tiers = c("quasi", "partial", "full"),
		declare_body_references_optional = TRUE
	),
	CountLikelihoodPlumbing = list(
		status = "active",
		source_name = "InferenceCountLikelihood",
		file = "inference_all_abstract_count_likelihood.R",
		dependencies = "LikelihoodTests",
		provides_capabilities = "count_likelihood_plumbing",
		allowed_likelihood_tiers = c("quasi", "full"),
		declare_body_references_optional = TRUE
	),
	KKGEE = list(
		status = "active",
		source_name = "InferenceMixinKKGEEShared",
		file = "inference_mixin_kk_gee_shared.R",
		dependencies = character(),
		owns_state = c("m", "use_rcpp", "max_abs_reasonable_coef", "kk_gee_engine"),
		requires_state = c("any_censoring", "cached_values", "harden", "n", "y"),
		requires_public_methods = character(),
		requires_private_methods = c(
			"cache_nonestimable_estimate", "cache_nonestimable_se", "clear_nonestimable_state",
			"compute_z_or_t_ci_from_s_and_df", "compute_z_or_t_two_sided_pval_from_s_and_df",
			"create_design_matrix", "expand_subject_or_block_weights_to_row_weights",
			"fit_with_hardened_qr_column_dropping", "gee_family", "gee_response_type",
			"get_fit_warm_start_fisher", "get_fit_warm_start_for_length",
			"get_fit_warm_start_weights", "set_fit_warm_start", "shared_gee_dispatch"
		),
		optional_public_methods = c(
			"compute_jackknife_wald_confidence_interval",
			"compute_jackknife_wald_two_sided_pval"
		),
		optional_private_methods = character(),
		provides_capabilities = c("kk_gee", "wald"),
		allowed_likelihood_tiers = c("quasi"),
		conflicts = character()
	),
	KKGLMM = list(
		status = "active",
		source_name = "InferenceMixinKKGLMMShared",
		file = "inference_mixin_kk_glmm_shared.R",
		dependencies = character(),
		owns_state = c(
			"m", "optimization_alg", "skip_glmm_pkg_check",
			"max_abs_reasonable_coef", "kk_glmm_engine"
		),
		requires_state = c("any_censoring", "cached_values", "harden", "n", "y"),
		requires_public_methods = c("get_testing_type", "num_cores"),
		requires_private_methods = c(
			"cache_nonestimable_estimate", "clear_nonestimable_state",
			"compute_standard_error_from_information_matrix",
			"compute_z_or_t_ci_from_s_and_df", "compute_z_or_t_two_sided_pval_from_s_and_df",
			"create_design_matrix", "expand_subject_or_block_weights_to_row_weights",
			"fit_with_hardened_qr_column_dropping", "glmm_family", "glmm_response_type",
			"shared"
		),
		optional_public_methods = character(),
		optional_private_methods = character(),
		requires_super_methods = c("compute_asymp_confidence_interval", "compute_asymp_two_sided_pval"),
		provides_capabilities = c("kk_glmm", "wald"),
		allowed_likelihood_tiers = c("partial", "full"),
		conflicts = character()
	),
	KKPassThrough = list(
		status = "active",
		source_name = "InferenceMixinKKPassThrough",
		file = "inference_mixin_kk_passthrough.R",
		dependencies = character(),
		owns_state = c(
			"m", "kk_passthrough", "y_temp", "dead", "w", "X",
			"any_censoring", "best_par", "optimization_alg", "cached_mod",
			"best_X_colnames", "best_Xmm_colnames"
		),
		requires_state = c("cached_values", "des_obj_priv_int", "has_match_structure", "n", "y"),
		requires_public_methods = c("duplicate", "num_cores"),
		requires_private_methods = c(
			"assert_valid_bootstrap_type", "cache_nonestimable_estimate",
			"effective_parallel_cores",
			"expand_subject_or_block_weights_to_row_weights", "get_X", "has_private_method",
			"object_has_private_method", "par_lapply", "supports_likelihood_tests"
		),
		optional_public_methods = character(),
		optional_private_methods = c("compute_fast_bootstrap_distr", "compute_weighted_estimate_ivwc"),
		requires_super_methods = "approximate_bootstrap_distribution_beta_hat_T",
		provides_capabilities = c("kk_passthrough", "nonparametric_bootstrap"),
		allowed_likelihood_tiers = EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS,
		conflicts = character()
	),
	KKCompound = list(
		status = "active",
		source_name = "InferenceMixinKKPassThroughCompound",
		file = "inference_mixin_kk_passthrough_compound.R",
		dependencies = "KKPassThrough",
		owns_state = "kk_passthrough_compound",
		requires_state = c("cached_values", "has_match_structure"),
		requires_public_methods = character(),
		requires_private_methods = c("cache_nonestimable_estimate", "compute_basic_kk_match_data_impl"),
		optional_public_methods = character(),
		optional_private_methods = character(),
		provides_capabilities = "kk_compound",
		allowed_likelihood_tiers = EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS,
		conflicts = character()
	),
	OffOptimumLikelihoodEval = list(
		status = "active",
		source_name = "InferenceMixinOffOptimumLikelihoodEval",
		file = "inference_mixin_off_optimum_likelihood_eval.R",
		dependencies = character(),
		owns_state = character(),
		requires_state = character(),
		requires_public_methods = character(),
		requires_private_methods = c("get_default_information_source", "get_likelihood_test_spec"),
		optional_public_methods = character(),
		optional_private_methods = character(),
		provides_capabilities = "off_optimum_likelihood_eval",
		allowed_likelihood_tiers = c("partial", "full"),
		conflicts = character()
	),
	BartlettApproximation = list(
		status = "active",
		source_name = "InferenceExtBartlettApprox",
		file = "inference_ext_bartlett_approx.R",
		dependencies = "ParametricLikelihoodBootstrap",
		owns_state = c(
			"bartlett_factor_mc_min_usable_fraction",
			"bartlett_factor_mc_max_attempts_per_replicate"
		),
		requires_state = "active_resampling_operation",
		requires_public_methods = character(),
		requires_private_methods = c(
			"supports_lik_ratio_param_bootstrap",
			"run_param_bootstrap_replicates",
			"param_bootstrap_lr_extreme"
		),
		optional_public_methods = character(),
		optional_private_methods = character(),
		requires_super_methods = character(),
		requires_capabilities = "parametric_likelihood_bootstrap",
		provides_capabilities = "bartlett_approximation",
		allowed_likelihood_tiers = c("partial", "full"),
		conflicts = character()
	)
)

InferenceComponent = function(
	name,
	status = c("active", "scaffold"),
	source_name = name,
	file,
	public = list(),
	private = list(),
	dependencies = character(),
	provides_public_methods = names(public) %||% character(),
	provides_private_methods = names(private) %||% character(),
	owns_state = character(),
	requires_state = character(),
	requires_public_methods = character(),
	requires_private_methods = character(),
	optional_public_methods = character(),
	optional_private_methods = character(),
	requires_super_methods = character(),
	requires_capabilities = character(),
	provides_capabilities = character(),
	allowed_likelihood_tiers = EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS,
	conflicts = character(),
	allowed_host_overrides = list(public = character(), private = character()),
	forbidden_refs = list(private = character(), self = character(), super = character())
) {
	status = match.arg(status)
	component = list(
		name = name,
		status = status,
		source_name = source_name,
		file = file,
		public = public,
		private = private,
		dependencies = dependencies,
		provides_public_methods = provides_public_methods,
		provides_private_methods = provides_private_methods,
		owns_state = owns_state,
		requires_state = requires_state,
		requires_public_methods = requires_public_methods,
		requires_private_methods = requires_private_methods,
		optional_public_methods = optional_public_methods,
		optional_private_methods = optional_private_methods,
		requires_super_methods = requires_super_methods,
		requires_capabilities = requires_capabilities,
		provides_capabilities = provides_capabilities,
		allowed_likelihood_tiers = allowed_likelihood_tiers,
		conflicts = conflicts,
		allowed_host_overrides = allowed_host_overrides,
		forbidden_refs = forbidden_refs
	)
	validate_inference_component(component)
	component
}

validate_inference_component = function(component) {
	required = c(
		"name", "status", "source_name", "file", "public", "private",
		"dependencies", "provides_public_methods", "provides_private_methods",
		"owns_state", "requires_state", "requires_public_methods",
		"requires_private_methods", "optional_public_methods",
		"optional_private_methods", "requires_super_methods",
		"requires_capabilities", "provides_capabilities",
		"allowed_likelihood_tiers", "conflicts", "allowed_host_overrides",
		"forbidden_refs"
	)
	missing = setdiff(required, names(component))
	if (length(missing) > 0L) {
		stop(sprintf(
			"Inference component %s is missing required field(s): %s",
			component$name %||% "<unknown>",
			paste(missing, collapse = ", ")
		), call. = FALSE)
	}
	if (!is.character(component$name) || length(component$name) != 1L || !nzchar(component$name)) {
		stop("Inference component field `name` must be a non-empty character scalar.", call. = FALSE)
	}
	if (!(component$status %in% EDI_COMPONENT_ALLOWED_STATUSES)) {
		stop(sprintf("Inference component %s has invalid status.", component$name), call. = FALSE)
	}
	if (!is.list(component$public) || !is.list(component$private)) {
		stop(sprintf("Inference component %s must provide public/private lists.", component$name), call. = FALSE)
	}
	for (field in setdiff(required, c("public", "private", "allowed_host_overrides", "forbidden_refs"))) {
		if (!is.character(component[[field]])) {
			stop(sprintf("Inference component %s has non-character `%s`.", component$name, field), call. = FALSE)
		}
	}
	if (!all(component$allowed_likelihood_tiers %in% EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS)) {
		stop(sprintf("Inference component %s has invalid likelihood tier.", component$name), call. = FALSE)
	}
	if (!identical(sort(component$provides_public_methods), sort(component_public_names(component)))) {
		stop(sprintf("Inference component %s has stale public method metadata.", component$name), call. = FALSE)
	}
	if (!identical(sort(component$provides_private_methods), sort(component_private_names(component)))) {
		stop(sprintf("Inference component %s has stale private method metadata.", component$name), call. = FALSE)
	}
	invisible(TRUE)
}

register_inference_component = function(component) {
	validate_inference_component(component)
	if (exists(component$name, envir = EDI_INFERENCE_COMPONENTS, inherits = FALSE)) {
		stop(sprintf("Inference component already registered for %s.", component$name), call. = FALSE)
	}
	assign(component$name, component, envir = EDI_INFERENCE_COMPONENTS)
	invisible(component)
}

clear_inference_component_registry = function() {
	rm(list = ls(EDI_INFERENCE_COMPONENTS), envir = EDI_INFERENCE_COMPONENTS)
	invisible(TRUE)
}

inference_component_registry_as_list = function() {
	mget(ls(EDI_INFERENCE_COMPONENTS), envir = EDI_INFERENCE_COMPONENTS, inherits = FALSE)
}

get_inference_component = function(name) {
	if (!exists(name, envir = EDI_INFERENCE_COMPONENTS, inherits = FALSE)) {
		stop(sprintf("No inference component registered for %s.", name), call. = FALSE)
	}
	get(name, envir = EDI_INFERENCE_COMPONENTS, inherits = FALSE)
}

component_public_names = function(component) {
	names(component$public) %||% character()
}

component_private_names = function(component) {
	names(component$private) %||% character()
}

inference_component_source_parts = function(source) {
	if (inherits(source, "R6ClassGenerator")) {
		public = as.list(source$public_methods)
		public$clone = NULL
		private = c(as.list(source$private_methods), as.list(source$private_fields))
		return(list(public = public, private = private))
	}
	if (is.list(source) && all(c("public", "private") %in% names(source))) {
		return(list(
			public = as.list(source$public),
			private = as.list(source$private)
		))
	}
	stop("Inference component source must be an R6 generator or public/private list.", call. = FALSE)
}

complete_component_reference_contract = function(component) {
	refs = component_body_references(component)
	declared = component_declared_reference_names(component)
	missing_private = setdiff(refs$private, declared$private)
	missing_self = setdiff(refs$self, declared$self)
	missing_super = setdiff(refs$super, declared$super)
	component$optional_private_methods = sort(unique(c(
		component$optional_private_methods,
		missing_private
	)))
	component$optional_public_methods = sort(unique(c(
		component$optional_public_methods,
		missing_self
	)))
	component$requires_super_methods = sort(unique(c(
		component$requires_super_methods,
		missing_super
	)))
	validate_inference_component(component)
	component
}

populate_inference_component_registry = function(
	ns = environment(populate_inference_component_registry),
	component_names = names(EDI_COMPONENT_SPECS)
) {
	clear_inference_component_registry()
	for (name in component_names) {
		spec = EDI_COMPONENT_SPECS[[name]]
		source_name = spec$source_name %||% name
		declare_body_references_optional = isTRUE(spec$declare_body_references_optional)
		spec$source_name = NULL
		spec$declare_body_references_optional = NULL
		source = get(source_name, envir = ns, inherits = TRUE)
		parts = inference_component_source_parts(source)
		component = do.call(InferenceComponent, c(
			list(
				name = name,
				source_name = source_name,
				public = parts$public,
				private = parts$private,
				provides_public_methods = names(parts$public) %||% character(),
				provides_private_methods = names(parts$private) %||% character()
			),
			spec
		))
		if (declare_body_references_optional) {
			component = complete_component_reference_contract(component)
		}
		register_inference_component(component)
	}
	invisible(EDI_INFERENCE_COMPONENTS)
}

component_body_references = function(component) {
	refs = list(private = character(), self = character(), super = character())
	collect_from_expr = function(expr) {
		if (!is.call(expr) && !is.expression(expr)) return(invisible(NULL))
		if (is.call(expr) &&
				identical(as.character(expr[[1L]]), "$") &&
				length(expr) >= 3L &&
				is.symbol(expr[[2L]])) {
			lhs = as.character(expr[[2L]])
			if (lhs %in% names(refs)) {
				refs[[lhs]] <<- c(refs[[lhs]], as.character(expr[[3L]])[1L])
			}
		}
		for (i in seq_along(expr)) {
			collect_from_expr(expr[[i]])
		}
		invisible(NULL)
	}
	for (slot in c("public", "private")) {
		for (fn in component[[slot]]) {
			if (is.function(fn)) collect_from_expr(body(fn))
		}
	}
	lapply(refs, function(x) sort(unique(x)))
}

component_declared_reference_names = function(component) {
	list(
		private = sort(unique(c(
			component$provides_private_methods,
			component$owns_state,
			component$requires_state,
			component$requires_private_methods,
			component$optional_private_methods,
			component$forbidden_refs$private %||% character()
		))),
		self = sort(unique(c(
			component$provides_public_methods,
			component$requires_public_methods,
			component$optional_public_methods,
			component$forbidden_refs$self %||% character()
		))),
		super = sort(unique(c(
			component$requires_super_methods,
			component$forbidden_refs$super %||% character()
		)))
	)
}

validate_component_body_references = function(component) {
	refs = component_body_references(component)
	declared = component_declared_reference_names(component)
	offenders = character()
	for (receiver in names(refs)) {
		missing = setdiff(refs[[receiver]], declared[[receiver]])
		if (length(missing) > 0L) {
			offenders = c(offenders, sprintf(
				"%s has undeclared %s reference(s): %s",
				component$name,
				receiver,
				paste(missing, collapse = ", ")
			))
		}
	}
	if (length(offenders) > 0L) {
		stop(paste(offenders, collapse = "\n"), call. = FALSE)
	}
	invisible(TRUE)
}

validate_no_scaffold_effective_components = function(class_names = ls(EDI_INFERENCE_CLASS_REGISTRY)) {
	scaffold_components = names(Filter(function(component) {
		identical(component$status, "scaffold")
	}, inference_component_registry_as_list()))
	for (class_name in class_names) {
		components = get_effective_components(class_name)
		used_scaffolds = intersect(components, scaffold_components)
		if (length(used_scaffolds) > 0L) {
			stop(sprintf(
				"%s uses scaffold component(s): %s",
				class_name,
				paste(used_scaffolds, collapse = ", ")
			), call. = FALSE)
		}
	}
	invisible(TRUE)
}

resolve_component_dependencies = function(component_names, satisfied_components = character()) {
	duplicated_direct = sort(unique(component_names[duplicated(component_names)]))
	if (length(duplicated_direct) > 0L) {
		stop(sprintf(
			"Duplicate direct component(s): %s",
			paste(duplicated_direct, collapse = ", ")
		), call. = FALSE)
	}
	component_registry_names = ls(EDI_INFERENCE_COMPONENTS)
	unknown_components = setdiff(component_names, component_registry_names)
	if (length(unknown_components) > 0L) {
		stop(sprintf(
			"Unknown component(s): %s",
			paste(unknown_components, collapse = ", ")
		), call. = FALSE)
	}
	scaffold_components = component_names[vapply(component_names, function(component_name) {
		identical(get_inference_component(component_name)$status, "scaffold")
	}, logical(1L))]
	if (length(scaffold_components) > 0L) {
		stop(sprintf(
			"Scaffold component(s) cannot be resolved: %s",
			paste(scaffold_components, collapse = ", ")
		), call. = FALSE)
	}

	direct_dependency_hits = character()
	for (component_name in component_names) {
		deps = get_inference_component(component_name)$dependencies
		direct_dependency_hits = c(direct_dependency_hits, intersect(component_names, deps))
	}
	if (length(direct_dependency_hits) > 0L) {
		stop(sprintf(
			"Direct component list duplicates transitive dependency component(s): %s",
			paste(sort(unique(direct_dependency_hits)), collapse = ", ")
		), call. = FALSE)
	}

	resolved = character()
	visiting = character()
	expanded_dependencies = character()

	visit = function(component_name, path = character()) {
		if (component_name %in% satisfied_components) return(invisible(NULL))
		if (component_name %in% resolved) return(invisible(NULL))
		if (component_name %in% visiting) {
			cycle = c(path, component_name)
			stop(sprintf(
				"Component dependency cycle detected: %s",
				paste(cycle, collapse = " -> ")
			), call. = FALSE)
		}
		if (!(component_name %in% component_registry_names)) {
			stop(sprintf("Unknown component(s): %s", component_name), call. = FALSE)
		}
		component = get_inference_component(component_name)
		if (identical(component$status, "scaffold")) {
			stop(sprintf(
				"Scaffold component(s) cannot be resolved: %s",
				component_name
			), call. = FALSE)
		}
		visiting <<- c(visiting, component_name)
		for (dep in component$dependencies) {
			expanded_dependencies <<- c(expanded_dependencies, dep)
			visit(dep, c(path, component_name))
		}
		visiting <<- setdiff(visiting, component_name)
		resolved <<- c(resolved, component_name)
		invisible(NULL)
	}

	for (component_name in component_names) {
		visit(component_name)
	}

	duplicated_dependencies = sort(unique(expanded_dependencies[
		duplicated(expanded_dependencies) &
			!(expanded_dependencies %in% satisfied_components)
	]))
	if (length(duplicated_dependencies) > 0L) {
		stop(sprintf(
			"Component dependency graph duplicates transitive component(s): %s",
			paste(duplicated_dependencies, collapse = ", ")
		), call. = FALSE)
	}

	resolved
}

normalise_inference_overrides = function(overrides = list()) {
	if (is.null(overrides)) overrides = list()
	utils::modifyList(
		list(public = character(), private = character(), state = character(), public_private = character()),
		overrides
	)
}

entry_kinds = function(entries) {
	if (length(entries) == 0L) return(character())
	stats::setNames(
		ifelse(vapply(entries, is.function, logical(1L)), "method", "state"),
		names(entries)
	)
}

combine_component_slot = function(target_name, component_names, slot, host_entries = list(), overrides = list(), resolve = TRUE) {
	overrides = normalise_inference_overrides(overrides)
	if (isTRUE(resolve)) {
		component_names = resolve_component_dependencies(component_names)
	}
	allowed = unique(c(overrides[[slot]], overrides$state))
	combined = list()
	combined_kinds = character()
	for (component_name in component_names) {
		component_entries = as.list(get_inference_component(component_name)[[slot]])
		incoming_names = names(component_entries) %||% character()
		incoming_kinds = entry_kinds(component_entries)
		collisions = intersect(names(combined), incoming_names)
		kind_collisions = collisions[combined_kinds[collisions] != incoming_kinds[collisions]]
		undeclared_kind_collisions = setdiff(kind_collisions, allowed)
		if (length(undeclared_kind_collisions) > 0L) {
			stop(sprintf(
				"%s has undeclared %s method/state collision(s): %s",
				target_name,
				slot,
				paste(undeclared_kind_collisions, collapse = ", ")
			), call. = FALSE)
		}
		undeclared = setdiff(collisions, allowed)
		if (length(undeclared) > 0L) {
			stop(sprintf(
				"%s has undeclared %s component collision(s): %s",
				target_name,
				slot,
				paste(undeclared, collapse = ", ")
			), call. = FALSE)
		}
		combined = utils::modifyList(combined, component_entries)
		combined_kinds = entry_kinds(combined)
	}
	host_names = names(host_entries) %||% character()
	host_collisions = intersect(names(combined), host_names)
	undeclared_host_collisions = setdiff(host_collisions, allowed)
	if (length(undeclared_host_collisions) > 0L) {
		stop(sprintf(
			"%s overrides component %s member(s) without declaration: %s",
			target_name,
			slot,
			paste(undeclared_host_collisions, collapse = ", ")
		), call. = FALSE)
	}
	host_kinds = entry_kinds(host_entries)
	kind_collisions = host_collisions[combined_kinds[host_collisions] != host_kinds[host_collisions]]
	undeclared_kind_collisions = setdiff(kind_collisions, allowed)
	if (length(undeclared_kind_collisions) > 0L) {
		stop(sprintf(
			"%s overrides component %s method/state member(s) without declaration: %s",
			target_name,
			slot,
			paste(undeclared_kind_collisions, collapse = ", ")
		), call. = FALSE)
	}
	utils::modifyList(combined, host_entries)
}

assemble_public = function(target_name, component_names = character(), public = list(), overrides = list(), resolve = TRUE) {
	combine_component_slot(
		target_name = target_name,
		component_names = component_names,
		slot = "public",
		host_entries = public,
		overrides = overrides,
		resolve = resolve
	)
}

assemble_private = function(target_name, component_names = character(), private = list(), overrides = list(), resolve = TRUE) {
	combine_component_slot(
		target_name = target_name,
		component_names = component_names,
		slot = "private",
		host_entries = private,
		overrides = overrides,
		resolve = resolve
	)
}

r6_inherited_public_names = function(inherit) {
	if (is.null(inherit)) return(character())
	names(inherit$public_methods) %||% character()
}

r6_inherited_private_names = function(inherit) {
	if (is.null(inherit)) return(character())
	names(inherit$private_methods) %||% character()
}

validate_inference_class_definition = function(
	classname,
	inherit = NULL,
	component_names = character(),
	public = list(),
	private = list(),
	metadata = list(),
	overrides = list(),
	public_methods_for_capability = NULL,
	capability_requires = NULL,
	resolve_components = TRUE
) {
	overrides = normalise_inference_overrides(overrides)
	if (is.null(public_methods_for_capability)) {
		public_methods_for_capability = get("public_methods_for_capability", envir = parent.env(environment()))
	}
	if (is.null(capability_requires)) {
		capability_requires = get("capability_requires", envir = parent.env(environment()))
	}
	if (isTRUE(resolve_components)) {
		component_names = resolve_component_dependencies(component_names)
	}
	components = lapply(component_names, get_inference_component)
	likelihood_tier = metadata$likelihood_tier %||% "none"
	capabilities = unique(c(
		as.character(unlist(lapply(components, `[[`, "provides_capabilities"), use.names = FALSE)),
		metadata$capabilities %||% character()
	))
	public_names = unique(c(r6_inherited_public_names(inherit), names(public) %||% character()))
	private_names = unique(c(r6_inherited_private_names(inherit), names(private) %||% character()))

	public_private_overlap = intersect(public_names, private_names)
	undeclared_public_private_overlap = setdiff(public_private_overlap, overrides$public_private)
	if (length(undeclared_public_private_overlap) > 0L) {
		stop(sprintf(
			"%s has undeclared public/private name duplication: %s",
			classname,
			paste(undeclared_public_private_overlap, collapse = ", ")
		), call. = FALSE)
	}

	for (capability in intersect(capabilities, names(capability_requires))) {
		requirements = capability_requires[[capability]]
		allowed_tiers = requirements$likelihood_tier %||% EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS
		if (!(likelihood_tier %in% allowed_tiers)) {
			stop(sprintf(
				"%s advertises capability %s with disallowed likelihood tier `%s`.",
				classname,
				capability,
				likelihood_tier
			), call. = FALSE)
		}
		missing_public = setdiff(requirements$public_methods %||% character(), public_names)
		if (length(missing_public) > 0L) {
			stop(sprintf(
				"%s advertises capability %s without required public method(s): %s",
				classname,
				capability,
				paste(missing_public, collapse = ", ")
			), call. = FALSE)
		}
		missing_private = setdiff(requirements$private_methods %||% character(), private_names)
		if (length(missing_private) > 0L) {
			stop(sprintf(
				"%s advertises capability %s without required private method(s): %s",
				classname,
				capability,
				paste(missing_private, collapse = ", ")
			), call. = FALSE)
		}
		missing_state = setdiff(requirements$private_state %||% character(), private_names)
		if (length(missing_state) > 0L) {
			stop(sprintf(
				"%s advertises capability %s without required private state: %s",
				classname,
				capability,
				paste(missing_state, collapse = ", ")
			), call. = FALSE)
		}
		missing_capabilities = setdiff(requirements$capabilities %||% character(), capabilities)
		if (length(missing_capabilities) > 0L) {
			stop(sprintf(
				"%s advertises capability %s without required capability/capabilities: %s",
				classname,
				capability,
				paste(missing_capabilities, collapse = ", ")
			), call. = FALSE)
		}
		missing_metadata = setdiff(requirements$metadata %||% character(), names(metadata))
		if (length(missing_metadata) > 0L) {
			stop(sprintf(
				"%s advertises capability %s without required metadata field(s): %s",
				classname,
				capability,
				paste(missing_metadata, collapse = ", ")
			), call. = FALSE)
		}
	}

	for (component in components) {
		if (!(likelihood_tier %in% component$allowed_likelihood_tiers)) {
			stop(sprintf(
				"%s uses component %s with disallowed likelihood tier `%s`.",
				classname,
				component$name,
				likelihood_tier
			), call. = FALSE)
		}
		missing_public = setdiff(component$requires_public_methods, public_names)
		if (length(missing_public) > 0L) {
			stop(sprintf(
				"%s is missing public method(s) required by %s: %s",
				classname,
				component$name,
				paste(missing_public, collapse = ", ")
			), call. = FALSE)
		}
		missing_private = setdiff(component$requires_private_methods, private_names)
		if (length(missing_private) > 0L) {
			stop(sprintf(
				"%s is missing private method(s) required by %s: %s",
				classname,
				component$name,
				paste(missing_private, collapse = ", ")
			), call. = FALSE)
		}
		missing_state = setdiff(component$requires_state, private_names)
		if (length(missing_state) > 0L) {
			stop(sprintf(
				"%s is missing private state required by %s: %s",
				classname,
				component$name,
				paste(missing_state, collapse = ", ")
			), call. = FALSE)
		}
		missing_capabilities = setdiff(component$requires_capabilities, capabilities)
		if (length(missing_capabilities) > 0L) {
			stop(sprintf(
				"%s is missing capability required by %s: %s",
				classname,
				component$name,
				paste(missing_capabilities, collapse = ", ")
			), call. = FALSE)
		}
	}

	for (capability in intersect(names(public_methods_for_capability), capabilities)) {
		missing_methods = setdiff(public_methods_for_capability[[capability]], public_names)
		if (length(missing_methods) > 0L) {
			stop(sprintf(
				"%s advertises capability %s without public method(s): %s",
				classname,
				capability,
				paste(missing_methods, collapse = ", ")
			), call. = FALSE)
		}
	}
	invisible(TRUE)
}

define_inference_class = function(
	classname,
	inherit = NULL,
	components = character(),
	public = list(),
	private = list(),
	active = NULL,
	metadata = list(),
	overrides = list(),
	public_methods_for_capability = NULL,
	lock_objects = FALSE,
	...
) {
	if (!identical(lock_objects, FALSE)) {
		stop("define_inference_class() must keep `lock_objects = FALSE` until the R6 tree is stable.", call. = FALSE)
	}
	component_names = resolve_component_dependencies(components)
	assembled_public = assemble_public(classname, component_names, public, overrides, resolve = FALSE)
	assembled_private = assemble_private(classname, component_names, private, overrides, resolve = FALSE)
	validate_inference_class_definition(
		classname = classname,
		inherit = inherit,
		component_names = component_names,
		public = assembled_public,
		private = assembled_private,
		metadata = metadata,
		overrides = overrides,
		public_methods_for_capability = public_methods_for_capability,
		resolve_components = FALSE
	)
	R6::R6Class(
		classname = classname,
		lock_objects = FALSE,
		inherit = inherit,
		public = assembled_public,
		private = assembled_private,
		active = active,
		...
	)
}

assert_valid_mixin_composition = function(target_name, mixin_names, public_overrides = character(), private_overrides = character()) {
	duplicated_mixins = sort(unique(mixin_names[duplicated(mixin_names)]))
	if (length(duplicated_mixins) > 0L) {
		stop(sprintf(
			"%s composes duplicate mixin component(s): %s",
			target_name,
			paste(duplicated_mixins, collapse = ", ")
		), call. = FALSE)
	}
	for (mixin_name in mixin_names) {
		deps = EDI_MIXIN_DEPENDENCIES[[mixin_name]]
		missing_deps = setdiff(deps, mixin_names)
		if (length(missing_deps) > 0L) {
			stop(sprintf(
				"%s composes %s without required component(s): %s",
				target_name,
				mixin_name,
				paste(missing_deps, collapse = ", ")
			), call. = FALSE)
		}
	}
	for (slot in c("public", "private")) {
		slot_names = as.character(unlist(lapply(mixin_names, function(mixin_name) {
			mixin = get(mixin_name, envir = parent.frame(2L), inherits = TRUE)
			names(mixin[[slot]])
		}), use.names = FALSE))
		collisions = sort(unique(slot_names[duplicated(slot_names)]))
		allowed = EDI_MIXIN_ALLOWED_COLLISIONS[[target_name]][[slot]]
		undeclared = setdiff(collisions, allowed)
		if (length(undeclared) > 0L) {
			stop(sprintf(
				"%s has undeclared %s mixin collision(s): %s",
				target_name,
				slot,
				paste(undeclared, collapse = ", ")
			), call. = FALSE)
		}
		override_names = if (slot == "public") public_overrides else private_overrides
		override_collisions = intersect(slot_names, override_names)
		allowed_overrides = EDI_MIXIN_ALLOWED_OVERRIDES[[target_name]][[slot]]
		undeclared_overrides = setdiff(override_collisions, allowed_overrides)
		if (length(undeclared_overrides) > 0L) {
			stop(sprintf(
				"%s overrides inherited %s method(s) without declaration: %s",
				target_name,
				slot,
				paste(undeclared_overrides, collapse = ", ")
			), call. = FALSE)
		}
	}
	invisible(TRUE)
}

compose_inference_mixins = function(target_name, mixin_names, public = list(), private = list()) {
	assert_valid_mixin_composition(
		target_name = target_name,
		mixin_names = mixin_names,
		public_overrides = names(public),
		private_overrides = names(private)
	)
	component_names = unname(EDI_LEGACY_MIXIN_COMPONENT_NAMES[mixin_names])
	missing_component_names = mixin_names[is.na(component_names)]
	if (length(missing_component_names) > 0L) {
		stop(sprintf(
			"%s uses mixin(s) without canonical component mapping: %s",
			target_name,
			paste(missing_component_names, collapse = ", ")
		), call. = FALSE)
	}
	if (length(ls(EDI_INFERENCE_COMPONENTS)) == 0L) {
		populate_inference_component_registry(ns = parent.frame(), component_names = component_names)
	}
	allowed_collisions = EDI_MIXIN_ALLOWED_COLLISIONS[[target_name]] %||% list(public = character(), private = character())
	allowed_overrides = EDI_MIXIN_ALLOWED_OVERRIDES[[target_name]] %||% list(public = character(), private = character())
	overrides = list(
		public = unique(c(allowed_collisions$public, allowed_overrides$public, names(public))),
		private = unique(c(allowed_collisions$private, allowed_overrides$private, names(private)))
	)
	list(
		public = assemble_public(target_name, component_names, public, overrides = overrides, resolve = FALSE),
		private = assemble_private(target_name, component_names, private, overrides = overrides, resolve = FALSE)
	)
}
