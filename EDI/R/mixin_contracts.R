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
	composed_public = list()
	composed_private = list()
	for (mixin_name in mixin_names) {
		mixin = get(mixin_name, envir = parent.frame(), inherits = TRUE)
		composed_public = utils::modifyList(composed_public, as.list(mixin$public))
		composed_private = utils::modifyList(composed_private, as.list(mixin$private))
	}
	list(
		public = utils::modifyList(composed_public, public),
		private = utils::modifyList(composed_private, private)
	)
}
