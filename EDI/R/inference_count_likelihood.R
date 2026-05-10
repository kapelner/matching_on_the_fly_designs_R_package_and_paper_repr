#' Count Likelihood Inference Base
#'
#' Shared base for count models that are likelihood-backed and use the standard
#' treatment-effect caching / bootstrap / randomization-inference flow.
#'
#' @keywords internal
InferenceCountLikelihood = R6::R6Class("InferenceCountLikelihood",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = list(
		#' @description
		#' Computes the treatment estimate using the underlying model.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		}
	),

	private = list(
		supports_likelihood_tests = function(){
			TRUE
		},

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		generate_mod = function(estimate_only = FALSE) stop(class(self)[1], " must implement generate_mod()"),

		create_bootstrap_worker_state = function(){
			private$create_design_backed_bootstrap_worker_state()
		},

		load_bootstrap_sample_into_worker = function(worker_state, indices){
			private$load_bootstrap_sample_into_design_backed_worker(worker_state, indices)
		},

		compute_bootstrap_worker_estimate = function(worker_state){
			private$compute_bootstrap_worker_estimate_via_compute_treatment_estimate(worker_state)
		},

		get_standard_error = function(){
			private$shared(estimate_only = FALSE)
			if (isTRUE(private$supports_information_preference())) {
				se = private$compute_standard_error_from_information_matrix()
				if (is.finite(se)) return(se)
			}
			private$cached_values$s_beta_hat_T
		},

		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df
		},

		compute_score_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "score")
		},

		compute_gradient_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "gradient")
		},

		compute_lik_ratio_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "lik_ratio")
		},

		get_likelihood_test_spec = function(){
			NULL
		},

		make_warm_fit_null_wrapper = function(spec, cache_key){
			last_start = NULL
			fit_null_formals = tryCatch(names(formals(spec$fit_null)), error = function(e) character())
			accepts_start = "start" %in% fit_null_formals
			function(delta){
				warm_enabled = isTRUE(private$null_fit_warm_start_enabled)
				cache_state = if (warm_enabled) private$get_likelihood_null_warm_state(cache_key) else NULL
				start = if (warm_enabled) last_start else NULL
				if (warm_enabled && is.null(start) && !is.null(cache_state)) start = cache_state$start
				fit = tryCatch(
					if (accepts_start) spec$fit_null(delta, start = start) else spec$fit_null(delta),
					error = function(e) NULL
				)
				extract_start = spec$extract_start %||% function(fit_obj) NULL
				last_start <<- if (warm_enabled && accepts_start) tryCatch(extract_start(fit), error = function(e) NULL) else NULL
				if (warm_enabled && accepts_start) {
					private$set_likelihood_null_warm_state(cache_key, delta = delta, start = last_start)
				}
				fit
			}
		},

		compute_likelihood_test_two_sided_pval = function(delta, testing_type){
			private$get_memoized_likelihood_test_pval(
				delta = delta,
				testing_type = testing_type,
				spec = private$get_likelihood_test_spec(),
				warm_cache_key = paste0("likelihood_test:", testing_type)
			)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			has_cached_se = !is.null(private$cached_values$s_beta_hat_T) &&
				length(private$cached_values$s_beta_hat_T) == 1L &&
				isTRUE(is.finite(private$cached_values$s_beta_hat_T))
			if (isTRUE(!is.null(private$cached_values$beta_hat_T) && (estimate_only || has_cached_se))) return(invisible(NULL))
			model_output = private$generate_mod(estimate_only = estimate_only)
			private$cached_mod = model_output
			if (is.null(model_output)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}
			private$cached_values$beta_hat_T = model_output$b[2]
			if (estimate_only) return(invisible(NULL))

			ssq = model_output$ssq_b_2
			if (!is.null(ssq) && !is.na(ssq) && ssq > 0) {
				private$cached_values$s_beta_hat_T = sqrt(ssq)
			} else {
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$df = model_output$df %||% NA_real_
		}
	)
)
