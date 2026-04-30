#' GLM and Kaplan-Meier Inference
#'
#' Abstract class providing MLE/KM-based inference methods for GLM and survival models.
#'
#' @keywords internal
InferenceMLEorKMforGLMs = R6::R6Class("InferenceMLEorKMforGLMs",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
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

		compute_lik_ratio_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "lik_ratio")
		},

		get_likelihood_test_spec = function(){
			NULL
		},

		compute_likelihood_test_two_sided_pval = function(delta, testing_type){
			spec = private$get_likelihood_test_spec()
			if (is.null(spec)) {
				stop(class(self)[1], " does not expose a likelihood-test specification.", call. = FALSE)
			}

			j = as.integer(spec$j)
			if (length(j) != 1L || !is.finite(j) || j < 1L) return(NA_real_)

			null_fit = tryCatch(spec$fit_null(delta), error = function(e) NULL)
			null_params = null_fit$b %||% null_fit$params
			if (is.null(null_fit) || is.null(null_params) || length(null_params) < j || !is.finite(null_params[j])) {
				return(NA_real_)
			}

			if (testing_type == "score") {
				score = tryCatch(spec$score(null_fit), error = function(e) NULL)
				information = tryCatch(private$get_information_matrix(spec = spec, fit = null_fit), error = function(e) NULL)
				if (is.null(score) || is.null(information)) return(NA_real_)
				res = score_test_from_score_information_cpp(as.numeric(score), as.matrix(information), j)
				return(as.numeric(res$p_value %||% res))
			}

			if (testing_type == "lik_ratio") {
				full_negloglik = tryCatch(spec$neg_loglik(spec$full_fit), error = function(e) NA_real_)
				null_negloglik = tryCatch(spec$neg_loglik(null_fit), error = function(e) NA_real_)
				if (!is.finite(full_negloglik) || !is.finite(null_negloglik)) return(NA_real_)
				res = likelihood_ratio_test_from_negloglik_cpp(full_negloglik, null_negloglik, df = 1L)
				return(as.numeric(res$p_value %||% res))
			}

			stop("Unsupported testing_type: ", testing_type, call. = FALSE)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			has_cached_se = !is.null(private$cached_values$s_beta_hat_T) &&
				length(private$cached_values$s_beta_hat_T) == 1L &&
				isTRUE(is.finite(private$cached_values$s_beta_hat_T))
			if (isTRUE(!is.null(private$cached_values$beta_hat_T) && (estimate_only || has_cached_se))) return(invisible(NULL))
			model_output = private$generate_mod(estimate_only = estimate_only) #abstract function implemented by daughter classes. Should return a list with 'b' and 'ssq_b_2'.
			private$cached_mod = model_output
			# generate_mod may return NULL when the fit failed for all column subsets
			# (tryCatch in fit_with_hardened_qr_column_dropping catches convergence errors
			# and returns NULL). NULL$b[2] evaluates to NULL in R without error, which
			# would propagate as NULL from compute_estimate(). Guard here.
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
