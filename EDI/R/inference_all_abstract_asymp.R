#' Asymptotic Inference
#'
#' Abstract class for asymptotic inference.
#'
#' @keywords internal
InferenceAsymp = R6::R6Class("InferenceAsymp",
	lock_objects = FALSE,
	inherit = InferenceBoot,
	public = list(
		#' @description
		#' Computes an asymptotic confidence interval using the configured test.
		#'
		#' @param alpha					Significance level 1 - \code{alpha}. Default 0.05.
		#'
		#' @return 	A Wald-type confidence interval.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			switch(
				private$testing_type,
				wald = private$compute_wald_confidence_interval(alpha),
				score = private$invert_test_pval_confidence_interval(alpha),
				lik_ratio = private$invert_test_pval_confidence_interval(alpha),
				stop("Unsupported testing_type: ", private$testing_type)
			)
		},

		#' @description
		#' Computes an asymptotic two-sided p-value for the treatment effect.
		#'
		#' @param delta					Null treatment effect to test against. Default 0.
		#'
		#' @return 	The asymptotic p-value.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			switch(
				private$testing_type,
				wald = private$compute_wald_two_sided_pval_impl(delta),
				score = private$compute_score_two_sided_pval_impl(delta),
				lik_ratio = private$compute_lik_ratio_two_sided_pval_impl(delta),
				stop("Unsupported testing_type: ", private$testing_type)
			)
		},

		#' @description
		#' Sets the asymptotic testing method used by p-values and CIs.
		#'
		#' @param testing_type One of \code{"wald"}, \code{"score"}, or \code{"lik_ratio"}.
		#'
		#' @return The inference object, invisibly.
		set_testing_type = function(testing_type = c("wald", "score", "lik_ratio")){
			testing_type = private$normalize_testing_type(testing_type)
			supported = private$get_supported_testing_types_impl()
			if (!testing_type %in% supported) {
				stop(
					class(self)[1], " does not support testing_type = \"", testing_type,
					"\". Supported values are: ", paste(supported, collapse = ", "),
					call. = FALSE
				)
			}
			private$testing_type = testing_type
			invisible(self)
		},

		#' @description
		#' Gets the asymptotic testing method used by p-values and CIs.
		get_testing_type = function(){
			private$testing_type
		},

		#' @description
		#' Gets the asymptotic testing methods supported by this inference object.
		get_supported_testing_types = function(){
			private$get_supported_testing_types_impl()
		},

		#' @description
		#' Computes the Wald two-sided p-value regardless of configured testing type.
		compute_wald_two_sided_pval = function(delta = 0){
			private$compute_wald_two_sided_pval_impl(delta)
		},

		#' @description
		#' Computes the score two-sided p-value regardless of configured testing type.
		compute_score_two_sided_pval = function(delta = 0){
			private$compute_score_two_sided_pval_impl(delta)
		},

		#' @description
		#' Computes the likelihood-ratio two-sided p-value regardless of configured testing type.
		compute_lik_ratio_two_sided_pval = function(delta = 0){
			private$compute_lik_ratio_two_sided_pval_impl(delta)
		},

		#' @description
		#' Computes the likelihood-ratio two-sided p-value regardless of configured testing type.
		compute_likelihood_ratio_two_sided_pval = function(delta = 0){
			private$compute_lik_ratio_two_sided_pval_impl(delta)
		},

		#' @description
		#' Abstract method to compute the treatment estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		#' @return 	A scalar treatment estimate.
		compute_estimate = function(estimate_only = FALSE){
			stop("Must be implemented by concrete class.")
		},

		#' @description
		#' Returns the model object from the last call that produced the treatment
		#' estimate and SE. Calls \code{compute_estimate()} first if needed.
		#'
		#' @return The cached model object (type depends on the concrete class).
		get_mod = function(){
			if (is.null(private$cached_mod)) self$compute_estimate()
			private$cached_mod
		},

		#' @description
		#' Prints a summary of the model from the last call that produced the
		#' treatment estimate and SE.
		get_summary = function(){
			mod = self$get_mod()
			if (is.null(mod)) {
				cat("No model available (call compute_estimate() first).\n")
				return(invisible(NULL))
			}
			if (identical(class(mod), "list")) {
				if (!is.null(private$cached_values$summary_table)) {
					print(private$cached_values$summary_table)
				} else {
					print(mod)
				}
			} else {
				print(summary(mod))
			}
			invisible(NULL)
		}
	),

	private = list(
		cached_mod = NULL,
		testing_type = "wald",
		get_standard_error = function() stop("Must be implemented by concrete class or shared helper."),
		get_degrees_of_freedom = function() NA_real_,

		get_supported_testing_types_impl = function(){
			if (isTRUE(private$supports_likelihood_tests())) c("wald", "score", "lik_ratio") else "wald"
		},

		supports_likelihood_tests = function(){
			FALSE
		},

		normalize_testing_type = function(testing_type){
			if (length(testing_type) != 1L) testing_type = testing_type[1L]
			testing_type = tolower(as.character(testing_type))
			switch(
				testing_type,
				wald = "wald",
				score = "score",
				lr = "lik_ratio",
				lrt = "lik_ratio",
				lik_ratio = "lik_ratio",
				likelihood_ratio = "lik_ratio",
				stop("testing_type must be one of: wald, score, lik_ratio", call. = FALSE)
			)
		},

		compute_wald_confidence_interval = function(alpha){
			est = self$compute_estimate()
			se = private$get_standard_error()
			df = private$get_degrees_of_freedom()
			
			if (!is.finite(est) || !is.finite(se) || se <= 0) return(c(NA_real_, NA_real_))
			
			critical_val = if (is.finite(df)) stats::qt(1 - alpha / 2, df = df) else stats::qnorm(1 - alpha / 2)
			
			ci = c(est - critical_val * se, est + critical_val * se)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},

		compute_wald_two_sided_pval_impl = function(delta){
			est = self$compute_estimate()
			se = private$get_standard_error()
			df = private$get_degrees_of_freedom()
			
			if (!is.finite(est) || !is.finite(se) || se <= 0) return(NA_real_)
			
			t_stat = (est - delta) / se
			
			if (is.finite(df)) {
				2 * stats::pt(-abs(t_stat), df = df)
			} else {
				2 * stats::pnorm(-abs(t_stat))
			}
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
			if (is.null(null_fit) || is.null(null_fit$b) && is.null(null_fit$params)) {
				return(NA_real_)
			}

			if (testing_type == "score") {
				score = tryCatch(spec$score(null_fit), error = function(e) NULL)
				information = tryCatch({
					if (!is.null(spec$fisher_information)) {
						spec$fisher_information(null_fit)
					} else if (!is.null(null_fit$fisher_information)) {
						null_fit$fisher_information
					} else if (!is.null(spec$observed_information)) {
						spec$observed_information(null_fit)
					} else if (!is.null(null_fit$observed_information)) {
						null_fit$observed_information
					} else {
						spec$information(null_fit)
					}
				}, error = function(e) NULL)
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

		invert_test_pval_confidence_interval = function(alpha){
			est = self$compute_estimate()
			if (!is.finite(est)) return(c(NA_real_, NA_real_))

			target = function(delta){
				pval = self$compute_asymp_two_sided_pval(delta)
				if (!is.finite(pval)) return(NA_real_)
				pval - alpha
			}

			se = private$get_standard_error()
			step = if (is.finite(se) && se > 0) se else max(abs(est), 1)
			step = max(step, 1e-4)

			find_bound = function(direction){
				inner = tryCatch(target(est), error = function(e) NA_real_)
				if (!is.finite(inner)) return(NA_real_)
				if (inner < 0) return(est)

				lower = est
				upper = est
				outer = NA_real_
				for (i in seq_len(60L)){
					candidate = est + direction * step * 2^(i - 1L)
					outer = tryCatch(target(candidate), error = function(e) NA_real_)
					if (is.finite(outer) && outer <= 0) {
						if (direction < 0) {
							lower = candidate
							upper = est
						} else {
							lower = est
							upper = candidate
						}
						break
					}
				}
				if (!is.finite(outer) || outer > 0) return(NA_real_)
				tryCatch(stats::uniroot(target, lower = lower, upper = upper, tol = 1e-6)$root, error = function(e) NA_real_)
			}

			ci = c(find_bound(-1), find_bound(1))
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},

		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		create_bootstrap_worker_state = function(){
			private$create_design_backed_bootstrap_worker_state()
		},

		load_bootstrap_sample_into_worker = function(worker_state, indices){
			private$load_bootstrap_sample_into_design_backed_worker(worker_state, indices)
		},

		compute_bootstrap_worker_estimate = function(worker_state){
			private$compute_bootstrap_worker_estimate_via_compute_treatment_estimate(worker_state)
		},
		
		# Shared helpers for z/t tests
		compute_z_or_t_ci_from_s_and_df = function(alpha){
			beta_hat_T = private$cached_values$beta_hat_T
			s_beta_hat_T = private$cached_values$s_beta_hat_T
			is_z = private$cached_values$is_z
			df = private$cached_values$df
			
			if (length(beta_hat_T) != 1L || length(s_beta_hat_T) != 1L) return(c(NA_real_, NA_real_))
			if (!is.finite(beta_hat_T) || !is.finite(s_beta_hat_T) || s_beta_hat_T <= 0) return(c(NA_real_, NA_real_))
			
			mult = if (isTRUE(is_z) || !is.finite(df)) stats::qnorm(1 - alpha / 2) else stats::qt(1 - alpha / 2, df = df)
			ci = c(beta_hat_T - mult * s_beta_hat_T, beta_hat_T + mult * s_beta_hat_T)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},

		compute_z_or_t_two_sided_pval_from_s_and_df = function(delta){
			beta_hat_T = private$cached_values$beta_hat_T
			s_beta_hat_T = private$cached_values$s_beta_hat_T
			is_z = private$cached_values$is_z
			df = private$cached_values$df
			
			if (length(beta_hat_T) != 1L || length(s_beta_hat_T) != 1L) return(NA_real_)
			if (!is.finite(beta_hat_T) || !is.finite(s_beta_hat_T) || s_beta_hat_T <= 0) return(NA_real_)
			
			val = (beta_hat_T - delta) / s_beta_hat_T
			if (isTRUE(is_z) || !is.finite(df)) 2 * stats::pnorm(-abs(val)) else 2 * stats::pt(-abs(val), df = df)
		},

		invert_ci_to_find_two_sided_pval_for_treatment_effect = function(delta = 0){
			# Use bisection to find alpha such that delta is on the boundary of the CI.
			# The p-value is the largest alpha such that delta is outside the (1-alpha) CI.
			
			# Check if delta is already outside a very tight CI (alpha=0.99)
			# or inside a very wide CI (alpha=1e-10)
			
			f = function(alpha) {
				ci = self$compute_asymp_confidence_interval(alpha = alpha)
				if (any(!is.finite(ci))) return(NA_real_)
				# Distance to the nearest boundary. If delta is inside, return positive.
				# If delta is outside, return negative.
				# Actually, easier: return 0 if delta is on boundary.
				min(abs(ci - delta))
			}
			
			# More robust: check if delta is within the range of the estimate
			est = self$compute_estimate()
			if (!is.finite(est)) return(NA_real_)
			
			# We want to find alpha such that ci_boundary(alpha) == delta
			# This is monotonic in alpha.
			
			target_fn = function(alpha) {
				ci = self$compute_asymp_confidence_interval(alpha = alpha)
				if (est > delta) {
					# Null is below estimate, we care about the lower bound
					ci[1] - delta
				} else {
					# Null is above estimate, we care about the upper bound
					ci[2] - delta
				}
			}
			
			# p-value is usually between 0 and 1.
			lower = .Machine$double.eps
			upper = 1 - .Machine$double.eps
			
			# Check if target_fn has different signs at boundaries
			tl = tryCatch(target_fn(lower), error = function(e) NA_real_)
			tu = tryCatch(target_fn(upper), error = function(e) NA_real_)
			
			if (!is.finite(tl) || !is.finite(tu)) return(NA_real_)
			if (tl * tu > 0) {
				# If both are same sign, delta is either always outside or always inside
				if (abs(est - delta) < .Machine$double.eps) return(1)
				if (abs(tl) < abs(tu)) return(lower) else return(upper)
			}
			
			stats::uniroot(target_fn, lower = lower, upper = upper, tol = 1e-6)$root
		}
	)
)
