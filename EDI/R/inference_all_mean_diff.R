#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring)
#' sequential experimental design estimation and test object
#' after the sequential design is completed.
#' For incidence responses coded as 0/1, the asymptotic inference is numerically
#' equivalent to Welch-Satterthwaite unpooled two-sample inference.
#'
#'
#' @export
#' @examples
#' \dontrun{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 6, response_type = "continuous")
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
#'
#' seq_des_inf = InferenceAllSimpleMeanDiff$
#'   new(seq_des)
#' seq_des_inf$compute_treatment_estimate()
#' seq_des_inf$compute_asymp_confidence_interval()
#' seq_des_inf$compute_asymp_two_sided_pval_for_treatment_effect()
#' }
InferenceAllSimpleMeanDiff = R6::R6Class("InferenceAllSimpleMeanDiff",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize a simple mean-difference inference object.
		#' @param des_obj A DesignSeqOneByOne object whose entire n subjects are assigned
		#'   and response y is recorded within.
		#' @param verbose Whether to print progress messages. Default \code{FALSE}.
		#' @param max_resample_attempts Maximum number of times a single bootstrap replicate
		#'   may be redrawn when the drawn sample fails validity screening. If all attempts
		#'   fail the replicate is recorded as \code{NA}, silently reducing the effective \code{B}.
		#'   Must be a positive integer. Default \code{50L}.
		initialize = function(des_obj, verbose = FALSE, max_resample_attempts = 50L){
			assertCount(max_resample_attempts, positive = TRUE)
			super$initialize(des_obj, verbose)
			private$max_resample_attempts = max_resample_attempts
		},

		#' @description
		#' Computes the appropriate estimate for mean difference
		#'
		#' @return	The setting-appropriate (see description) numeric estimate of the treatment effect
		#'
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			if (is.null(private$cached_values$beta_hat_T)){
				private$cached_values$yTs = private$y[private$w == 1]
				private$cached_values$yCs = private$y[private$w == 0]

				# Check for empty groups in bootstrap samples
				if (length(private$cached_values$yTs) == 0 || length(private$cached_values$yCs) == 0) {
					return(NA_real_) # Return NA if either group is empty
				}

				private$cached_values$beta_hat_T = mean(private$cached_values$yTs) - mean(private$cached_values$yCs)
			}
			private$cached_values$beta_hat_T
		},




		#' @description
		#' Computes a 1-alpha level frequentist confidence interval for the randomization test
		#'
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		#' @param	r		The number of randomization vectors. The default is 501.
		#' @param	pval_epsilon			The bisection algorithm tolerance. The default is 0.005.
		#' @param	show_progress		Show a text progress indicator.
		#' @param ci_search_control Optional randomization-CI search control list passed through 
		#'   to the base method.
		#' @return	A 1 - alpha sized frequentist confidence interval
		compute_confidence_interval_rand = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE, ci_search_control = NULL){
			if (private$des_obj_priv_int$response_type %in% c("proportion", "count", "survival")) {
				stop("Randomization confidence intervals are not supported for InferenceAllSimpleMeanDiff with proportion, count, or survival response types due to inconsistent estimator units on the transformed scale.")
			}
			super$compute_confidence_interval_rand(alpha = alpha, r = r, pval_epsilon = pval_epsilon, show_progress = show_progress, ci_search_control = ci_search_control)
		}
	),

	private = list(
		max_resample_attempts = 50L,
		get_standard_error = function(){
			if (is.null(private$cached_values$s_beta_hat_T)) private$shared()
			private$cached_values$s_beta_hat_T
		},

		get_degrees_of_freedom = function(){
			if (is.null(private$cached_values$df)) private$shared()
			private$cached_values$df
		},

		compute_fast_bootstrap_distr = function(B, ...) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			# KK designs use design-aware resampling not available via these args; fall back to R loop.
			if (private$is_KK) return(NULL)

			# Simple (non-KK) bootstrap: args = (n, y, dead, w)
			args = list(...)
			max_resample_attempts = private$max_resample_attempts
			n = args[[1]]
			y = args[[2]]
			dead = args[[3]]
			w = args[[4]]

			y_mat = matrix(NA_real_, nrow = n, ncol = B)
			w_mat = matrix(NA_integer_, nrow = n, ncol = B)

			for (b in 1:B) {
				attempt = 1
				repeat {
					i_b = sample(n, n, replace = TRUE)
					w_b = w[i_b]
					if (any(w_b == 1, na.rm = TRUE) && any(w_b == 0, na.rm = TRUE)) {
						if (!private$any_censoring) break
						dead_b_temp = dead[i_b]
						if (any(dead_b_temp[w_b == 1] == 1) && any(dead_b_temp[w_b == 0] == 1) && min(y[i_b]) > 0) break
					}
					attempt = attempt + 1
					if (attempt > max_resample_attempts) break
				}
				if (attempt <= max_resample_attempts) {
					y_mat[, b] = y[i_b]
					w_mat[, b] = w_b
				}
			}

			# Vectorized mean-diff per bootstrap sample: mean(y_b[w_b==1]) - mean(y_b[w_b==0])
			is_T = (w_mat == 1L)
			has_val = !is.na(w_mat)
			nT = colSums(is_T, na.rm = TRUE)
			nC = colSums(!is_T & has_val, na.rm = TRUE)
			sum_yT = colSums(y_mat * is_T, na.rm = TRUE)
			sum_yC = colSums(y_mat * (!is_T) * has_val, na.rm = TRUE)
			res = ifelse(nT > 0 & nC > 0, sum_yT / nT - sum_yC / nC, NA_real_)
			return(res)
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)

			# Optimization: w_mat is already pre-computed in generate_permutations
			w_mat = permutations$w_mat
			res = compute_simple_mean_diff_parallel_cpp(as.numeric(y), w_mat, as.numeric(delta), private$n_cpp_threads(ncol(w_mat)))
			return(res)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

					if (is.null(private$cached_values$beta_hat_T)){
						self$compute_treatment_estimate()
					}

	                # Check for insufficient samples for variance calculation
					nT = length(private$cached_values$yTs)
					nC = length(private$cached_values$yCs)

	                if (nT <= 1 || nC <= 1) { # Need at least 2 samples for variance
	                    private$cached_values$s_beta_hat_T = NA_real_
	                    private$cached_values$df = NA_real_
	                    private$cached_values$is_z = FALSE
	                    return() # Exit early
	                }

					s_1_sq = var(private$cached_values$yTs) / nT
					s_2_sq = var(private$cached_values$yCs) / nC
					private$cached_values$s_beta_hat_T = sqrt(s_1_sq + s_2_sq)
					private$cached_values$df = (s_1_sq + s_2_sq)^2 / (
													s_1_sq^2 / (nT - 1) + s_2_sq^2 / (nC - 1)
												) #Welch-Satterthwaite formula
					private$cached_values$is_z = FALSE
				}
			))
