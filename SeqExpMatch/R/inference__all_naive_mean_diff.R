#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @export
SeqDesignInferenceAllSimpleMeanDiff = R6::R6Class("SeqDesignInferenceAllSimpleMeanDiff",
	inherit = SeqDesignInference,
	public = list(
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference 
		#' 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs), 
		#' 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)	
			assertNoCensoring(private$any_censoring)
		},
		
		#' @description
		#' Computes the appropriate estimate for mean difference
		#' 
		#' @return 	The setting-appropriate (see description) numeric estimate of the treatment effect
		#' 
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignCRD$new(n = 6, response_type = "continuous")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceAllSimpleMeanDiff$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		#' 	
		compute_treatment_estimate = function(){
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
		#' Computes a 1-alpha level frequentist confidence interval differently for all response types, estimate types and test types.
		#' 
		#' Here we use the theory that MLE's computed for GLM's are asymptotically normal. 
		#' Hence these confidence intervals are asymptotically valid and thus approximate for any sample size.
		#' 
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#' 
		#' @return 	A (1 - alpha)-sized frequentist confidence interval for the treatment effect
		#' 
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignCRD$new(n = 6)
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceAllSimpleMeanDiff$new(seq_des)
		#' seq_des_inf$compute_mle_confidence_interval()
		#' }
		#'		
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)	
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$shared()
			}			
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		
		#' @description

		
		#' Computes a 2-sided p-value for all types of inferential settings written about in the initializer
		#'
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#' 
		#' @return 	The approximate frequentist p-value
		#' 
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignCRD$new(n = 6)
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceAllSimpleMeanDiff$new(seq_des)
		#' seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
		#' }
		#' 				
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			if (is.null(private$cached_values$df)){
				private$shared()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},

		#' @description
		#' Computes a 1-alpha level frequentist confidence interval for the randomization test
		#'
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#' @param nsim_exact_test		The number of randomization vectors. The default is 501.
		#' @param pval_epsilon			The bisection algorithm tolerance. The default is 0.005.
		#' @param show_progress		Show a text progress indicator.
		#' @return 	A 1 - alpha sized frequentist confidence interval
		compute_confidence_interval_rand = function(alpha = 0.05, nsim_exact_test = 501, pval_epsilon = 0.005, show_progress = TRUE){
			if (private$seq_des_obj_priv_int$response_type %in% c("proportion", "count", "survival")) {
				stop("Randomization confidence intervals are not supported for SeqDesignInferenceAllSimpleMeanDiff with proportion, count, or survival response types due to inconsistent estimator units on the transformed scale.")
			}
			super$compute_confidence_interval_rand(alpha = alpha, nsim_exact_test = nsim_exact_test, pval_epsilon = pval_epsilon, show_progress = show_progress)
		}
	),
	
			private = list(
		compute_fast_bootstrap_distr = function(B, max_resample_attempts, n, y, dead, w) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)

			y_mat = matrix(0.0, nrow = n, ncol = B)
			w_mat = matrix(0L, nrow = n, ncol = B)

			for (b in 1:B) {
				attempt = 1
				i_b = NULL
				w_b = NULL
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

				if (attempt > max_resample_attempts) {
					w_mat[, b] = rep(0L, n)
				} else {
					y_mat[, b] = y[i_b]
					w_mat[, b] = w_b
				}
			}

			res = compute_simple_mean_diff_parallel_cpp(y_mat, w_mat, private$num_cores)
			return(res)
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)

			nsim = length(permutations)
			n = length(y)
			
			w_mat = matrix(0L, nrow = n, ncol = nsim)
			for (i in 1:nsim) {
				w_mat[, i] = permutations[[i]]$w
			}

			y_shifted = copy(y)
			if (delta != 0) return(NULL)
			
			y_mat = matrix(y_shifted, nrow = n, ncol = nsim)

			res = compute_simple_mean_diff_parallel_cpp(y_mat, w_mat, private$num_cores)
			return(res)
		},

		shared = function(){
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
