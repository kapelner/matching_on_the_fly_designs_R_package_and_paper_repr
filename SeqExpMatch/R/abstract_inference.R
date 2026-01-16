#' Inference for A Sequential Design
#' 
#' @description
#' An abstract R6 Class that estimates, tests and provides intervals for a treatment effect in a sequential design.
#' This class takes a \code{SeqDesign} object as an input where this object
#' contains data for a fully completed sequential experiment (i.e. all treatment
#' assignments were allocated and all responses were collected). Then the user
#' specifies the type of estimation (mean_difference-or-medians or default_regression) and the type
#' of sampling assumption (i.e. the superpopulation assumption leading to MLE-or-KM-based inference or 
#' the finite population assumption implying randomization-exact-based inference) and then can query the
#' estimate and pval for the test. If the test is normal-theory based it is 
#' testing the population H_0: beta_T = 0 and if the test is a randomization test, 
#' it is testing the sharp null that H_0: Y_T_i = Y_C_i for all subjects. Confidence
#' interval construction is available for normal-theory based test type as well.
SeqDesignInference = R6::R6Class("SeqDesignInference",
	public = list(
		#' Begin Inference
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' 
		#' 
		#' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.

		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference 
		#' 							(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 							for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		#'
		#' @return A new `SeqDesignInference` object.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertClass(seq_des_obj, "SeqDesign")
			assertCount(num_cores, positive = TRUE)
			assertFlag(verbose)
			seq_des_obj$assert_experiment_completed()
			
			private$cached_values = list()
			private$any_censoring = seq_des_obj$any_censoring()
			private$seq_des_obj = seq_des_obj
			private$seq_des_obj_priv_int = seq_des_obj$.__enclos_env__$private
			private$y = private$seq_des_obj_priv_int$y
			private$y_temp = private$y
			private$w = private$seq_des_obj_priv_int$w
			private$dead = private$seq_des_obj_priv_int$dead
			private$is_KK = is(seq_des_obj, "SeqDesignKK14") #SeqDesignKK14 is the base class of all KK designs
			private$n = seq_des_obj$get_n()
			private$num_cores = num_cores
			private$verbose = verbose
			if (private$verbose){
				cat(paste0(
					"Initialized inference methods for a ",
					class(seq_des_obj)[1],
					" design and response type ",
					seq_des_obj$get_response_type(),
					".\n"
				))
			}
		},
		
		#' Set Custom Randomization Statistic Computation
		#'
		#' @description
		#' For advanced users only. This allows changing the default estimate inside randomization tests and interval construction.
		#' For example, when the response is continuous, instead of using ybarT - ybarC, you may want to use the studentized version
		#' i.e., (ybarT - ybarC) / sqrt(s^2_T / n_T + s^2_C / n_C). Work by Chung & Romano (2013, Annals of Statistics), 
		#' Janssen (1997), and others shows that studentized permutation tests are asymptotically valid and often asymptotically optimal.
		#' In finite samples, the studentized test often approximates a pivotal distribution better, leading to higher power.
		#' 
		#' 
		#' @param custom_randomization_statistic_function	A function that is run that returns a scalar value representing the statistic of interest
		#'													which is computed during each iteration sampling from the null distribution as w is drawn
		#'													drawn from the design. This function is embedded into this class and has write access to all of 
		#'													its data and functions (both public and private) so be careful! Setting this to NULL removes 
		#'													whatever function was set previously essentially. When there is no custom function, the default 
		#'													\code{self$compute_treatment_estimate()} will be run.
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
		#' #now let's set the statistic during randomization tests and intervals to the 
		#' #studentized average difference (the t stat). 
		#' seq_des_inf$set_custom_randomization_statistic_function(function(){
		#'    yTs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 1]
		#'	  yCs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 0]
		#'	  (mean(yTs) - mean(yCs)) / sqrt(var(yTs) / length(yTs) + var(yCs) / length(yCs))
		#' })
		#' }
		set_custom_randomization_statistic_function = function(custom_randomization_statistic_function){
			assertFunction(custom_randomization_statistic_function, null.ok = TRUE)			
			# Embed the function into this class as a private function
		    private[["custom_randomization_statistic_function"]] = custom_randomization_statistic_function
		    if (!is.null(custom_randomization_statistic_function)){
		   	 	# Make sure the function's environment is the class instance so it can access all the data
				environment(private[["custom_randomization_statistic_function"]]) = environment(self$initialize)
			}		    						
		},
		
		#' @description
		#' Creates the boostrap distribution of the estimate for the treatment effect
		#' 
		#' @param B						Number of bootstrap samples. The default is 501.
		#' 
		#' @return 	A vector of length \code{B} with the bootstrap values of the estimates of the treatment effect
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
		#' seq_des_inf = SeqDesignInference$new(seq_des)
		#' beta_hat_T_bs = seq_des_inf$approximate_bootstrap_distribution_beta_hat_T()
		#' ggplot(data.frame(beta_hat_T_bs = beta_hat_T_bs)) + 
		#'   geom_histogram(aes(x = beta_hat_T_bs))
		#' }
		#' 			
				approximate_bootstrap_distribution_beta_hat_T = function(B = 501){
					assertCount(B, positive = TRUE)
					
					n = private$seq_des_obj$get_n()
					y = private$seq_des_obj_priv_int$y
					dead = private$seq_des_obj_priv_int$dead
					w = private$seq_des_obj_priv_int$w
					X = private$get_X()
		
					duplicate_inference_fn <- function(){
						# Pass verbose flag to the duplicated object for logging internal to compute_treatment_estimate
						duplicated_obj = self$duplicate(verbose = private$verbose) 
                        if (private$verbose) {
                            cat("        duplicate_inference_fn: returned R6 object\n")
                        }
                        duplicated_obj
					}
		
					compute_estimate_from_bootstrap_data <- function(bootstrap_data){
						if (private$verbose) {
						  cat("      DEBUG (R): Entering compute_estimate_from_bootstrap_data\n")
						  cat("      DEBUG (R): bootstrap_data keys:", paste(names(bootstrap_data), collapse = ", "), "\n")
						}
						thread_inf_obj = bootstrap_data$inf_obj
		
						# Check if thread_inf_obj is NULL before accessing its private members
						if (is.null(thread_inf_obj)) {
							if (private$verbose) {
								cat("      DEBUG (R): thread_inf_obj is NULL. Returning NA_real_.\n")
							}
							return(NA_real_)
						}

						thread_inf_obj$.__enclos_env__$private$y = bootstrap_data$y
						thread_inf_obj$.__enclos_env__$private$dead = bootstrap_data$dead
						thread_inf_obj$.__enclos_env__$private$X = bootstrap_data$X
						thread_inf_obj$.__enclos_env__$private$w = bootstrap_data$w
						thread_inf_obj$.__enclos_env__$private$cached_values = list()
		
						result = tryCatch({
							if (private$verbose) {
								cat("      DEBUG (R): Calling compute_treatment_estimate() on duplicated object.\n")
							}
							thread_inf_obj$compute_treatment_estimate()
						}, error = function(e) {
							if (private$verbose) {
							  cat("      DEBUG (R): Bootstrap sample compute_treatment_estimate failed (R tryCatch):", e$message, "\n")
							}
							NA_real_ # Return NA for problematic bootstrap samples
						})
                        if (private$verbose) { # Use the verbose flag of the *original* object
                            cat("        DEBUG (R): compute_estimate_from_bootstrap_data: returning result:", result, "\n")
                        }
						result
					}
		
					if (private$verbose) {
					  cat("    Starting base_bootstrap_loop_cpp for B =", B, "\n")
					}
					indices = bootstrap_indices_cpp(n, B)
		
					beta_hat_T_bs = base_bootstrap_loop_cpp(
						indices,
						y,
						dead,
						X,
						w,
						duplicate_inference_fn,
						compute_estimate_from_bootstrap_data, # Corrected: use the actual function
						private$num_cores
					)
					if (private$verbose) {
					  cat("    Finished base_bootstrap_loop_cpp.\n")
					}
					beta_hat_T_bs		
				},
				#' @description
		#' Computes a 1-alpha level frequentist bootstrap confidence interval differently for all response types, estimate types and test types.
		#' 
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#' @param B						Number of bootstrap samples. The default is NA which corresponds to B=501.
		#' @param na.rm 				Should we remove beta_hat_T's that are NA's? Default is \code{FALSE}.
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
		#' seq_des_inf = SeqDesignInference$new(seq_des)
		#' seq_des_inf$compute_mle_confidence_interval()
		#' }
		#' 					
		compute_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, na.rm = FALSE){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertLogical(na.rm)
			quantile(private$get_or_cache_bootstrap_samples(B), c(alpha / 2, 1 - alpha / 2), na.rm = na.rm)
		},
			
		#' @description
		#' Computes a bootstrap two-sided p-value for H_0: betaT = delta. 
		#' It does so differently for all response types, estimate types and test types.
		#' 
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#' @param B						Number of bootstrap samples. The default is NA which corresponds to B=501.
		#' @param na.rm 				Should we remove beta_hat_T's that are NA's? Default is \code{FALSE}.
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
		#' seq_des_inf = SeqDesignInference$new(seq_des)
		#' seq_des_inf$compute_bootstrap_two_sided_pval()
		#' }
		#' 					
		compute_bootstrap_two_sided_pval = function(delta = 0, B = 501, na.rm = FALSE){
			assertNumeric(delta)
			assertLogical(na.rm)

			beta_hat_T_bs = private$get_or_cache_bootstrap_samples(B)
			
			2 * min(
				mean(delta < beta_hat_T_bs, na.rm = na.rm), 
				mean(delta > beta_hat_T_bs, na.rm = na.rm)
			)
		},				
		
		#' @description
		#' Under the sharp null of 
		#' forall i H_0: y_i_T - y_i_C = delta 
		#' there will be a distribution of the estimates of the treatment effect (over many realizations of assignments)
		#'
		#' @param nsim_exact_test		The number of randomization vectors to use. The default is 501.
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#' @param transform_responses	"none" for no transformation (default), "log" for log (your option when response type is survival) and 
		#'								"logit" for logit (your option when response type is proportion). This is mostly an 
		#'								internal parameter set to something besides "none" when computing randomization confidence intervals
		#'								for survival and proportion response types. 
		#' @return 	A vector of size \code{nsim_exact_test} that has the values of beta_hat_T over many w draws.
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
		#' seq_des_inf = SeqDesignInference$new(seq_des)
		#' beta_hat_T_diff_ws = seq_des_inf$compute_beta_hat_T_randomization_distr_under_sharp_null()
		#' ggplot(data.frame(beta_hat_T_diff_ws = beta_hat_T_diff_ws)) + 
		#'   geom_histogram(aes(x = beta_hat_T_diff_ws))
		#' }
		compute_beta_hat_T_randomization_distr_under_sharp_null = function(nsim_exact_test = 501, delta = 0, transform_responses = "none"){
			assertNumeric(delta)
			assertCount(nsim_exact_test, positive = TRUE)
			
			# Create a template design object with transformations applied
			# This will be duplicated by each thread
			seq_des_template = private$seq_des_obj$duplicate()
			
			# Apply response transformation if necessary
			if (transform_responses == "log"){
				seq_des_template$.__enclos_env__$private$y = log(copy(private$y))
			} else if (transform_responses == "logit"){
				seq_des_template$.__enclos_env__$private$y = logit(copy(private$y))
			}			
			
			# Apply delta adjustment if nonzero
			if (delta != 0){
				if (private$seq_des_obj_priv_int$response_type %in% c("count", "incidence")){
					stop("randomization tests with delta nonzero are not supported for count or incidence response types")
				}
				# Testing H_0: y_T_i - y_C_i = delta <=> (y_T_i - delta) - y_C_i = 0 
				# So adjust treatment responses downward by delta
				seq_des_template$.__enclos_env__$private$y[private$w == 1] = 
					seq_des_template$.__enclos_env__$private$y[private$w == 1] - delta
			}
			
			# Store is_KK flag for use in callback
			is_KK = private$is_KK
			
			# Iteration function that uses thread-local objects
			run_randomization_iteration <- function(thread_des_obj, thread_inf_obj){
				# Redraw treatment assignments
				thread_des_obj$.__enclos_env__$private$redraw_w_according_to_design()
				
				# Update inference object with new assignments
				thread_inf_obj$.__enclos_env__$private$seq_des_obj_priv_int = thread_des_obj$.__enclos_env__$private
				thread_inf_obj$.__enclos_env__$private$cached_values = list()
				
				# For KK designs, recompute match data
				if (is_KK && "compute_basic_match_data" %in% names(thread_inf_obj$.__enclos_env__$private)){
					thread_inf_obj$.__enclos_env__$private$compute_basic_match_data()
					if ("compute_reservoir_and_match_statistics" %in% names(thread_inf_obj$.__enclos_env__$private)){
						thread_inf_obj$.__enclos_env__$private$compute_reservoir_and_match_statistics()
					}
				}
				
				# Compute and return treatment estimate
				thread_inf_obj$.__enclos_env__$private$compute_treatment_estimate_during_randomization_inference()
			}

			# Pure R implementation for robustness
			# We use a single pair of objects if cores = 1
			# Parallelizing randomization tests in R can be done via mclapply if needed
			if (private$num_cores == 1) {
				thread_des_obj = seq_des_template$duplicate()
				thread_inf_obj = self$duplicate()
				
				beta_hat_T_diff_ws = replicate(nsim_exact_test, run_randomization_iteration(thread_des_obj, thread_inf_obj))
			} else {
				# Use parallel backend if more than 1 core requested
				beta_hat_T_diff_ws = unlist(parallel::mclapply(1:nsim_exact_test, function(i) {
					thread_des_obj = seq_des_template$duplicate()
					thread_inf_obj = self$duplicate()
					run_randomization_iteration(thread_des_obj, thread_inf_obj)
				}, mc.cores = private$num_cores))
			}
			
			beta_hat_T_diff_ws
		},
		
		#' @description
		#' Fisher's randomization test which means that H_0: y_i_T - y_i_C = delta for all subjects
		#' either the classic different-in-means estimate of the additive treatment effect, 
		#' i.e. ybar_T - ybar_C or the default_regression estimate of the additive treatment effect linearly i.e. 
		#' the treatment different adjusted linearly for the p covariates.
		#'
		#' @param nsim_exact_test		The number of randomization vectors to use in the randomization test (ignored if \code{test_type}
		#' 								is not "randomization-exact"). The default is 501 providing pvalue resolution to a fifth of a percent.
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#' @param transform_responses	"none" for no transformation (default), "log" for log and "logit" for logit. This is mostly an 
		#'								internal parameter set something besides "none" when computing randomization confidence intervals
		#'								for non-continuous responses.
		#' @param na.rm 				Should we remove beta_hat_T's that are NA's? Default is \code{FALSE}.
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
		#' seq_des_inf = SeqDesignInference$new(seq_des)
		#' seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
		#' }
		#' 		
		compute_two_sided_pval_for_treatment_effect_rand = function(nsim_exact_test = 501, delta = 0, transform_responses = "none", na.rm = FALSE){
			assertLogical(na.rm)
			#approximate the null distribution by computing estimates on many draws of w
			t0s = self$compute_beta_hat_T_randomization_distr_under_sharp_null(nsim_exact_test, delta, transform_responses)
			#this calculates the actual estimate to compare against the null distribution
			t = private$compute_treatment_estimate_during_randomization_inference()
			#finally compute the p-value
			if (private$verbose) {
			  cat("        Randomization Pval - t0s (length", length(t0s), ", NAs:", sum(is.na(t0s)), "): ", head(t0s), "\n")
			  cat("        Randomization Pval - t:", t, "\n")
			}

			if (na.rm){
				nsim_exact_test_adjusted = sum(!is.na(t0s))
			} else {
				nsim_exact_test_adjusted = length(t0s)
			}
            
            if (nsim_exact_test_adjusted == 0) return(NA_real_) # Avoid division by zero

			2 * min(
				sum(t0s >= t, na.rm = na.rm) / nsim_exact_test_adjusted, 
				sum(t0s <= t, na.rm = na.rm) / nsim_exact_test_adjusted
			)
		},

		
		#' @description
		#' Fisher's randomization test which means that H_0: y_i_T - y_i_C = delta for all subjects

		#' @description
		#' Computes a 1-alpha level frequentist confidence interval for the randomization test
		#' 
		#' Here we invert the randomization test that tests the strong null H_0: y_T_i - y_C_i = delta <=> (y_T_i - delta) - y_C_i = 0 so 
		#' we adjust the treatment responses downward by delta. We then find the set of all delta values that is above 1 - alpha/2 (i.e. two-sided)
		#' This is accomplished via a bisection algorithm (algorithm 1 of Glazer and Stark, 2025 available at
		#' https://arxiv.org/abs/2405.05238). These confidence intervals are exact to within tolerance \code{pval_epsilon}.
		#' As far as we know, this only works for response type continuous, uncensored survival (where we work in log-time and then
		#' return to natural time when finished) and proportion (where we work in logit-rate and return to rate when finished).
		#' 
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#' @param nsim_exact_test		The number of randomization vectors (applicable for test type "randomization-exact" only). 
		#' 								The default is 1000 providing good resolutions to confidence intervals.
		#' @param pval_epsilon			The bisection algorithm tolerance for the test inversion (applicable for test type "randomization-exact" only). 
		#' 								The default is to find a CI accurate to within 0.005.
		#' 
		#' @return 	A 1 - alpha sized frequentist confidence interval for the treatment effect
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
		#' seq_des_inf = SeqDesignInference$new(seq_des)
		#' seq_des_inf$compute_mle_confidence_interval()
		#' }
		#' 		
		compute_confidence_interval_rand = function(alpha = 0.05, nsim_exact_test = 501, pval_epsilon = 0.005){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertCount(nsim_exact_test, positive = TRUE)
			assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
						
			switch(private$seq_des_obj_priv_int$response_type,
				continuous = {
					lower_upper_ci_bounds = self$compute_bootstrap_confidence_interval(alpha / 100, na.rm = TRUE) #ensure a wider CI to be the starting position then pare down				
					
					ci = c(
							private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test, 
								l = lower_upper_ci_bounds[1], 
								u = self$compute_treatment_estimate(), 
								pval_th = alpha / 2, 
								tol = pval_epsilon, 
								transform_responses = "none",
								lower = TRUE),
							private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test, 
								l = self$compute_treatment_estimate(), 
								u = lower_upper_ci_bounds[2],  
								pval_th = alpha / 2, 
								tol = pval_epsilon, 
								transform_responses = "none", 
								lower = FALSE)
						)					
				},
				incidence =  stop("Confidence intervals are not supported for randomization tests for mean difference in incidence outomes"),
				count =      stop("Confidence intervals are not supported for randomization tests for mean difference in count outomes"),
				proportion = {
					lower_upper_ci_bounds = self$compute_bootstrap_confidence_interval(alpha / 100, na.rm = TRUE) #ensure a wider CI to be the starting position then pare down
					est = self$compute_treatment_estimate()
					
					#the bootstrap may give a bound outside of (0,1) which cannot be logitized, so correct it
					if (is.na(lower_upper_ci_bounds[1]) || lower_upper_ci_bounds[1] <= 0){
						lower_upper_ci_bounds[1] = min(.Machine$double.eps, est / 2)
					}
					if (is.na(lower_upper_ci_bounds[2]) || lower_upper_ci_bounds[2] >= 1){
						lower_upper_ci_bounds[2] = max(1 - .Machine$double.eps, (1 + est) / 2)
					}
					#now we work in logit rate so we can freely move around by +/-delta
					l_logit = logit(lower_upper_ci_bounds[1])
					u_logit = logit(lower_upper_ci_bounds[2])
					est_logit = logit(est)
					
					ci = c(
							private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test, 
								l = l_logit, 
								u = est_logit, 
								pval_th = alpha / 2, 
								tol = pval_epsilon,
								transform_responses = "logit", 
								lower = TRUE),
							private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test, 
								l = est_logit, 
								u = u_logit, 
								pval_th = alpha / 2, 
								tol = pval_epsilon, 
								transform_responses = "logit", 
								lower = FALSE)
						)					
					ci = inv_logit(ci)								
				},
				survival =   {
					lower_upper_ci_bounds = self$compute_bootstrap_confidence_interval(alpha / 100, na.rm = TRUE) #ensure a wider CI to be the starting position then pare down	
					est = self$compute_treatment_estimate()
					#the bootstrap may give a negative lower bound which cannot be logged, so correct it
					if (is.na(lower_upper_ci_bounds[1]) || lower_upper_ci_bounds[1] <= 0){
						lower_upper_ci_bounds[1] = min(.Machine$double.eps, est / 2)
					}
					#now we work in log time so we can freely move around by +/-delta
					l_log = log(lower_upper_ci_bounds[1])
					u_log = log(lower_upper_ci_bounds[2])
					est_log = log(est)
					
					ci = c(
							private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test, 
								l = l_log, 
								u = est_log, 
								pval_th = alpha / 2, 
								tol = pval_epsilon,
								transform_responses = "log", 
								lower = TRUE),
							private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test, 
								l = est_log, 
								u = u_log, 
								pval_th = alpha / 2, 
								tol = pval_epsilon, 
								transform_responses = "log",
								lower = FALSE)
						)
					#now we return to the natural time unit
					ci = exp(ci)					
				}			
			)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, sep = "%")
			ci
		},

		#' @description
		#' Duplicate this inference object
		#'
		#' @param verbose 	A flag indicating whether messages should be displayed to the user. Default is \code{FALSE}
		#' @return 			A new `SeqDesignInference` object with the same data
		duplicate = function(verbose = FALSE){
			i = tryCatch({
				do.call(get(class(self)[1])$new, args = list(
					seq_des_obj = private$seq_des_obj,
					num_cores = private$num_cores,
					verbose = verbose # Pass verbose here
				))
			}, error = function(e) {
				if (verbose) {
					cat("      Error during R6 object duplication:", e$message, "\n")
				}
				NULL # Return NULL on error
			})
			if (is.null(i)) return(NULL) # If construction failed, return NULL
			
			i$.__enclos_env__$private$seq_des_obj_priv_int = 					private$seq_des_obj_priv_int
			i$.__enclos_env__$private$any_censoring = 							private$any_censoring
			i$.__enclos_env__$private$n = 										private$n
			i$.__enclos_env__$private$p = 										private$p
			i$.__enclos_env__$private$X = 										private$X
			i$.__enclos_env__$private$custom_randomization_statistic_function = private$custom_randomization_statistic_function
			i
		}			
	),
	
	private = list(
		seq_des_obj = NULL,
		seq_des_obj_priv_int = NULL, 
		is_KK = NULL,
		any_censoring = NULL,
		num_cores = NULL,
		verbose = FALSE,
		n = NULL,
		p = NULL,
		y = NULL,
		w = NULL,
		dead = NULL,
		y_temp = NULL,
		X = NULL, #get_X is defined later as it needs some logic dependent on the design type
		custom_randomization_statistic_function = NULL,
		cached_values = list(),
		
		compute_treatment_estimate_during_randomization_inference = function(){
			if (is.null(private$custom_randomization_statistic_function)){ #i.e., the default
				self$compute_treatment_estimate()
			} else {
				private$custom_randomization_statistic_function()
			}
		},
				
		create_design_matrix = function(){
			cbind(1, private$seq_des_obj_priv_int$w, private$get_X())
		},
		
		get_X = function(){
			if (is.null(private$X)){
				if (is.null(private$seq_des_obj_priv_int$X)){
					private$seq_des_obj_priv_int$covariate_impute_if_necessary_and_then_create_model_matrix()
				}	
				private$X = private$seq_des_obj_priv_int$compute_all_subject_data()$X_all
			}
			private$X
		},	
		
		get_or_cache_bootstrap_samples = function(B){			
			if (is.null(private$cached_values$beta_hat_T_bs)){
				private$cached_values$beta_hat_T_bs = self$approximate_bootstrap_distribution_beta_hat_T(B)
				if (private$verbose) {
					cat("        Bootstrap samples summary (B=", B, "): \n")
					print(summary(private$cached_values$beta_hat_T_bs))
				}
			} else {
				B_0 = length(private$cached_values$beta_hat_T_bs)
				if (B_0 > B){
					return (private$cached_values$beta_hat_T_bs[1 : B]) #send back what we need but don't reduce the cache
				} else if (B_0 < B){ 
					new_samples = self$approximate_bootstrap_distribution_beta_hat_T(B - B_0)
					private$cached_values$beta_hat_T_bs = c(private$cached_values$beta_hat_T_bs, new_samples)
					if (private$verbose) {
					    cat("        Bootstrap samples summary (B=", B, ", added", length(new_samples), "): \n")
					    print(summary(private$cached_values$beta_hat_T_bs))
				    }
				}
			}
			private$cached_values$beta_hat_T_bs			
		},
		
		compute_z_or_t_ci_from_s_and_df = function(alpha){
			one_minus_alpha_over_two = 1 - alpha / 2
			z_or_t_val = 	if (private$cached_values$is_z){
								qnorm(one_minus_alpha_over_two)
							} else {
								qt(one_minus_alpha_over_two, private$cached_values$df)
							}
			moe = z_or_t_val * private$cached_values$s_beta_hat_T
			ci = private$cached_values$beta_hat_T + c(-moe, moe)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, sep = "%")
			ci
		},
		
		compute_z_or_t_two_sided_pval_from_s_and_df = function(delta){
			z_or_t_stat = (private$cached_values$beta_hat_T - delta) / private$cached_values$s_beta_hat_T
			z_or_t_stats = c(-z_or_t_stat, z_or_t_stat)
			probs = if (private$cached_values$is_z){ 
						stats::pnorm(z_or_t_stats)
					} else {
						stats::pt(z_or_t_stats, private$cached_values$df)
					}
			2 * min(probs)
		},			
		
		compute_ci_by_inverting_the_randomization_test_iteratively = function(nsim_exact_test, l, u, pval_th, tol, transform_responses, lower){
			# Pure R bisection loop for robustness with R6 objects and callbacks
			pval_l = as.numeric(self$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test, delta = l, transform_responses = transform_responses))
			pval_u = as.numeric(self$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test, delta = u, transform_responses = transform_responses))

			# Bisection loop
			repeat {
				# Check convergence
				if ((abs(u - l)) <= tol || (abs(pval_u - pval_l)) <= tol) {
					return(ifelse(lower, l, u))
				}

				# Compute midpoint
				m = (l + u) / 2.0
				pval_m = as.numeric(self$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test, delta = m, transform_responses = transform_responses))

				# Update bounds based on bisection logic
				if (pval_m >= pval_th && lower) {
					u = m
					pval_u = pval_m
				} else if (pval_m >= pval_th && !lower) {
					l = m
					pval_l = pval_m
				} else if (lower) {
					l = m
					pval_l = pval_m
				} else { # !lower
					u = m
					pval_u = pval_m
				}
			}
		}
	)
)
