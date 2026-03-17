# Inference for A Sequential Design
#
# @description
# An abstract R6 Class that estimates, tests and provides intervals for a treatment effect in a sequential design.
# This class takes a \code{SeqDesign} object as an input where this object
# contains data for a fully completed sequential experiment (i.e. all treatment
# assignments were allocated and all responses were collected). Then the user
# specifies the type of estimation (mean_difference-or-medians or default_regression) and the type
# of sampling assumption (i.e. the superpopulation assumption leading to MLE-or-KM-based inference or
# the finite population assumption implying randomization-exact-based inference) and then can query the
# estimate and pval for the test. If the test is normal-theory based it is
# testing the population H_0: beta_T = 0 and if the test is a randomization test,
# it is testing the sharp null that H_0: Y_T_i = Y_C_i for all subjects. Confidence
# interval construction is available for normal-theory based test type as well.
#
# @keywords internal
SeqDesignInference = R6::R6Class("SeqDesignInference",
	public = list(
		# Begin Inference
		# @description
		# Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#
		#
		# @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.

		# @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		# 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference
		# 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs),
		# 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
		# @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		#
		# @return A new `SeqDesignInference` object.
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
			private$cached_values$rand_distr_cache = list()
			private$cached_values$permutations_cache = list()
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

		# Set Custom Randomization Statistic Computation
		#
		# @description
		# For advanced users only. This allows changing the default estimate inside randomization tests and interval construction.
		# For example, when the response is continuous, instead of using ybarT - ybarC, you may want to use the studentized version
		# i.e., (ybarT - ybarC) / sqrt(s^2_T / n_T + s^2_C / n_C). Work by Chung & Romano (2013, Annals of Statistics),
		# Janssen (1997), and others shows that studentized permutation tests are asymptotically valid and often asymptotically optimal.
		# In finite samples, the studentized test often approximates a pivotal distribution better, leading to higher power.
		#
		#
		# @param custom_randomization_statistic_function	A function that is run that returns a scalar value representing the statistic of interest
		#													which is computed during each iteration sampling from the null distribution as w is drawn
		#													drawn from the design. This function is embedded into this class and has write access to all of
		#													its data and functions (both public and private) so be careful! Setting this to NULL removes
		#													whatever function was set previously essentially. When there is no custom function, the default
		#													\code{self$compute_treatment_estimate()} will be run.
		# @examples
		# \dontrun{
		# seq_des = SeqDesignCRD$new(n = 6)
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		# seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#
		# seq_des_inf = SeqDesignInferenceAllSimpleMeanDiff$new(seq_des)
		# #now let's set the statistic during randomization tests and intervals to the
		# #studentized average difference (the t stat).
		# seq_des_inf$set_custom_randomization_statistic_function(function(){
		#    yTs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 1]
		#	  yCs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 0]
		#	  (mean(yTs) - mean(yCs)) / sqrt(var(yTs) / length(yTs) + var(yCs) / length(yCs))
		# })
		# }
		set_custom_randomization_statistic_function = function(custom_randomization_statistic_function){
			assertFunction(custom_randomization_statistic_function, null.ok = TRUE)
			# Embed the function into this class as a private function
		    private[["custom_randomization_statistic_function"]] = custom_randomization_statistic_function
			if (!is.null(custom_randomization_statistic_function)){
		   	 	# Make sure the function's environment is the class instance so it can access all the data
				environment(private[["custom_randomization_statistic_function"]]) = environment(self$initialize)
			}
			private$cached_values$t0s_rand = NULL
			private$cached_values$rand_distr_cache = list()
			private$cached_values$custom_stat_analysis = NULL
		},

		# @description
		# Computes the appropriate estimate
		# @return 	The numeric estimate of the treatment effect
		compute_treatment_estimate = function(){
			stop(class(self)[1], " must implement compute_treatment_estimate()")
		},

		# @description
		# Computes a 1-alpha level frequentist confidence interval
		# @param alpha					The confidence level.
		# @return 	A confidence interval
		compute_asymp_confidence_interval = function(alpha = 0.05){
			stop(class(self)[1], " must implement compute_asymp_confidence_interval(alpha)")
		},

		# @description
		# Computes a 2-sided p-value
		# @param delta					The null difference to test against.
		# @return 	The frequentist p-value
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			stop(class(self)[1], " must implement compute_asymp_two_sided_pval_for_treatment_effect(delta)")
		},

		# @description
		# Computes an exact confidence interval implementation-specific to the concrete class.
		# @param ... Parameters passed by subclasses.
		compute_exact_confidence_interval = function(...){
			stop("This method is not implemented in ", class(self)[1])
		},

		# @description
		# Computes an exact two-sided p-value implementation-specific to the concrete class.
		# @param ... Parameters passed by subclasses.
		compute_exact_two_sided_pval_for_treatment_effect = function(...){
			stop("This method is not implemented in ", class(self)[1])
		},

		# @description
		# Creates the boostrap distribution of the estimate for the treatment effect
		#
		# @param B						Number of bootstrap samples. The default is 501.
		# @param max_resample_attempts	Maximum number of times to redraw a bootstrap sample that fails
		# 								validity checks (both treatment arms present; for censored responses,
		# 								each arm must also contain at least one observed event and all survival
		# 								times must be strictly positive to avoid numerical underflow in the
		# 								Weibull log-likelihood). If all attempts are exhausted the sample is
		# 								recorded as \code{NA}. The default is 50.
		#
		# @return 	A vector of length \code{B} with the bootstrap values of the estimates of the treatment effect
		#
		# @examples
		# \dontrun{
		# seq_des = SeqDesignCRD$new(n = 6)
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		# seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#
		# seq_des_inf = SeqDesignInference$new(seq_des)
		# beta_hat_T_bs = seq_des_inf$approximate_bootstrap_distribution_beta_hat_T()
		# ggplot(data.frame(beta_hat_T_bs = beta_hat_T_bs)) +
		#   geom_histogram(aes(x = beta_hat_T_bs))
		# }
		#
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, max_resample_attempts = 50){
			assertCount(B, positive = TRUE)
			assertCount(max_resample_attempts, positive = TRUE)

			n = private$seq_des_obj$get_n()
			y = private$y
			dead = private$dead
			w = private$w
			X = private$get_X()

			# Check if subclass provides a C++ OpenMP dispatcher to bypass the slow R loop
			if (private$has_private_method("compute_fast_bootstrap_distr")) {
				fast_distr = private$compute_fast_bootstrap_distr(B, max_resample_attempts, n, y, dead, w)
				if (!is.null(fast_distr)) {
					return(fast_distr)
				}
			}

			# Pure R bootstrap implementation (fallback for GLMs or when C++ fast-path is unavailable)
			if (private$num_cores == 1) {
				beta_hat_T_bs = rep(NA_real_, B)
				for (b in 1:B) {
					# Generate bootstrap indices (sample with replacement)
					attempt = 1
					repeat {
						i_b = sample_int_replace_cpp(n, n)
						w_b = w[i_b]
						if (any(w_b == 1, na.rm = TRUE) && any(w_b == 0, na.rm = TRUE)) {
							if (!private$any_censoring) break
							dead_b_temp = dead[i_b]
							if (any(dead_b_temp[w_b == 1] == 1) && any(dead_b_temp[w_b == 0] == 1) &&
								min(y[i_b]) > 0) break
						}
						attempt = attempt + 1
						if (attempt > max_resample_attempts) break
					}

					# Create bootstrap sample
					y_b = y[i_b]
					dead_b = dead[i_b]
					X_b = X[i_b, , drop = FALSE]

					# Create a duplicate inference object for this bootstrap sample
					boot_inf_obj = self$duplicate()

					# Update the bootstrap object with resampled data
					boot_inf_obj$.__enclos_env__$private$y = y_b
					boot_inf_obj$.__enclos_env__$private$dead = dead_b
					boot_inf_obj$.__enclos_env__$private$w = w_b
					boot_inf_obj$.__enclos_env__$private$X = X_b
					boot_inf_obj$.__enclos_env__$private$cached_values = list()

					# Compute treatment estimate on bootstrap sample
					if (attempt > max_resample_attempts) {
						if (private$verbose) {
							cat("      Bootstrap sample", b, "failed: unable to draw both treatment groups after", max_resample_attempts, "attempts\n")
						}
						beta_hat_T_bs[b] = NA_real_
					} else {
						beta_hat_T_bs[b] = tryCatch({
							boot_inf_obj$compute_treatment_estimate()
						}, error = function(e) {
							if (private$verbose) {
								cat("      Bootstrap sample", b, "failed:", e$message, "\n")
							}
							NA_real_
						})
					}
				}
				return(beta_hat_T_bs)
			} else {
				# Parallel bootstrap execution for R loop fallback
				beta_hat_T_bs = unlist(parallel::mclapply(1:B, function(b) {
					attempt = 1
					repeat {
						i_b = sample_int_replace_cpp(n, n)
						w_b = w[i_b]
						if (any(w_b == 1, na.rm = TRUE) && any(w_b == 0, na.rm = TRUE)) {
							if (!private$any_censoring) break
							dead_b_temp = dead[i_b]
							if (any(dead_b_temp[w_b == 1] == 1) && any(dead_b_temp[w_b == 0] == 1) &&
								min(y[i_b]) > 0) break
						}
						attempt = attempt + 1
						if (attempt > max_resample_attempts) break
					}

					if (attempt > max_resample_attempts) {
						return(NA_real_)
					}

					boot_inf_obj = self$duplicate()
					boot_inf_obj$.__enclos_env__$private$y = y[i_b]
					boot_inf_obj$.__enclos_env__$private$dead = dead[i_b]
					boot_inf_obj$.__enclos_env__$private$w = w[i_b]
					boot_inf_obj$.__enclos_env__$private$X = X[i_b, , drop = FALSE]
					boot_inf_obj$.__enclos_env__$private$cached_values = list()

					tryCatch({
						boot_inf_obj$compute_treatment_estimate()
					}, error = function(e) {
						NA_real_
					})
				}, mc.cores = private$num_cores))
				return(beta_hat_T_bs)
			}
		},

		# @description
		# Computes a 1-alpha level frequentist bootstrap confidence interval differently for all response types, estimate types and test types.
		#
		# @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		# @param B						Number of bootstrap samples. The default is NA which corresponds to B=501.
		# @param na.rm 				Should we remove beta_hat_T's that are NA's? Default is \code{FALSE}.
		#
		# @return 	A (1 - alpha)-sized frequentist confidence interval for the treatment effect
		#
		# @examples
		# \dontrun{
		# seq_des = SeqDesignCRD$new(n = 6)
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		# seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#
		# seq_des_inf = SeqDesignInference$new(seq_des)
		# seq_des_inf$compute_asymp_confidence_interval()
		# }
		#
		compute_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, na.rm = FALSE){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertLogical(na.rm)
			beta_hat_T_bs = private$get_or_cache_bootstrap_samples(B)
			na_bs = !is.finite(beta_hat_T_bs)
			if (!na.rm && any(na_bs)) {
				warning("Bootstrap samples contain NA/NaN/Inf; dropping non-finite values for confidence interval computation.")
			}
			beta_hat_T_bs = beta_hat_T_bs[!na_bs]
			if (length(beta_hat_T_bs) == 0) {
				return(c(NA_real_, NA_real_))
			}
			quantile(beta_hat_T_bs, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
		},

		# @description
		# Computes a bootstrap two-sided p-value for H_0: betaT = delta.
		# It does so differently for all response types, estimate types and test types.
		#
		# @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		# @param B						Number of bootstrap samples. The default is NA which corresponds to B=501.
		# @param na.rm 				Should we remove beta_hat_T's that are NA's? Default is \code{FALSE}.
		#
		# @return 	The approximate frequentist p-value
		#
		# @examples
		# \dontrun{
		# seq_des = SeqDesignCRD$new(n = 6)
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		# seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#
		# seq_des_inf = SeqDesignInference$new(seq_des)
		# seq_des_inf$compute_bootstrap_two_sided_pval()
		# }
		#
		compute_bootstrap_two_sided_pval = function(delta = 0, B = 501, na.rm = FALSE){
			assertNumeric(delta)
			assertLogical(na.rm)

			beta_hat_obs = self$compute_treatment_estimate()
			if (is.na(beta_hat_obs)) return(NA_real_)
			beta_hat_T_bs = private$get_or_cache_bootstrap_samples(B)
			na_bs = !is.finite(beta_hat_T_bs)
			if (!na.rm && any(na_bs)) {
				warning("Bootstrap samples contain NA/NaN/Inf; dropping non-finite values for p-value computation.")
			}
			beta_hat_T_bs = beta_hat_T_bs[!na_bs]
			if (length(beta_hat_T_bs) == 0) return(NA_real_)

			# Percentile bootstrap two-sided p-value -- consistent with compute_bootstrap_confidence_interval,
			# which uses the percentile CI [Q_{alpha/2}, Q_{1-alpha/2}].
			# delta lies inside the (1-alpha) CI if and only if this p-value >= alpha.
			# Floor at 2/B (the minimum resolvable p-value with B bootstrap samples), avoiding
			# exact 0 which is a discretization artifact and not a meaningful probability.
			n_bs = length(beta_hat_T_bs)
			p_left  = mean(beta_hat_T_bs <= delta, na.rm = TRUE)
			p_right = mean(beta_hat_T_bs >= delta, na.rm = TRUE)
			min(1, max(2 / n_bs, 2 * min(p_left, p_right)))
		},

		# @description
		# Under the sharp null of
		# forall i H_0: y_i_T - y_i_C = delta
		# there will be a distribution of the estimates of the treatment effect (over many realizations of assignments)
		#
		# @param nsim_exact_test		The number of randomization vectors to use. The default is 501.
		# @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		# @param transform_responses	"none" for no transformation (default), "log" for log (your option when response type is survival), "logit" for logit (your option when response type is proportion) and "log1p" for log(y+1) (your option when response type is count). This is mostly an
		#								internal parameter set to something besides "none" when computing randomization confidence intervals
		#								for non-continuous response types.
		# @param show_progress		Show a text progress bar when running in serial. Ignored for parallel execution.
		# @param permutations		(Optional) A list of pre-computed randomization permutations to use for the exact test.
		# 							If NULL (default), permutations will be drawn on the fly.
		# @return 	A vector of size \code{nsim_exact_test} that has the values of beta_hat_T over many w draws.
		#
		# @examples
		# \dontrun{
		# seq_des = SeqDesignCRD$new(n = 6)
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		# seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#
		# seq_des_inf = SeqDesignInference$new(seq_des)
		# beta_hat_T_diff_ws = seq_des_inf$compute_beta_hat_T_randomization_distr_under_sharp_null()
		# ggplot(data.frame(beta_hat_T_diff_ws = beta_hat_T_diff_ws)) +
		#   geom_histogram(aes(x = beta_hat_T_diff_ws))
		# }
		compute_beta_hat_T_randomization_distr_under_sharp_null = function(nsim_exact_test = 501, delta = 0, transform_responses = "none", show_progress = TRUE, permutations = NULL){
			assertNumeric(delta)
			assertCount(nsim_exact_test, positive = TRUE)

			if (is.null(permutations)) {
				permutations = private$generate_permutations(nsim_exact_test)
			}

			bypass_checks = FALSE
			if (transform_responses == "already_transformed"){
				transform_responses = "none"
				bypass_checks = TRUE
			}

			# Create a template design object with transformations applied
			# This will be duplicated by each thread
			seq_des_template = private$seq_des_obj$duplicate()

			# Apply response transformation if necessary. Count models that operate on the
			# original response scale still use transform_responses = "log" to indicate a
			# multiplicative null, but the responses themselves must remain on the count scale.
			if (transform_responses == "log"){
				if (private$seq_des_obj_priv_int$response_type != "count") {
					seq_des_template$.__enclos_env__$private$y = log(copy(private$y))
				}
			} else if (transform_responses == "logit"){
				seq_des_template$.__enclos_env__$private$y = logit(copy(private$y))
			} else if (transform_responses == "log1p"){
				seq_des_template$.__enclos_env__$private$y = log1p(copy(private$y))
			}

			# Apply delta adjustment if nonzero
			if (delta != 0 && !bypass_checks){
				if (private$seq_des_obj_priv_int$response_type == "incidence"){
					stop("randomization tests with delta nonzero are not supported for incidence response type")
				}
				if (private$seq_des_obj_priv_int$response_type == "count" &&
					!(transform_responses %in% c("log1p", "already_transformed", "log"))){
					stop("randomization tests with delta nonzero are not supported for count response type without log1p transform")
				}
				if (private$seq_des_obj_priv_int$response_type == "proportion" && transform_responses != "logit"){
					stop("randomization tests with delta nonzero are not supported for proportion response type without logit transform (values must remain in (0,1))")
				}
				if (private$seq_des_obj_priv_int$response_type == "survival" && transform_responses != "log"){
					stop("randomization tests with delta nonzero are not supported for survival response type without log transform (values must be positive)")
				}
				if (transform_responses == "log" &&
						private$seq_des_obj_priv_int$response_type == "count") {
					# Multiplicative sharp null: H_0: y_T = y_C * exp(delta)
					# Multiply treatment responses by exp(-delta) so template holds y_C for everyone.
					# This keeps counts non-negative (unlike additive subtraction which can yield y < 0).
					seq_des_template$.__enclos_env__$private$y[private$w == 1] =
						seq_des_template$.__enclos_env__$private$y[private$w == 1] * exp(-delta)
				} else {
					# Testing H_0: y_T_i - y_C_i = delta <=> (y_T_i - delta) - y_C_i = 0
					# So adjust treatment responses downward by delta on the (possibly transformed) scale
					seq_des_template$.__enclos_env__$private$y[private$w == 1] =
						seq_des_template$.__enclos_env__$private$y[private$w == 1] - delta
					# Inverse-transform back to the original scale so downstream models receive
					# valid inputs: survreg/coxph require y > 0; betareg requires y in (0,1).
					# exp(log(y_T) - delta) = y_T * exp(-delta) > 0 always
					# inv_logit(logit(y_T) - delta) is always in (0,1)
					if (transform_responses == "log" &&
							private$seq_des_obj_priv_int$response_type != "count"){
						seq_des_template$.__enclos_env__$private$y = exp(seq_des_template$.__enclos_env__$private$y)
					} else if (transform_responses == "logit"){
						seq_des_template$.__enclos_env__$private$y = inv_logit(seq_des_template$.__enclos_env__$private$y)
					} else if (transform_responses == "log1p"){
						seq_des_template$.__enclos_env__$private$y = expm1(seq_des_template$.__enclos_env__$private$y)
					}
				}
			}

			# Store is_KK flag for use in callback
			is_KK = private$is_KK
			custom_stat_analysis = private$analyze_custom_randomization_statistic()
			use_lightweight_custom_stat = isTRUE(custom_stat_analysis$can_use_lightweight_yw_only)
			custom_stat_needs_match_data = isTRUE(custom_stat_analysis$needs_match_data)
			lightweight_custom_context = NULL
			base_template_y = seq_des_template$.__enclos_env__$private$y
			base_template_dead = seq_des_template$.__enclos_env__$private$dead
			if (use_lightweight_custom_stat) {
				lightweight_custom_context = private$build_lightweight_custom_randomization_context()
			}

			use_perms = !is.null(permutations) && length(permutations) >= nsim_exact_test

			# Iteration function that uses thread-local objects
			run_randomization_iteration <- function(thread_des_obj, thread_inf_obj, perm_idx = NULL){
				if (use_lightweight_custom_stat && use_perms && !is.null(perm_idx)) {
					perm_data = permutations[[perm_idx]]
					w_sim = perm_data$w
					y_sim = base_template_y
					if (delta != 0) {
						if (transform_responses == "log") {
							y_sim[w_sim == 1] = y_sim[w_sim == 1] * exp(delta)
						} else if (transform_responses == "logit") {
							y_sim[w_sim == 1] = inv_logit(logit(y_sim[w_sim == 1]) + delta)
						} else if (transform_responses == "log1p") {
							y_sim[w_sim == 1] = (y_sim[w_sim == 1] + 1) * exp(delta) - 1
						} else {
							y_sim[w_sim == 1] = y_sim[w_sim == 1] + delta
						}
					}

					estimate = private$evaluate_lightweight_custom_randomization_statistic(
						context = lightweight_custom_context,
						y = y_sim,
						w = w_sim,
						dead = base_template_dead
					)
					if (identical(estimate, NA_real_)) return(NA_real_)
					if (is.numeric(estimate) && length(estimate) == 1 && is.finite(estimate)) {
						return(estimate)
					}
					warning("run_randomization_iteration returned unexpected type or length.")
					return(NA_real_)
				}

				if (use_perms && !is.null(perm_idx)) {
					# Use cached w
					perm_data = permutations[[perm_idx]]

					# Update design object directly with cached assignments so custom stats
					# accessing seq_des_obj_priv_int$w get the permuted w.
					thread_des_obj$.__enclos_env__$private$w = perm_data$w
					if (is_KK && private$object_has_private_method(thread_des_obj, "match_indic")){
						thread_des_obj$.__enclos_env__$private$match_indic = perm_data$match_indic
					}

					# Update inference object directly with cached assignments
					thread_inf_obj$.__enclos_env__$private$w = perm_data$w
					if (is_KK && private$object_has_private_method(thread_inf_obj, "match_indic")){
						thread_inf_obj$.__enclos_env__$private$match_indic = perm_data$match_indic
					}

					# We still need to link the design object private environment if not already linked,
					# or ensure y/dead are correct.
					# thread_des_obj is not used for redraw here, but its 'y' and 'dead' are used below.
					# ensure thread_inf_obj has the correct y/dead from thread_des_obj (which has the delta shift).
					thread_inf_obj$.__enclos_env__$private$seq_des_obj_priv_int = thread_des_obj$.__enclos_env__$private
					thread_inf_obj$.__enclos_env__$private$y = thread_des_obj$.__enclos_env__$private$y
					thread_inf_obj$.__enclos_env__$private$dead = thread_des_obj$.__enclos_env__$private$dead

				} else {
					# Redraw treatment assignments
					thread_des_obj$.__enclos_env__$private$redraw_w_according_to_design()

					# Update inference object with new assignments
					thread_inf_obj$.__enclos_env__$private$seq_des_obj_priv_int = thread_des_obj$.__enclos_env__$private
					# Also update the direct private fields that compute_treatment_estimate uses
					thread_inf_obj$.__enclos_env__$private$y = thread_des_obj$.__enclos_env__$private$y
					thread_inf_obj$.__enclos_env__$private$w = thread_des_obj$.__enclos_env__$private$w
					thread_inf_obj$.__enclos_env__$private$dead = thread_des_obj$.__enclos_env__$private$dead
				}

				thread_inf_obj$.__enclos_env__$private$cached_values = list()

				# For delta != 0: the template holds y_C for everyone (treated responses were
				# shifted down by delta in the setup above).  Add delta back for the newly-
				# assigned treatment subjects so the null distribution is centred at delta
				# (not at 0).  This is equivalent to observing y_C + delta * w_new under the
				# sharp null H_0: y_{T,i} - y_{C,i} = delta.
				if (delta != 0) {
					w_sim = thread_inf_obj$.__enclos_env__$private$w
					y_sim = thread_inf_obj$.__enclos_env__$private$y  # y_C for all
					if (transform_responses == "log") {
						y_sim[w_sim == 1] = y_sim[w_sim == 1] * exp(delta)
					} else if (transform_responses == "logit") {
						y_sim[w_sim == 1] = inv_logit(logit(y_sim[w_sim == 1]) + delta)
					} else if (transform_responses == "log1p") {
						y_sim[w_sim == 1] = (y_sim[w_sim == 1] + 1) * exp(delta) - 1
					} else {
						y_sim[w_sim == 1] = y_sim[w_sim == 1] + delta
					}
					thread_inf_obj$.__enclos_env__$private$y = y_sim
				}

				# For KK designs, recompute match data and match statistics (the latter only if necessary)
				if (is_KK &&
						(is.null(private$custom_randomization_statistic_function) || custom_stat_needs_match_data) &&
						private$object_has_private_method(thread_inf_obj, "compute_basic_match_data")){
					# if using cache, match_indic is already set. If not, it comes from thread_des_obj
					if (!use_perms) {
						thread_inf_obj$.__enclos_env__$private$match_indic = thread_des_obj$.__enclos_env__$private$match_indic
					}
					thread_inf_obj$.__enclos_env__$private$compute_basic_match_data()
					if (private$object_has_private_method(thread_inf_obj, "compute_reservoir_and_match_statistics")){
						thread_inf_obj$.__enclos_env__$private$compute_reservoir_and_match_statistics()
					}
				}

				# Compute and return treatment estimate; catch errors from degenerate permutations
				estimate = tryCatch(
					thread_inf_obj$.__enclos_env__$private$compute_treatment_estimate_during_randomization_inference(),
					error = function(e) NA_real_
				)
				if (identical(estimate, NA_real_)) return(NA_real_)
				if (is.list(estimate) && "b" %in% names(estimate) && length(estimate$b) > 0) {
					estimate_val = as.numeric(estimate$b[1]) # Assuming first coef is the effect
				} else if (is.numeric(estimate) && length(estimate) == 1) {
					estimate_val = estimate
				} else {
					# This should ideally not happen if compute_treatment_estimate_during_randomization_inference
					# is well-behaved, but for robustness
					warning("run_randomization_iteration returned unexpected type or length.")
					return(NA_real_)
				}
				if (!is.finite(estimate_val)) {
					return(NA_real_)
				}
				estimate_val
			}

			# Check if subclass provides a C++ OpenMP dispatcher to bypass the slow R loop
			if (use_perms && private$has_private_method("compute_fast_randomization_distr")) {
				fast_distr = private$compute_fast_randomization_distr(
					seq_des_template$.__enclos_env__$private$y,
					permutations,
					delta,
					transform_responses
				)
				if (!is.null(fast_distr)) {
					return(fast_distr)
				}
			}

			# Pure R implementation for robustness
			# We use a single pair of objects if cores = 1
			# Parallelizing randomization tests in R can be done via mclapply if needed
			if (private$num_cores == 1) {
				thread_des_obj = NULL
				thread_inf_obj = NULL
				if (!use_lightweight_custom_stat) {
					thread_des_obj = seq_des_template$duplicate()
					thread_inf_obj = self$duplicate()
				}
				beta_hat_T_diff_ws = rep(NA_real_, nsim_exact_test)
				pb = NULL
				if (isTRUE(show_progress)) {
					pb = utils::txtProgressBar(min = 0, max = nsim_exact_test, style = 3)
					on.exit(close(pb), add = TRUE)
				}
				collected = 0
				total_attempts = 0
				max_attempts = nsim_exact_test * 10L

				# If using cache, we just iterate through indices
				if (use_perms) {
					for (i in 1:nsim_exact_test) {
						val = run_randomization_iteration(thread_des_obj, thread_inf_obj, perm_idx = i)
						beta_hat_T_diff_ws[i] = val
						if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
					}
				} else {
					while (collected < nsim_exact_test && total_attempts < max_attempts) {
						val = run_randomization_iteration(thread_des_obj, thread_inf_obj)
						total_attempts = total_attempts + 1L
						if (!is.na(val)) {
							collected = collected + 1L
							beta_hat_T_diff_ws[collected] = val
							if (!is.null(pb)) {
								utils::setTxtProgressBar(pb, collected)
							}
						}
					}
				}
			} else {
				# Use parallel backend if more than 1 core requested
				beta_hat_T_diff_ws = unlist(parallel::mclapply(1:nsim_exact_test, function(i) {
					thread_des_obj = NULL
					thread_inf_obj = NULL
					if (!use_lightweight_custom_stat) {
						thread_des_obj = seq_des_template$duplicate()
						thread_inf_obj = self$duplicate()
					}
					run_randomization_iteration(thread_des_obj, thread_inf_obj, perm_idx = if(use_perms) i else NULL)
				}, mc.cores = private$num_cores))
			}

			# Explicitly convert to numeric to catch non-numeric elements
			if (!is.numeric(beta_hat_T_diff_ws)) {
				if (private$verbose) {
					cat("        WARNING: beta_hat_T_diff_ws is not numeric (type:", typeof(beta_hat_T_diff_ws), "). Coercing to numeric.\n")
				}
				beta_hat_T_diff_ws = as.numeric(beta_hat_T_diff_ws)
			}

			beta_hat_T_diff_ws
		},

		# @description
		# Fisher's randomization test which means that H_0: y_i_T - y_i_C = delta for all subjects
		# either the classic different-in-means estimate of the additive treatment effect,
		# i.e. ybar_T - ybar_C or the default_regression estimate of the additive treatment effect linearly i.e.
		# the treatment different adjusted linearly for the p covariates.
		#
		# @param nsim_exact_test		The number of randomization vectors to use in the randomization test (ignored if \code{test_type}
		# 								is not "randomization-exact"). The default is 501 providing pvalue resolution to a fifth of a percent.
		# @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		# @param transform_responses	"none" for no transformation (default), "log" for log (your option when response type is survival), "logit" for logit (your option when response type is proportion) and "log1p" for log(y+1) (your option when response type is count). This is mostly an
		#								internal parameter set something besides "none" when computing randomization confidence intervals
		#								for non-continuous responses.
		# @param na.rm 				Should we remove beta_hat_T's that are NA's? Default is \code{TRUE}.
		# @param show_progress		Show a text progress bar when running in serial. Ignored for parallel execution.
		# @param permutations		(Optional) A list of pre-computed randomization permutations to use for the exact test.
		# 							If NULL (default), permutations will be drawn on the fly.
		#
		# @return 	The approximate frequentist p-value
		#
		# @examples
		# \dontrun{
		# seq_des = SeqDesignCRD$new(n = 6)
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		# seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#
		# seq_des_inf = SeqDesignInference$new(seq_des)
		# seq_des_inf$compute_asymp_two_sided_pval_for_treatment_effect()
		# }
		#
		compute_two_sided_pval_for_treatment_effect_rand = function(nsim_exact_test = 501, delta = 0, transform_responses = "none", na.rm = TRUE, show_progress = TRUE, permutations = NULL){
			assertLogical(na.rm)
			if (private$seq_des_obj_priv_int$response_type == "incidence"){
				stop("Randomization tests are not supported for incidence outcomes. Use SeqDesignInferenceIncidExactZhang for the Zhang method.")
			}
			if (is.null(permutations)) {
				permutations = private$generate_permutations(nsim_exact_test)
			}

			cache_key = private$build_randomization_distribution_cache_key(
				nsim_exact_test = nsim_exact_test,
				delta = delta,
				transform_responses = transform_responses,
				permutations = permutations
			)

			# Fast path for linear estimators (continuous, no custom statistic):
			# t0s(delta) = t0s(0) + delta exactly (mean diff, OLS, KK compound).
			# Reuses cached t0s from a prior delta=0 call, so no new simulations are needed.
			if (transform_responses == "none" &&
					is.null(private[["custom_randomization_statistic_function"]]) &&
					!is.null(private$cached_values$t0s_rand) &&
					length(private$cached_values$t0s_rand) >= nsim_exact_test) {
				t0s = private$cached_values$t0s_rand[seq_len(nsim_exact_test)] + delta
				t = private$compute_treatment_estimate_during_randomization_inference()
				if (length(t) != 1 || !is.finite(t)) return(NA_real_)
				na_t0s = !is.finite(t0s)
				if (length(t0s) > 0 && !na.rm && any(na_t0s)) {
					warning("Randomization distribution contains NA; dropping NA values for p-value computation.")
				}
				nsim_adj = sum(!na_t0s)
				if (nsim_adj == 0L) return(NA_real_)
				return(min(1, max(2 / nsim_adj, 2 * min(
					sum(t0s >= t, na.rm = TRUE) / nsim_adj,
					sum(t0s <= t, na.rm = TRUE) / nsim_adj
				))))
			}

			if (is.null(private$cached_values$rand_distr_cache)) {
				private$cached_values$rand_distr_cache = list()
			}

			if (!is.null(cache_key) &&
					!is.null(private$cached_values$rand_distr_cache[[cache_key]]) &&
					length(private$cached_values$rand_distr_cache[[cache_key]]) >= nsim_exact_test) {
				t0s = private$cached_values$rand_distr_cache[[cache_key]][seq_len(nsim_exact_test)]
			} else {
				#approximate the null distribution by computing estimates on many draws of w
				t0s = self$compute_beta_hat_T_randomization_distr_under_sharp_null(
					nsim_exact_test,
					delta,
					transform_responses,
					show_progress = show_progress,
					permutations = permutations
				)
				if (!is.null(cache_key)) {
					private$cached_values$rand_distr_cache[[cache_key]] = t0s
				}
			}

			# Cache t0s when delta=0 on the raw scale for future exact-shift fast-path use.
			if (delta == 0 && transform_responses == "none" && is.null(private[["custom_randomization_statistic_function"]])) {
				private$cached_values$t0s_rand = t0s
			}

			#this calculates the actual estimate to compare against the null distribution
			t = private$compute_treatment_estimate_during_randomization_inference()
			#finally compute the p-value
			if (private$verbose) {
			  cat("        Randomization Pval - t0s (length", length(t0s), ", NAs:", sum(is.na(t0s)), "): ", head(t0s), "\n")
			  cat("        Randomization Pval - t:", t, "\n")
			}

			# Handle case where t might be a vector or non-finite
			if (length(t) != 1 || !is.finite(t)) return(NA_real_)

			na_t0s = !is.finite(t0s)
			if (length(t0s) > 0 && !na.rm && any(na_t0s)) {
				warning("Randomization distribution contains NA; dropping NA values for p-value computation.")
			}
			nsim_exact_test_adjusted = sum(!na_t0s)

			if (nsim_exact_test_adjusted == 0) return(NA_real_) # Avoid division by zero

			# Two-sided p-value capped at 1, floored at 2/nsim to avoid exact 0 (discretization artifact).
			min(1, max(2 / nsim_exact_test_adjusted, 2 * min(
				sum(t0s >= t, na.rm = TRUE) / nsim_exact_test_adjusted,
				sum(t0s <= t, na.rm = TRUE) / nsim_exact_test_adjusted
			)))
		},


		# @description
		# Fisher's randomization test which means that H_0: y_i_T - y_i_C = delta for all subjects

		# @description
		# Computes a 1-alpha level frequentist confidence interval for the randomization test
		#
		# Here we invert the randomization test that tests the strong null H_0: y_T_i - y_C_i = delta <=> (y_T_i - delta) - y_C_i = 0 so
		# we adjust the treatment responses downward by delta. We then find the set of all delta values that is above 1 - alpha/2 (i.e. two-sided)
		# This is accomplished via a bisection algorithm (algorithm 1 of Glazer and Stark, 2025 available at
		# https://arxiv.org/abs/2405.05238). These confidence intervals are exact to within tolerance \code{pval_epsilon}.
		# As far as we know, this works for response types continuous, ordinal (cumulative-link or adjacent-category logit/probit estimates), uncensored survival (where we work in log-time and then
		# return to natural time when finished), count (where we work in log1p-count and then return to natural count when finished) and proportion (where we work in logit-rate and return to rate when finished).
		#
		# @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		# @param nsim_exact_test		The number of randomization vectors (applicable for test type "randomization-exact" only).
		# 								The default is 1000 providing good resolutions to confidence intervals.
		# @param pval_epsilon			The bisection algorithm tolerance for the test inversion (applicable for test type "randomization-exact" only).
		# 								The default is to find a CI accurate to within 0.005.
		#
		# @param show_progress		Show a text progress indicator for the bisection p-value span. Ignored for parallel execution.
		# @return 	A 1 - alpha sized frequentist confidence interval for the treatment effect
		#
		# @examples
		# \dontrun{
		# seq_des = SeqDesignCRD$new(n = 6)
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		# seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		# seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#
		# seq_des_inf = SeqDesignInference$new(seq_des)
		# seq_des_inf$compute_asymp_confidence_interval()
		# }
		#
		compute_confidence_interval_rand = function(alpha = 0.05, nsim_exact_test = 501, pval_epsilon = 0.005, show_progress = TRUE){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertCount(nsim_exact_test, positive = TRUE)
			assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
			assertLogical(show_progress)
			show_progress = isTRUE(show_progress) && private$num_cores == 1

			# TRUE for any model that operates on the original response scale (GLMs, GEE, GLMM).
			# These branches keep y untransformed and pass a transform_arg to the randomization test.
			# FALSE for simple mean-difference estimators that pre-transform y (e.g. logit, log).
			is_model_on_original_scale = inherits(self, "SeqDesignInferenceMLEorKMforGLMs") ||
				inherits(self, "SeqDesignInferenceAbstractKKGEE") ||
				inherits(self, "SeqDesignInferenceAbstractKKGLMM")

			build_randomization_ci_search_bounds = function(inf_obj,
					transform_arg,
					permutations,
					target_pval){
				obj_private = inf_obj$.__enclos_env__$private

				if (transform_arg == "none" &&
						(is.null(obj_private$cached_values$t0s_rand) ||
							length(obj_private$cached_values$t0s_rand) < nsim_exact_test)) {
					inf_obj$compute_two_sided_pval_for_treatment_effect_rand(
						nsim_exact_test = nsim_exact_test,
						delta = 0,
						transform_responses = transform_arg,
						show_progress = FALSE,
						permutations = permutations
					)
				}

				est = tryCatch(
					as.numeric(inf_obj$compute_treatment_estimate()),
					error = function(e) NA_real_
				)
				if (length(est) == 0L || !is.finite(est[1])) {
					est = NA_real_
				} else {
					est = est[1]
				}

				asym_ci = tryCatch(
					as.numeric(inf_obj$compute_asymp_confidence_interval(alpha = alpha * 2)),
					error = function(e) c(NA_real_, NA_real_)
				)
				if (length(asym_ci) < 2L || !all(is.finite(asym_ci[1:2]))) {
					asym_ci = c(NA_real_, NA_real_)
				} else {
					asym_ci = sort(asym_ci[1:2])
				}

				if (!is.finite(est) && all(is.finite(asym_ci))) {
					est = mean(asym_ci)
				}

				if (!all(is.finite(asym_ci))) {
					y_vals = obj_private$y
					scale_guess = suppressWarnings(stats::sd(y_vals, na.rm = TRUE))
					if (!is.finite(scale_guess) || scale_guess <= 0) {
						scale_guess = suppressWarnings(stats::IQR(y_vals, na.rm = TRUE) / 1.349)
					}
					if (!is.finite(scale_guess) || scale_guess <= 0) {
						scale_guess = 1
					}
					if (!is.finite(est)) {
						est = suppressWarnings(stats::median(y_vals, na.rm = TRUE))
						if (!is.finite(est)) {
							est = 0
						}
					}
					asym_ci = c(est - 2 * scale_guess, est + 2 * scale_guess)
				}

				l = asym_ci[1]
				u = asym_ci[2]
				if (l >= est) {
					l = est - max(abs(u - est), 1)
				}
				if (u <= est) {
					u = est + max(abs(est - l), 1)
				}

				evaluate_pval = function(delta){
					pval = tryCatch(
						as.numeric(inf_obj$compute_two_sided_pval_for_treatment_effect_rand(
							nsim_exact_test = nsim_exact_test,
							delta = delta,
							transform_responses = transform_arg,
							show_progress = FALSE,
							permutations = permutations
						)),
						error = function(e) NA_real_
					)
					if (length(pval) == 0L) {
						return(NA_real_)
					}
					pval[1]
				}

				expand_bound = function(bound, lower){
					pval_bound = evaluate_pval(bound)
					if (is.finite(pval_bound) && pval_bound < target_pval) {
						return(bound)
					}

					step = abs(est - bound)
					if (!is.finite(step) || step <= 0) {
						step = 1
					}

					for (iter in seq_len(12L)) {
						step = step * 2
						candidate = if (lower) est - step else est + step
						pval_candidate = evaluate_pval(candidate)
						bound = candidate
						if (is.finite(pval_candidate) && pval_candidate < target_pval) {
							break
						}
					}
					bound
				}

				list(
					est = est,
					l = expand_bound(l, lower = TRUE),
					u = expand_bound(u, lower = FALSE)
				)
			}

			switch(private$seq_des_obj_priv_int$response_type,
				continuous = {
					perms = private$generate_permutations(nsim_exact_test)
					bounds = build_randomization_ci_search_bounds(
						inf_obj = self,
						transform_arg = "none",
						permutations = perms,
						target_pval = alpha / 2
					)

					if (private$num_cores > 1) {
						ci = private$compute_ci_both_bounds_parallel(
							nsim_exact_test = nsim_exact_test,
							l_lower = bounds$l,
							u_lower = bounds$est,
							l_upper = bounds$est,
							u_upper = bounds$u,
							pval_th = alpha / 2,
							tol = pval_epsilon,
							transform_responses = "none",
							permutations = perms
						)
					} else {
						ci = c(
							private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test,
								l = bounds$l,
								u = bounds$est,
								pval_th = alpha / 2,
								tol = pval_epsilon,
								transform_responses = "none",
								lower = TRUE,
								show_progress = show_progress,
								permutations = perms),
							private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test,
								l = bounds$est,
								u = bounds$u,
								pval_th = alpha / 2,
								tol = pval_epsilon,
								transform_responses = "none",
								lower = FALSE,
								show_progress = show_progress,
								permutations = perms)
						)
					}
					ci
				},
				incidence =  stop("Confidence intervals are not supported for randomization tests for incidence outcomes. Use SeqDesignInferenceIncidExactZhang for the Zhang method."),
				ordinal =    {
					perms = private$generate_permutations(nsim_exact_test)
					bounds = build_randomization_ci_search_bounds(
						inf_obj = self,
						transform_arg = "none",
						permutations = perms,
						target_pval = alpha / 2
					)

					if (private$num_cores > 1) {
						ci = private$compute_ci_both_bounds_parallel(
							nsim_exact_test = nsim_exact_test,
							l_lower = bounds$l,
							u_lower = bounds$est,
							l_upper = bounds$est,
							u_upper = bounds$u,
							pval_th = alpha / 2,
							tol = pval_epsilon,
							transform_responses = "none",
							permutations = perms
						)
					} else {
						ci = c(
							private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test,
								l = bounds$l,
								u = bounds$est,
								pval_th = alpha / 2,
								tol = pval_epsilon,
								transform_responses = "none",
								lower = TRUE,
								show_progress = show_progress,
								permutations = perms),
							private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test,
								l = bounds$est,
								u = bounds$u,
								pval_th = alpha / 2,
								tol = pval_epsilon,
								transform_responses = "none",
								lower = FALSE,
								show_progress = show_progress,
								permutations = perms)
						)
					}
					names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, sep = "%")
					ci
				},
				count =      {
					# Create a temporary object with transformed responses
					temp_inf = self$duplicate()
					# Prevent fork-unsafe nested mclapply: bootstrap and bisection run serially on
					# temp_inf (num_cores=1). Both bounds are still computed in parallel via an outer
					# mclapply when private$num_cores > 1 (see parallel branch below).
					temp_inf$.__enclos_env__$private$num_cores = 1L
					# Cache permutations for speed
					perms = temp_inf$.__enclos_env__$private$generate_permutations(nsim_exact_test)

					transform_arg = "already_transformed"

					if (is_model_on_original_scale){
						transform_arg = "log"
						# y remains counts
					} else {
						# No need to clamp 0 for log1p, it handles 0 naturally
						temp_inf$.__enclos_env__$private$y = log1p(temp_inf$.__enclos_env__$private$y)
					}

					temp_inf$.__enclos_env__$private$cached_values = list()
					bounds = build_randomization_ci_search_bounds(
						inf_obj = temp_inf,
						transform_arg = transform_arg,
						permutations = perms,
						target_pval = alpha / 2
					)

					# Run inversion on temp_inf which has transformed data.
					# Each bound's bisection is serial (temp_inf$num_cores=1); both bounds may be
					# computed simultaneously via an outer fork when private$num_cores > 1.
					if (private$num_cores > 1) {
						bound_specs = list(
							list(l = bounds$l, u = bounds$est, lower = TRUE),
							list(l = bounds$est, u = bounds$u, lower = FALSE)
						)
						results = parallel::mclapply(bound_specs, function(spec) {
							temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(
								nsim_exact_test, l = spec$l, u = spec$u,
								pval_th = alpha / 2, tol = pval_epsilon,
								transform_responses = transform_arg, lower = spec$lower,
								show_progress = FALSE, permutations = perms)
						}, mc.cores = min(2L, private$num_cores))
						ci = c(results[[1]], results[[2]])
					} else {
						ci = c(
							temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test,
								l = bounds$l, u = bounds$est, pval_th = alpha / 2, tol = pval_epsilon,
								transform_responses = transform_arg, lower = TRUE,
								show_progress = show_progress, permutations = perms),
							temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test,
								l = bounds$est, u = bounds$u, pval_th = alpha / 2, tol = pval_epsilon,
								transform_responses = transform_arg, lower = FALSE,
								show_progress = show_progress, permutations = perms)
						)
					}
					# ci = exp(ci) # Ratio of (1+y)
				},
				proportion = {
					# Create a temporary object with transformed responses
					temp_inf = self$duplicate()
					# Prevent fork-unsafe nested mclapply: bootstrap and bisection run serially on
					# temp_inf (num_cores=1). Both bounds are still computed in parallel via an outer
					# mclapply when private$num_cores > 1 (see parallel branch below).
					temp_inf$.__enclos_env__$private$num_cores = 1L
					# Cache permutations for speed
					perms = temp_inf$.__enclos_env__$private$generate_permutations(nsim_exact_test)

					transform_arg = "already_transformed"

					# Clamp to (epsilon, 1-epsilon) to avoid 0/1 for both GLM (Beta Reg) and Logit
					y_clamped = pmax(.Machine$double.eps, pmin(1 - .Machine$double.eps, temp_inf$.__enclos_env__$private$y))

					if (is_model_on_original_scale){
						transform_arg = "logit"
						# y remains proportions but clamped
						temp_inf$.__enclos_env__$private$y = y_clamped
					} else {
						temp_inf$.__enclos_env__$private$y = logit(y_clamped)
					}

					temp_inf$.__enclos_env__$private$cached_values = list()
					bounds = build_randomization_ci_search_bounds(
						inf_obj = temp_inf,
						transform_arg = transform_arg,
						permutations = perms,
						target_pval = alpha / 2
					)

					# Run inversion on temp_inf which has transformed data.
					# Each bound's bisection is serial (temp_inf$num_cores=1); both bounds may be
					# computed simultaneously via an outer fork when private$num_cores > 1.
					if (private$num_cores > 1) {
						bound_specs = list(
							list(l = bounds$l, u = bounds$est, lower = TRUE),
							list(l = bounds$est, u = bounds$u, lower = FALSE)
						)
						results = parallel::mclapply(bound_specs, function(spec) {
							temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(
								nsim_exact_test, l = spec$l, u = spec$u,
								pval_th = alpha / 2, tol = pval_epsilon,
								transform_responses = transform_arg, lower = spec$lower,
								show_progress = FALSE, permutations = perms)
						}, mc.cores = min(2L, private$num_cores))
						ci = c(results[[1]], results[[2]])
					} else {
						ci = c(
							temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test,
								l = bounds$l, u = bounds$est, pval_th = alpha / 2, tol = pval_epsilon,
								transform_responses = transform_arg, lower = TRUE,
								show_progress = show_progress, permutations = perms),
							temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test,
								l = bounds$est, u = bounds$u, pval_th = alpha / 2, tol = pval_epsilon,
								transform_responses = transform_arg, lower = FALSE,
								show_progress = show_progress, permutations = perms)
						)
					}
					# ci = exp(ci) # Odds Ratio - Removed for consistency with logit estimate
				},
				survival =   {
					# Create a temporary object with transformed responses
					temp_inf = self$duplicate()
					# Prevent fork-unsafe nested mclapply: bootstrap and bisection run serially on
					# temp_inf (num_cores=1). Both bounds are still computed in parallel via an outer
					# mclapply when private$num_cores > 1 (see parallel branch below).
					temp_inf$.__enclos_env__$private$num_cores = 1L
					# Cache permutations for speed
					perms = temp_inf$.__enclos_env__$private$generate_permutations(nsim_exact_test)

					transform_arg = "already_transformed"

					if (is_model_on_original_scale){
						transform_arg = "log"
						# y remains time
					} else {
						# Clamp 0 to epsilon (though survival times should be >0)
						y_clamped = pmax(.Machine$double.eps, temp_inf$.__enclos_env__$private$y)
						temp_inf$.__enclos_env__$private$y = log(y_clamped)
					}

					temp_inf$.__enclos_env__$private$cached_values = list()
					bounds = build_randomization_ci_search_bounds(
						inf_obj = temp_inf,
						transform_arg = transform_arg,
						permutations = perms,
						target_pval = alpha / 2
					)

					# Run inversion on temp_inf which has transformed data.
					# Each bound's bisection is serial (temp_inf$num_cores=1); both bounds may be
					# computed simultaneously via an outer fork when private$num_cores > 1.
					if (private$num_cores > 1) {
						bound_specs = list(
							list(l = bounds$l, u = bounds$est, lower = TRUE),
							list(l = bounds$est, u = bounds$u, lower = FALSE)
						)
						results = parallel::mclapply(bound_specs, function(spec) {
							temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(
								nsim_exact_test, l = spec$l, u = spec$u,
								pval_th = alpha / 2, tol = pval_epsilon,
								transform_responses = transform_arg, lower = spec$lower,
								show_progress = FALSE, permutations = perms)
						}, mc.cores = min(2L, private$num_cores))
						ci = c(results[[1]], results[[2]])
					} else {
						ci = c(
							temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test,
								l = bounds$l, u = bounds$est, pval_th = alpha / 2, tol = pval_epsilon,
								transform_responses = transform_arg, lower = TRUE,
								show_progress = show_progress, permutations = perms),
							temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test,
								l = bounds$est, u = bounds$u, pval_th = alpha / 2, tol = pval_epsilon,
								transform_responses = transform_arg, lower = FALSE,
								show_progress = show_progress, permutations = perms)
						)
					}
					# ci = exp(ci) # Ratio - Removed for consistency with log estimate
				})
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, sep = "%")
			ci
		},

		# @description
		# Duplicate this inference object
		#
		# @param verbose 	A flag indicating whether messages should be displayed to the user. Default is \code{FALSE}
		# @return 			A new `SeqDesignInference` object with the same data
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
			i$.__enclos_env__$private$cached_values$permutations_cache = private$cached_values$permutations_cache
			if (!is.null(i$.__enclos_env__$private$custom_randomization_statistic_function)){
				environment(i$.__enclos_env__$private$custom_randomization_statistic_function) = environment(i$initialize)
			}
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

	build_randomization_distribution_cache_key = function(nsim_exact_test, delta, transform_responses, permutations){
		if (is.null(permutations)) {
			return(NULL)
		}
		stat_sig = private$custom_randomization_statistic_signature()
		perm_sig = private$permutations_signature(permutations, nsim_exact_test)
		paste(
			"nsim", as.integer(nsim_exact_test),
			"delta", formatC(delta, digits = 17, format = "fg", flag = "#"),
			"transform", transform_responses,
			"stat", stat_sig,
			"perm", perm_sig,
			sep = "|"
		)
	},

	custom_randomization_statistic_signature = function(){
		if (is.null(private$custom_randomization_statistic_function)) {
			return("default")
		}
		private$stable_signature(list(
			formals = formals(private$custom_randomization_statistic_function),
			body = body(private$custom_randomization_statistic_function)
		))
	},

	permutations_signature = function(permutations, nsim_exact_test){
		private$stable_signature(list(
			nsim_exact_test = as.integer(nsim_exact_test),
			permutations = permutations[seq_len(min(length(permutations), nsim_exact_test))]
		))
	},

	stable_signature = function(obj){
		raw_sig = serialize(obj, NULL, xdr = FALSE)
		ints = as.integer(raw_sig)
		if (length(ints) == 0L) {
			return("0:0:0")
		}

		modulus = 2147483647
		h1 = 0
		h2 = 0
		step = max(1L, floor(length(ints) / 64L))
		for (i in seq_along(ints)) {
			val = ints[i]
			h1 = (h1 * 131 + val) %% modulus
			if (i == 1L || i == length(ints) || (i %% step) == 0L) {
				h2 = (h2 * 65599 + val + i) %% modulus
			}
		}
		paste(length(ints), as.integer(h1), as.integer(h2), sep = ":")
	},

	extract_dollar_paths = function(expr){
		paths = list()
		if (is.call(expr)) {
			if (identical(expr[[1]], as.name("$")) && length(expr) == 3L) {
				path = private$resolve_dollar_path(expr)
				if (!is.null(path)) {
					paths = c(paths, list(path))
				}
			}
			for (i in seq_along(expr)[-1]) {
				paths = c(paths, private$extract_dollar_paths(expr[[i]]))
			}
		}
		paths
	},

	resolve_dollar_path = function(expr){
		if (is.symbol(expr)) {
			return(as.character(expr))
		}
		if (is.call(expr) && identical(expr[[1]], as.name("$")) && length(expr) == 3L) {
			parent_path = private$resolve_dollar_path(expr[[2]])
			child_name = if (is.symbol(expr[[3]])) as.character(expr[[3]]) else NULL
			if (is.null(parent_path) || is.null(child_name)) {
				return(NULL)
			}
			return(c(parent_path, child_name))
		}
		NULL
	},

	analyze_custom_randomization_statistic = function(){
		if (!is.null(private$cached_values$custom_stat_analysis)) {
			return(private$cached_values$custom_stat_analysis)
		}
		if (is.null(private$custom_randomization_statistic_function)) {
			analysis = list(
				can_use_lightweight_yw_only = FALSE,
				needs_match_data = TRUE
			)
			private$cached_values$custom_stat_analysis = analysis
			return(analysis)
		}

		dollar_paths = private$extract_dollar_paths(body(private$custom_randomization_statistic_function))
		path_strings = vapply(dollar_paths, paste, character(1), collapse = "$")

		allowed_lightweight_paths = c(
			"private$y",
			"private$w",
			"private$dead",
			"private$seq_des_obj_priv_int",
			"private$seq_des_obj_priv_int$y",
			"private$seq_des_obj_priv_int$w",
			"private$seq_des_obj_priv_int$dead"
		)
		references_self = any(vapply(dollar_paths, function(path) length(path) > 0L && identical(path[1], "self"), logical(1)))
		can_use_lightweight_yw_only =
			!references_self &&
			all(path_strings %in% allowed_lightweight_paths)

		match_tokens = c(
			"match", "reservoir", "pair", "stratum", "strata",
			"matched", "discordant", "concordant"
		)
		needs_match_data = TRUE
		if (can_use_lightweight_yw_only) {
			needs_match_data = FALSE
		} else if (length(path_strings) > 0L) {
			needs_match_data = any(vapply(match_tokens, function(token) {
				any(grepl(token, path_strings, fixed = TRUE))
			}, logical(1)))
		}

		analysis = list(
			can_use_lightweight_yw_only = can_use_lightweight_yw_only,
			needs_match_data = needs_match_data
		)
		private$cached_values$custom_stat_analysis = analysis
		analysis
	},

	build_lightweight_custom_randomization_context = function(){
		if (is.null(private$custom_randomization_statistic_function)) {
			return(NULL)
		}
		orig_env = environment(private$custom_randomization_statistic_function)
		eval_env = new.env(parent = orig_env)
		private_proxy = new.env(parent = emptyenv())
		seq_priv_proxy = new.env(parent = emptyenv())
		private_proxy$seq_des_obj_priv_int = seq_priv_proxy
		eval_env$private = private_proxy
		custom_fun = private$custom_randomization_statistic_function
		environment(custom_fun) = eval_env
		list(
			fun = custom_fun,
			private_proxy = private_proxy,
			seq_des_obj_priv_int_proxy = seq_priv_proxy
		)
	},

	evaluate_lightweight_custom_randomization_statistic = function(context, y, w, dead = NULL){
		if (is.null(context)) {
			return(NA_real_)
		}
		context$private_proxy$y = y
		context$private_proxy$w = w
		context$private_proxy$dead = dead
		context$seq_des_obj_priv_int_proxy$y = y
		context$seq_des_obj_priv_int_proxy$w = w
		context$seq_des_obj_priv_int_proxy$dead = dead
		tryCatch(
			context$fun(),
			error = function(e) NA_real_
		)
	},

	has_private_method = function(method_name){
		method_name %in% names(private)
	},

			object_has_private_method = function(obj, method_name){
			method_name %in% names(obj$.__enclos_env__$private)
		},

					generate_permutations = function(nsim_exact_test = 501){
						cache_key = as.character(as.integer(nsim_exact_test))
						if (is.null(private$cached_values$permutations_cache)) {
							private$cached_values$permutations_cache = list()
						}
						if (!is.null(private$cached_values$permutations_cache[[cache_key]]) &&
								length(private$cached_values$permutations_cache[[cache_key]]) >= nsim_exact_test) {
							return(private$cached_values$permutations_cache[[cache_key]][seq_len(nsim_exact_test)])
						}

						perms = vector("list", nsim_exact_test)
						# Use a template design to generate
						des_temp = private$seq_des_obj$duplicate()

						for (i in 1:nsim_exact_test){
							des_temp$.__enclos_env__$private$redraw_w_according_to_design()
							perm_data = list(
								w = des_temp$.__enclos_env__$private$w
							)
							if (private$is_KK){
								perm_data$match_indic = des_temp$.__enclos_env__$private$match_indic
							}
							perms[[i]] = perm_data
						}
						private$cached_values$permutations_cache[[cache_key]] = perms
						perms
					},
			compute_treatment_estimate_during_randomization_inference = function(){			if (is.null(private$custom_randomization_statistic_function)){ #i.e., the default
				self$compute_treatment_estimate()
			} else {
				private$custom_randomization_statistic_function()
			}
		},

		create_design_matrix = function(){
			cbind(1, private$w, private$get_X())
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
			if (is.null(private$cached_values$boot_samples)){
				beta_samples = self$approximate_bootstrap_distribution_beta_hat_T(B)

				if (!is.numeric(beta_samples)){
					if (private$verbose) {
						cat("        ERROR: approximate_bootstrap_distribution_beta_hat_T returned non-numeric. Type:", typeof(beta_samples), "\n")
					}
					return(rep(NA_real_, B)) # Return NAs if not numeric
				}

				private$cached_values$boot_samples = beta_samples
				if (private$verbose) {
					cat("        Bootstrap samples summary (B=", B, "): \n")
					print(summary(private$cached_values$boot_samples))
				}
			} else {
				B_0 = length(private$cached_values$boot_samples)
				if (B_0 > B){
					return (private$cached_values$boot_samples[1 : B]) #send back what we need but don't reduce the cache
				} else if (B_0 < B){
					new_samples = self$approximate_bootstrap_distribution_beta_hat_T(B - B_0)

					if (!is.numeric(new_samples)){
						if (private$verbose) {
							cat("        ERROR: approximate_bootstrap_distribution_beta_hat_T (for new samples) returned non-numeric. Type:", typeof(new_samples), "\n")
						}
						# Pad with NAs if new samples are not numeric
						new_samples_padded = rep(NA_real_, B - B_0)
						private$cached_values$boot_samples = c(private$cached_values$boot_samples, new_samples_padded)
					} else {
						private$cached_values$boot_samples = c(private$cached_values$boot_samples, new_samples)
					}

					if (private$verbose) {
					    cat("        Bootstrap samples summary (B=", B, ", added", length(new_samples), "): \n")
					    print(summary(private$cached_values$boot_samples))
				    }
				}
			}
			private$cached_values$boot_samples
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

		# Parallel computation of both CI bounds using OpenMP
		compute_ci_both_bounds_parallel = function(nsim_exact_test, l_lower, u_lower, l_upper, u_upper, pval_th, tol, transform_responses, permutations = NULL){
			# Compute both CI bounds in parallel using mclapply (fork-based).
			# Calling R functions from OpenMP threads is unsafe; mclapply gives each
			# child its own R interpreter (copy-on-write fork), so it is safe even when
			# private$cached_values$t0s_rand is read by both children simultaneously.
			bound_specs = list(
				list(l = l_lower, u = u_lower, lower = TRUE),
				list(l = l_upper, u = u_upper, lower = FALSE)
			)
			results = parallel::mclapply(bound_specs, function(spec) {
				private$num_cores = 1L  # each forked child runs its bound serially
				private$compute_ci_by_inverting_the_randomization_test_iteratively(
					nsim_exact_test,
					l                  = spec$l,
					u                  = spec$u,
					pval_th            = pval_th,
					tol                = tol,
					transform_responses = transform_responses,
					lower              = spec$lower,
					show_progress      = FALSE,
					permutations       = permutations
				)
			}, mc.cores = min(2L, private$num_cores))

			c(results[[1]], results[[2]])
		},

		compute_ci_by_inverting_the_randomization_test_iteratively = function(nsim_exact_test, l, u, pval_th, tol, transform_responses, lower, show_progress = TRUE, permutations = NULL){
			# Pure R bisection loop for robustness with R6 objects and callbacks
			pval_l = as.numeric(self$compute_two_sided_pval_for_treatment_effect_rand(
				nsim_exact_test,
				delta = l,
				transform_responses = transform_responses,
				show_progress = FALSE,
				permutations = permutations
			))
			pval_u = as.numeric(self$compute_two_sided_pval_for_treatment_effect_rand(
				nsim_exact_test,
				delta = u,
				transform_responses = transform_responses,
				show_progress = FALSE,
				permutations = permutations
			))

			# If p-values are empty (numeric(0)), treat as NA
			if (length(pval_l) == 0) pval_l = NA_real_
			if (length(pval_u) == 0) pval_u = NA_real_

			# If an extreme bound causes model failure (NA p-value), progressively tighten
			# it toward the other bound (up to 30 halvings) until convergence is achieved.
			for (k in seq_len(30L)) {
				if (!is.na(pval_l) && !is.na(pval_u)) break
				if (is.na(pval_l)) {
					l = (l + u) / 2
					pval_l = as.numeric(self$compute_two_sided_pval_for_treatment_effect_rand(
						nsim_exact_test, delta = l, transform_responses = transform_responses,
						show_progress = FALSE, permutations = permutations))
				}
				if (is.na(pval_u)) {
					u = (l + u) / 2
					pval_u = as.numeric(self$compute_two_sided_pval_for_treatment_effect_rand(
						nsim_exact_test, delta = u, transform_responses = transform_responses,
						show_progress = FALSE, permutations = permutations))
				}
			}
			if (is.na(pval_l) || is.na(pval_u)) return(NA_real_)
			if (!all(is.finite(c(l, u)))) return(NA_real_)

			iter = 0
			progress_label = if (lower) "CI lower" else "CI upper"

			# Bisection loop
			repeat {
				# Check convergence
				pval_span = abs(pval_u - pval_l)
				if ((abs(u - l)) <= tol || pval_span <= tol) {
					if (isTRUE(show_progress)) {
						cat(sprintf("\r%s iter=%d pval_span=%.6g (target<=%.6g) done\n", progress_label, iter, pval_span, tol))
						utils::flush.console()
					}
					return(ifelse(lower, l, u))
				}

				# Compute midpoint
				m = (l + u) / 2.0
				pval_m = as.numeric(self$compute_two_sided_pval_for_treatment_effect_rand(
					nsim_exact_test,
					delta = m,
					transform_responses = transform_responses,
					show_progress = FALSE,
					permutations = permutations
				))

				# NA p-value at midpoint: model fails at this delta, treat as outside the CI
				# (extreme side) and move the extreme bound to m.
				if (length(pval_m) == 0L || is.na(pval_m)) {
					if (lower) { l = m; pval_l = 0 } else { u = m; pval_u = 0 }
					iter = iter + 1
					next
				}

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

				iter = iter + 1
				if (isTRUE(show_progress)) {
					pval_span = abs(pval_u - pval_l)
					cat(sprintf("\r%s iter=%d pval_span=%.6g (target<=%.6g)", progress_label, iter, pval_span, tol))
					utils::flush.console()
				}
			}
			}

		)
	)
