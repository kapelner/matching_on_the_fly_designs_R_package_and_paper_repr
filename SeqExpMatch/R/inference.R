#' Inference for A Sequential Design
#' 
#' @description
#' An R6 Class that estimates, tests and provides intervals for a treatment effect in a sequential design.
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
#' 
#' @export
SeqDesignInference = R6::R6Class("SeqDesignInference",
	public = list(
		#' @field estimate_type		The type of estimate to compute (either "mean_difference-or-medians" or "default_regression").
		estimate_type = NULL,
		#' @field test_type			The type of test to run (either "MLE-or-KM-based" or "randomization-exact").
		test_type = NULL,
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' 
		#' 
		#' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param estimate_type		The type of estimate to compute (either "mean_difference", "median_difference" or "default_regression"). 
		#' 							The "mean_difference" option for response type continuous, count and proportion assumes the target parameter is the
		#' 							difference in outcome expectation between the treatment and control arms. For incidence, this is the difference in
		#' 							log odds of probability of positive class between the treatment and control arms assuming a logistic regression model additive 
		#' 							in treatment effect. For survival, this is the difference in log outcome between the treatment and control arms assuming 
		#' 							a weibull model additive in treatment effect.
		#' 							log odds of probability of positive class between the treatment and control.
		#' 							The default is "default_regression" as this provides higher power if you feel comfortable assuming the appropriate glm.
		#' @param test_type			The type of test to run (either "MLE-or-KM-based" implying your subject entrant sampling 
		#' 							assumption is from a superpopulation or "randomization-exact" implying a finite sampling
		#' 							assumption). The default option is "randomization-exact" as it provided properly-sized 
		#' 							tests in our simulations.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference 
		#' 							(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 							for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD")
		#' seq_des$add_subject_to_experiment(c(1, 38, 142, 71, 5.3, 0, 0, 0, 1, 0))
		#' seq_des$add_subject_to_experiment(c(0, 27, 127, 60, 5.5, 0, 0, 0, 1, 0))
		#' seq_des$add_subject_to_experiment(c(1, 42, 169, 74, 5.1, 0, 1, 0, 0, 0))
		#' seq_des$add_subject_to_experiment(c(0, 59, 105, 62, 5.9, 0, 0, 0, 1, 0))
		#' seq_des$add_subject_to_experiment(c(1, 32, 186, 66, 5.6, 1, 0, 0, 0, 0))
		#' seq_des$add_subject_to_experiment(c(1, 37, 178, 75, 6.5, 0, 0, 0, 0, 1))
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInference$new(seq_des)
		#'  
		#' @return A new `SeqDesignTest` object.
		initialize = function(seq_des_obj, estimate_type = "default_regression", test_type = "randomization-exact", num_cores = 1, verbose = TRUE){
			assertClass(seq_des_obj, "SeqDesign")
			seq_des_obj$assert_experiment_completed()
			assertChoice(estimate_type, c("mean_difference", "median_difference", "default_regression"))
			assertChoice(test_type, c("MLE-or-KM-based", "randomization-exact"))
			assertCount(num_cores, positive = TRUE)
			assertFlag(verbose)
			
			self$estimate_type = estimate_type
			self$test_type = test_type
			private$seq_des_obj = seq_des_obj
			private$isKK = seq_des_obj$.__enclos_env__$private$isKK
			private$n = private$seq_des_obj$.__enclos_env__$private$n
			private$p = private$seq_des_obj$.__enclos_env__$private$p
			private$yTs = private$seq_des_obj$y[private$seq_des_obj$w == 1]
			private$yCs = private$seq_des_obj$y[private$seq_des_obj$w == 0]
			private$deadTs = private$seq_des_obj$dead[private$seq_des_obj$w == 1]
			private$deadCs = private$seq_des_obj$dead[private$seq_des_obj$w == 0]
			private$ybarT_minus_ybarC = mean(private$yTs) - mean(private$yCs)	
#			private$has_censoring = sum(private$seq_des_obj$dead) < private$n #if there's censoring, there'd be still some subjects alive at the close of the study
			private$num_cores = num_cores
			private$verbose = verbose
#			private$response_type = private$seq_des_obj$response_type
#			if (private$has_censoring & private$response_type == "survival"){
#				private$response_type = "survival_censored"
#			} else if (!private$has_censoring & private$response_type == "survival"){
#				private$response_type = "survival_uncensored"
#			}
			private$prob_T = private$seq_des_obj$prob_T
			private$prob_T_fifty_fifty = private$prob_T == 0.5 
			
			response_type = private$seq_des_obj$response_type
			private$is_continuous = response_type == "continuous"
			private$is_incidence = response_type == "incidence"
			private$is_count = response_type == "count"
			private$is_proportion = response_type == "proportion"
			private$is_survival = response_type == "survival"
			
			if (estimate_type == "median_difference" & !private$is_survival){
				stop("estimate_type \"median_difference\" is only available for response_type \"survival\"")
			}
			
			if (private$verbose){
				cat(paste0("Intialized inference methods for a ", seq_des_obj$design, " design, response type ", response_type, ", estimation type ", estimate_type, " and test type: ", test_type, ".\n"))
			}	
		},
		
		#' @description
		#' Computes for estimate type "mean_difference-or-medians" either
		#' (1a) for incidence outcomes, the additive log odds treatment effect using logistic regression
		#' (1b) for survival outcomes, the median difference for suvival using the Kaplan-Meier estimates for both arms 
		#' (1c) for count outcomes, the additive treatment effect on log count using negative binomial regression
		#' (1d) for proportion and continous outcomes (where the latter is not under an equal allocation KK design), 
		#' the classic mean_difference estimate of the additive treatment effect, 
		#' (1e) for continuous outcome, equal allocation to arms and KK designs, there's a special match-reservoir weighted 
		#' classic mean_difference estimate
		#' 
		#' Computes for estimte type "default_regression" either
		#' (2a) for incidence outcomes, the additive log odds treatment effect using logistic regression controlled for all other covariates
		#' (2b) for survival outcomes, the additive treatment effect on log suvival using Weibull regression controlled for all other covariates
		#' (2c) for count outcomes, the additive treatment effect on log count using negative binomial regression controlled for all other covariates
		#' (2d) for proportion outcome, the additive treatment effect on proportion using beta regression controlled for all other covariates
		#' (2e) for continous outcomes but not under an equal allocation KK design, the additive treatment effect using OLS regression controlled for all other covariates
		#' (2f) for continuous outcome, equal allocation to arms and KK designs, there's a special match-reservoir weighted 
		#' OLS regression controlled for all other covariates
		#' 
		#' @return 	The setting-appropriate (see description) numeric estimate of the treatment effect
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD", response_type = "continuous")
		#' seq_des$add_subject_to_experiment(c(1, 38, 142, 71, 5.3, 0, 0, 0, 1, 0))
		#' seq_des$add_subject_to_experiment(c(0, 27, 127, 60, 5.5, 0, 0, 0, 1, 0))
		#' seq_des$add_subject_to_experiment(c(1, 42, 169, 74, 5.1, 0, 1, 0, 0, 0))
		#' seq_des$add_subject_to_experiment(c(0, 59, 105, 62, 5.9, 0, 0, 0, 1, 0))
		#' seq_des$add_subject_to_experiment(c(1, 32, 186, 66, 5.6, 1, 0, 0, 0, 0))
		#' seq_des$add_subject_to_experiment(c(1, 37, 178, 75, 6.5, 0, 0, 0, 0, 1))
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInference$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' 		
		compute_treatment_estimate = function(){
			if (is.null(private$beta_T_hat)){				
				private$beta_T_hat = 	if (self$estimate_type == "mean_difference"){
											if (private$isKK & private$prob_T_fifty_fifty & private$is_continuous){
												private$post_matching_data_KK_inference_helper()$beta_hat_T	
											} else if (private$is_incidence){
												private$compute_logistic_model_summary(use_covariates = FALSE)[2, 1]
											} else if (private$is_count){
												private$compute_count_model_summary(use_covariates = FALSE)[2, 1]
											} else if (private$is_survival){
												private$compute_survival_model_summary(use_covariates = FALSE)[2, 1]
											} else {
												private$ybarT_minus_ybarC #faster than calling lm.fit or betareg (and equivalent)
											}											
										} else if (self$estimate_type == "median_difference"){
											private$survival_diff_medians_censored_estimate(private$seq_des_obj$y, private$seq_des_obj$dead, private$seq_des_obj$w)
										} else if (self$estimate_type == "default_regression"){
											if (private$isKK & private$prob_T_fifty_fifty & private$is_continuous){
												private$compute_MLE_based_inference_ols_KK()$beta_hat_T
											} else if (private$is_continuous){
												private$compute_continuous_model_summary(use_covariates = TRUE, estimate_only = TRUE)
											} else if (private$is_incidence){
												private$compute_logistic_model_summary(use_covariates = TRUE)[2, 1]
											} else if (private$is_count){
												private$compute_count_model_summary(use_covariates = TRUE)[2, 1]
											} else if (private$is_proportion){
												private$compute_proportion_model_summary(use_covariates = TRUE)[2, 1]
											} else if (private$is_survival){
												private$compute_survival_model_summary(use_covariates = TRUE)[2, 1]
											}
										}
			}
			private$beta_T_hat
		},
		
		#' @description
		#' Computes for estimate type "mean_difference-or-medians" either
		#' (1a) for incidence outcomes (ignoring the KK design structure), 
		#' the p-value for the test of the additive log odds treatment effect being zero using logistic regression's MLE normal approximation
		#' (1b) for survival outcomes (ignoring the KK design structure), the median difference for survival using the Kaplan-Meier estimates for both arms 
		#' (1c) for count, proportion and continous outcomes (all ignoring the KK design structure), 
		#' the classic mean_difference estimate of the additive treatment effect, 
		#' (1d) for continuous outcome, equal allocation to arms and KK designs, there's a special match-reservoir weighted 
		#' classic mean_difference estimate
		#' 
		#' Computes for estimte type "default_regression" either
		#' (2a) for incidence outcomes, the additive log odds treatment effect using logistic regression controlled for all other covariates
		#' (2b) for survival outcomes, the additive treatment effect on log suvival using Weibull regression controlled for all other covariates
		#' (2c) for count outcomes, the additive treatment effect on log count using negative binomial regression controlled for all other covariates
		#' (2d) for proportion outcome, the additive treatment effect on proportion using beta regression controlled for all other covariates
		#' (2e) for continous outcomes but not under an equal allocation KK design, the additive treatment effect using OLS regression controlled for all other covariates
		#' (2f) for continuous outcome, equal allocation to arms and KK designs, there's a special match-reservoir weighted 
		#' OLS regression controlled for all other covariates
		#' 
		#' 1. classic frequentist confidence interval (CI) of the additive treatment effect
		#' employing the normal theory approximation for both the
		#' (a) difference in means estimator i.e. [ybar_T - ybar_C +/- t_{alpha/2, n_T + n_C - 2} s_{ybar_T - ybar_C}] or
		#' (b) the default_regression estimator i.e. [beta_hat_T +/- t_{alpha/2, n + p - 2} s_{beta_hat_T}]
		#' where the z approximation is employed in lieu of the t is the design is a KK design or
		#' 
		#' 2. a randomization-based CI of an additive shift effect of the potential outcomes under treatment and control
		#' by an inversion of the randomization test at level alpha (this feature is incomplete).
		#' 
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#' @param nsim_exact_test		The number of randomization vectors. 
		#' 								The default is 1000 providing good resolutions to confidence intervals.
		#' @param B						Number of bootstrap samples for the survival response where \code{estimate_type} is "median_difference"
		#' 								(see \link{\code{controlTest::quantileControlTest}}. The default is 501 providing pvalue resolution 
		#' 								to a fifth of a percent.
		#' 
		#' @return 	A 1 - alpha sized frequentist confidence interval for the treatment effect
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD")
		#' seq_des$add_subject_to_experiment(c(1, 38, 142, 71, 5.3, 0, 0, 0, 1, 0))
		#' seq_des$add_subject_to_experiment(c(0, 27, 127, 60, 5.5, 0, 0, 0, 1, 0))
		#' seq_des$add_subject_to_experiment(c(1, 42, 169, 74, 5.1, 0, 1, 0, 0, 0))
		#' seq_des$add_subject_to_experiment(c(0, 59, 105, 62, 5.9, 0, 0, 0, 1, 0))
		#' seq_des$add_subject_to_experiment(c(1, 32, 186, 66, 5.6, 1, 0, 0, 0, 0))
		#' seq_des$add_subject_to_experiment(c(1, 37, 178, 75, 6.5, 0, 0, 0, 0, 1))
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInference$new(seq_des, test_type = "MLE-or-KM-based")
		#' seq_des_inf$compute_confidence_interval()
		#' 		
		compute_confidence_interval = function(alpha = 0.05, nsim_exact_test = 501, B = 501){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertCount(nsim_exact_test, positive = TRUE)
			one_minus_alpha_over_two = 1 - alpha / 2
			z_one_minus_alpha_over_two = qnorm(one_minus_alpha_over_two)
			
			if (test_type == "MLE-or-KM-based"){
				moe = 	if (self$estimate_type == "mean_difference"){
							if (private$isKK & private$prob_T_fifty_fifty & private$is_continuous){
								z_one_minus_alpha_over_two *
									private$post_matching_data_KK_inference_helper()$s_beta_hat_T
							} else if (private$is_incidence){
								z_one_minus_alpha_over_two *
									private$compute_logistic_model_summary(use_covariates = FALSE)[2, 2]
							} else if (private$is_count){
								z_one_minus_alpha_over_two *
									private$compute_count_model_summary(use_covariates = FALSE)[2, 2]
							} else if (private$is_survival){
								z_one_minus_alpha_over_two *
									private$compute_survival_model_summary(use_covariates = FALSE)[2, 2]
							} else { #proportion or continuous
								qt(one_minus_alpha_over_two, length(private$yTs) + length(private$yCs) - 2) * 
									private$simple_t_stat_calculations_for_difference_in_means_estimate()$s_beta_hat_T
							}											
						} else if (self$estimate_type == "median_difference"){
							test_obj = controlTest::quantileControlTest(private$yTs, private$deadTs, private$yCs, private$deadCs, B = B)
							z_one_minus_alpha_over_two *
								test_obj$se
						} else if (self$estimate_type == "default_regression"){
							if (private$isKK & private$prob_T_fifty_fifty & private$is_continuous){
								z_one_minus_alpha_over_two *
									private$compute_MLE_based_inference_ols_KK()$s_beta_hat_T
							} else if (private$is_continuous){
								qt(one_minus_alpha_over_two, private$n - private$p - 2) * #subtract two for intercept and allocation vector 
									private$compute_continuous_model_summary(use_covariates = TRUE)[2, 2]
							} else if (private$is_incidence){
								z_one_minus_alpha_over_two *
									private$compute_logistic_model_summary(use_covariates = TRUE)[2, 2]
							} else if (private$is_count){
								z_one_minus_alpha_over_two *
									private$compute_count_model_summary(use_covariates = TRUE)[2, 2]
							} else if (private$is_proportion){
								z_one_minus_alpha_over_two *
									private$compute_count_model_summary(use_covariates = TRUE)[2, 2]
							} else if (private$is_survival){
								z_one_minus_alpha_over_two *
									private$compute_survival_model_summary(use_covariates = TRUE)[2, 2]
							}
						}
				#return the CI
				private$beta_T_hat + c(-moe, moe)
			} else { #randomization
				assertCount(nsim_exact_test, positive = TRUE)
				c(
					private$compute_ci_lower_by_inverting_the_randomization_test(nsim_exact_test, private$num_cores),
					private$compute_ci_upper_by_inverting_the_randomization_test(nsim_exact_test, private$num_cores)
				)
			}
			
#			if (self$estimate_type == "mean_difference-or-medians" & self$test_type == "MLE-or-KM-based" & private$isKK){
#				private$common_MLE_based_ci(private$compute_MLE_or_KM_based_inference_difference_in_means_or_medians_KK(), alpha, use_Z = TRUE)
#			} else if (self$estimate_type == "mean_difference-or-medians" & self$test_type == "MLE-or-KM-based" & !private$isKK){
#				private$common_MLE_based_ci(private$compute_MLE_or_KM_based_inference_difference_in_means_or_medians(), alpha, use_Z = FALSE)
#			} else if (self$estimate_type == "default_regression" & self$test_type == "MLE-or-KM-based" & private$isKK){
#				private$common_MLE_based_ci(private$compute_MLE_based_inference_ols_KK(), alpha, use_Z = TRUE)	
#			} else if (self$estimate_type == "default_regression" & self$test_type == "MLE-or-KM-based" & !private$isKK){
#				private$common_MLE_based_ci(private$compute_MLE_based_inference_ols(), alpha, use_Z = FALSE)
#			} else if (self$estimate_type == "mean_difference-or-medians" & self$test_type == "randomization-exact"){
#				stop("This computation is not implemented yet.")
#				#TODO
#			} else if (self$estimate_type == "default_regression" & self$test_type == "randomization-exact"){
#				stop("This computation is not implemented yet.")
#				#TODO
#			}
		},
		
		#' @description
		#' Computes a 2-sided p-value for all types of inferential settings written about in the initializer
		#' (1) estimate type "mean_difference-or-medians" and test type "MLE-or-KM-based"
		#' This implies the classic mean_difference estimator which means that 
		#'   (a) For continous and proportion outcomes, H_0: E[Y_T] - E[Y_C] = delta,
		#'   (b) For incidence outcomes, H_0: log(Odds(P(Y_T = 1)) - log(Odds(P(Y_C = 1) = delta,
		#'   (c) For count outcomes, H_0: E[ln(Y_T)] - E[ln(Y_C)] = delta or
		#'   (d) For survival outcomes, H_0: MED[Y_T] - MED[Y_C] = delta
		#' (2) Fisher's randomization test which means that H_0: y_i_T - y_i_C = delta for all subjects
		#' either the classic different-in-means estimate of the additive treatment effect, 
		#' i.e. ybar_T - ybar_C or the default_regression estimate of the additive treatment effect linearly i.e. 
		#' the treatment different adjusted linearly for the p covariates.
		#' @param nsim_exact_test		The number of randomization vectors to use in the randomization test (ignored if \code{test_type}
		#' 								is not "randomization-exact"). The default is 501 providing pvalue resolution to a fifth of a percent.
		#' @param B						Number of bootstrap samples for the survival response where \code{estimate_type} is "median_difference"
		#' 								(see \link{\code{controlTest::quantileControlTest}}. The default is 501 providing pvalue resolution 
		#' 								to a fifth of a percent.
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#' 
		#' @return 	The approximate frequentist p-value
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD")
		#' seq_des$add_subject_to_experiment(c(1, 38, 142, 71, 5.3, 0, 0, 0, 1, 0))
		#' seq_des$add_subject_to_experiment(c(0, 27, 127, 60, 5.5, 0, 0, 0, 1, 0))
		#' seq_des$add_subject_to_experiment(c(1, 42, 169, 74, 5.1, 0, 1, 0, 0, 0))
		#' seq_des$add_subject_to_experiment(c(0, 59, 105, 62, 5.9, 0, 0, 0, 1, 0))
		#' seq_des$add_subject_to_experiment(c(1, 32, 186, 66, 5.6, 1, 0, 0, 0, 0))
		#' seq_des$add_subject_to_experiment(c(1, 37, 178, 75, 6.5, 0, 0, 0, 0, 1))
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInference$new(seq_des)
		#' seq_des_inf$compute_two_sided_pval_for_treatment_effect()
		#' 		
		compute_two_sided_pval_for_treatment_effect = function(nsim_exact_test = 501, B = 501, delta = 0){
			assertNumeric(delta)
			if (delta != 0){
				stop("nonzero treatment effect tests not yet supported")
			}
			assertCount(B, positive = TRUE)
			
			if (test_type == "MLE-or-KM-based"){
				if (self$estimate_type == "mean_difference"){
					if (private$isKK & private$prob_T_fifty_fifty & private$is_continuous){
						private$post_matching_data_KK_inference_helper()$p_val	
					} else if (private$is_incidence){
						private$compute_logistic_model_summary(use_covariates = FALSE)[2, 4]
					} else if (private$is_count){
						private$compute_count_model_summary(use_covariates = FALSE)[2, 4]
					} else if (private$is_survival){
						private$compute_survival_model_summary(use_covariates = FALSE)[2, 4]
					} else {
						private$simple_t_stat_calculations_for_difference_in_means_estimate()$p_val
					}											
				} else if (self$estimate_type == "median_difference"){
					controlTest::quantileControlTest(private$yTs, private$deadTs, private$yCs, private$deadCs, B = B)$pval
				} else if (self$estimate_type == "default_regression"){
					if (private$isKK & private$prob_T_fifty_fifty & private$is_continuous){
						private$compute_MLE_based_inference_ols_KK()$p_val
					} else if (private$is_continuous){
						private$compute_continuous_model_summary(use_covariates = TRUE)[2, 4]
					} else if (private$is_incidence){
						private$compute_logistic_model_summary(use_covariates = TRUE)[2, 4]
					} else if (private$is_count){
						private$compute_count_model_summary(use_covariates = TRUE)[2, 4]
					} else if (private$is_proportion){
						private$compute_proportion_model_summary(use_covariates = TRUE)[2, 4]
					} else if (private$is_survival){
						private$compute_survival_model_summary(use_covariates = TRUE)[2, 4]
					}
				}	
			} else { #randomization
				assertCount(nsim_exact_test, positive = TRUE)
				private$compute_randomization_test_p_val(nsim_exact_test, private$num_cores)				
			}
		}
	),
	
	private = list(
		seq_des_obj = NULL,
		isKK = FALSE,
		rand_inf_b_T_sims = NULL,
		num_cores = NULL,
		verbose = FALSE,
		n = NULL,
		p = NULL,
		yTs = NULL,
		yCs = NULL,
		deadTs = NULL,
		deadCs = NULL,
		ybarT_minus_ybarC = NULL,
#		has_censoring = NULL,
		beta_T_hat = NULL,
		prob_T = NULL,
		prob_T_fifty_fifty = NULL,
		KKstats = NULL,
		
		is_continuous = FALSE,
		is_incidence = FALSE,
		is_count = FALSE,
		is_proportion = FALSE,
		is_survival = FALSE,
		
		############# TESTING HELPER FUNCTIONS
		
		compute_ci_lower_by_inverting_the_randomization_test = function(nsim_exact_test, num_cores){
			NA
		},
		compute_ci_upper_by_inverting_the_randomization_test = function(nsim_exact_test, num_cores){
			NA
		},
			
		compute_randomization_test_p_val = function(nsim_exact_test, num_cores, delta = 0){
			if (is.null(private$rand_inf_b_T_sims)){
				cl = makeCluster(private$num_cores)
				#ensure to get all variables copied into scope (saves cycles)
				X = private$seq_des_obj$X
				y = private$seq_des_obj$y
				if (delta != 0){
					#we are testing against H_0: y_T_i - y_C_i = delta <=> (y_T_i - delta) - y_C_i = 0
					y[private$seq_des_obj$w == 1] = y[private$seq_des_obj$w == 1] - delta
				}
				dead = private$seq_des_obj$dead
				design = private$seq_des_obj$design
				estimate_type = self$estimate_type
				if (private$num_cores == 1){
					#easier on the OS I think...
					b_T_sims = array(NA, nsim_exact_test)
					for (r in 1 : nsim_exact_test){
						#cat("		r =", r, "/", nsim_exact_test, "\n")
						#do a fake run of the experiment to come up with a different w based on the design
						seq_des_r = SeqDesign$new(n = private$n, p = private$p, design = design, response_type = private$seq_des_obj$response_type, verbose = FALSE)
						for (t in 1 : private$n){
							seq_des_r$add_subject_to_experiment(X[t, ])
							seq_des_r$add_subject_response(t, y[t])
						}				
						b_T_sims[r] = SeqDesignInference$new(seq_des_r, estimate_type = estimate_type, verbose = FALSE)$compute_treatment_estimate()
					}
				} else {
					registerDoParallel(cl)			
					#now copy them to each core's memory
					clusterExport(cl, list("X", "y", "dead", "n", "p", "design", "response_type", "estimate_type"), envir = environment())
					#now do the parallelization
					b_T_sims = foreach(r = 1 : nsim_exact_test, .inorder = FALSE, .combine = c) %dopar% {
					#cat("		r =", r, "/", nsim_exact_test, "\n")
						#do a fake run of the experiment to come up with a different w based on the design
						seq_des_r = SeqDesign$new(n = n, p = p, design = design, response_type = response_type, verbose = FALSE)
						for (t in 1 : n){
							seq_des_r$add_subject_to_experiment(X[t, ])
							seq_des_r$add_subject_response(t, y[t], dead[t])
						}				
						SeqDesignInference$new(seq_des_r, estimate_type = estimate_type, verbose = FALSE)$compute_treatment_estimate()				
					}
					stopCluster(cl)
				}
				private$rand_inf_b_T_sims = b_T_sims
			}
			sum(abs(self$compute_treatment_estimate()) < abs(private$rand_inf_b_T_sims)) / nsim_exact_test		
		},
		
		############# CONFIDENCE INTERVALS HELPER FUNCTIONS
		
#		common_MLE_based_ci = function(inference_obj, alpha, use_Z){
#			qu = 1 - alpha / 2
#			moe = ifelse(use_Z, qnorm(qu), qt(qu, inference_obj$t_df)) * inference_obj$s_beta_hat_T
#			inference_obj$beta_hat_T + c(-moe, moe)
#		},
		
		######### ALL INFERENCE HELPER FUNCTIONS	
		
		
		
		simple_t_stat_calculations_for_difference_in_means_estimate = function(){
			nT = length(private$yTs)
			nC = length(private$yCs)
			s_beta_hat_T = sqrt(var(private$yTs) / nT + var(private$yCs) / nC)
			list(
					s_beta_hat_T = s_beta_hat_T, 
					p_val = 2 * pt(-abs(private$ybarT_minus_ybarC / s_beta_hat_T), nT + nC - 2)	
			)			
		},
		
		compute_continuous_model_summary = function(use_covariates, estimate_only = FALSE){
			if (estimate_only){
				ols_regr_mod = 	if (use_covariates){
									lm.fit(cbind(1, private$seq_des_obj$w, private$seq_des_obj$X), private$seq_des_obj$y)
								} else {
									lm.fit(cbind(1, private$seq_des_obj$w), private$seq_des_obj$y)
								}
				ols_regr_mod$coefficients[2]
			} else {
				ols_regr_mod = if (use_covariates){
									lm(private$seq_des_obj$y ~ ., data = cbind(data.frame(w = private$seq_des_obj$w), private$seq_des_obj$X))
								} else {
									lm(private$seq_des_obj$y ~ private$seq_des_obj$w)
								}			
				coef(summary(ols_regr_mod))
			}			
		},
		
		compute_logistic_model_summary = function(use_covariates, estimate_only = FALSE){
			if (estimate_only){
				stop("this speedup is not implemented yet... use estimate_only = FALSE")
			} else {	
				logistic_regr_mod = if (use_covariates){
										suppressWarnings(glm(private$seq_des_obj$y ~ cbind(private$seq_des_obj$w, private$seq_des_obj$X), family = "binomial"))
									} else {
										suppressWarnings(glm(private$seq_des_obj$y ~ private$seq_des_obj$w, family = "binomial"))
									}
				coef(summary(logistic_regr_mod))
			}
		},
		
		compute_count_model_summary = function(use_covariates, estimate_only = FALSE){
			if (estimate_only){
				stop("this speedup is not implemented yet... use estimate_only = FALSE")
			} else {
				negbin_regr_mod = 	if (use_covariates){
										suppressWarnings(MASS::glm.nb(y ~ ., data = cbind(data.frame(y = private$seq_des_obj$y, w = private$seq_des_obj$w), private$seq_des_obj$X)))									
									} else {
										suppressWarnings(MASS::glm.nb(y ~ ., data = data.frame(y = private$seq_des_obj$y, w = private$seq_des_obj$w)))
									}			
				coef(summary(negbin_regr_mod))
			}
		},
		
		compute_proportion_model_summary = function(use_covariates, estimate_only = FALSE){
			if (estimate_only){
				stop("this speedup is not implemented yet... use estimate_only = FALSE")
			} else {			
				beta_regr_mod = if (use_covariates){
									suppressWarnings(betareg::betareg(y ~ ., data = cbind(data.frame(y = private$seq_des_obj$y, w = private$seq_des_obj$w), private$seq_des_obj$X)))
								} else {
									suppressWarnings(betareg::betareg(y ~ ., data = data.frame(y = private$seq_des_obj$y, w = private$seq_des_obj$w)))
								}			
				coef(summary(beta_regr_mod))$mean
			}
		},
		
		compute_survival_model_summary = function(use_covariates, estimate_only = FALSE){
			if (estimate_only){
				stop("this speedup is not implemented yet... use estimate_only = FALSE")
			} else {			
				surv_obj = survival::Surv(private$seq_des_obj$y, private$seq_des_obj$dead)
				#warning!!!! sometimes weibull fails to converge and blows up here!!!
				surv_regr_mod = if (use_covariates){
									suppressWarnings(survival::survreg(surv_obj ~ ., data = cbind(data.frame(w = private$seq_des_obj$w), private$seq_des_obj$X), dist = "weibull"))
								} else {
									suppressWarnings(survival::survreg(surv_obj ~ ., data = data.frame(w = private$seq_des_obj$w), dist = "weibull"))
								}
				summary(surv_regr_mod)$table	
			}
		},
		
		survival_diff_medians_censored_estimate = function(y, dead, w){
			survival_obj = survival::Surv(y, dead)
			survival_fit_obj = survival::survfit(survival_obj ~ w) 
			survival_fit_res = summary(survival_fit_obj)$table
			survival_fit_res[2, 7] - survival_fit_res[1, 7]	
		},
		
		## compute_MLE_or_KM_based_inference_difference_in_means_or_medians = function(estimate_only = FALSE){
		## 
		##     if (private$seq_des_obj$response_type == "survival_censored"){
		##         survival_obj = survival::Surv(private$seq_des_obj$y, private$seq_des_obj$dead)
		##     }
		##     beta_hat_T = 	if (private$seq_des_obj$response_type == "survival_censored"){								
		##                         survival_fit_obj = survival::survfit2(survival_obj ~ private$seq_des_obj$w) 
		##                         survival_fit_res = summary(survival_fit_obj)$table
		##                         survival_fit_res[2, 7] - survival_fit_res[1, 7]
		##                     } else {
		##                         mean(yTs) - mean(yCs)
		##                     }
		## 
		##     if (estimate_only){
		##         list(beta_hat_T = beta_hat_T)
		##     } else if (private$seq_des_obj$response_type == "survival_censored"){
		##         #https://cran.r-project.org/web/packages/bpcp/index.html
		##         list(
		##             beta_hat_T = beta_hat_T,
		##             s_beta_hat_T = NA,
		##             t_df = NA,
		##             p_val = survival::survdiff(survival_obj ~ private$seq_des_obj$w)$pvalue
		##         )
		##     } else {
		##         nT = length(yTs)
		##         nC = length(yCs)
		##         s_beta_hat_T = sqrt(var(yTs) / nT + var(yCs) / nC)
		##         t_df = nT + nC - 2	
		##         list(
		##             beta_hat_T = beta_hat_T,
		##             s_beta_hat_T = s_beta_hat_T,
		##             t_df = t_df,
		##             p_val = 2 * pt(-abs(beta_hat_T / s_beta_hat_T), t_df)
		##         )			
		##     }			
		## },
		
#		compute_MLE_based_inference_difference_in_means_KK = function(){
#			KKstats = private$compute_post_matching_data_KK()
#			KKstats
#		},
		
#		compute_MLE_based_inference_ols = function(estimate_only = FALSE){
#			if (estimate_only){
#
#			} else {
#				mod_results = coef(summary(lm(private$seq_des_obj$y ~ ., data.frame(cbind(private$seq_des_obj$w, private$seq_des_obj$X)))))
#				list(
#					beta_hat_T = mod_results[2, 1],
#					s_beta_hat_T = mod_results[2, 2],
#					t_df = private$seq_des_obj$.__enclos_env__$private$n - private$seq_des_obj$.__enclos_env__$private$p - 1,
#					p_val = mod_results[2, 4]
#				)				
#			}			
#		},
		
		compute_MLE_based_inference_ols_KK = function(){
			KKstats = private$compute_post_matching_data_KK()
			
			if (KKstats$nRT <= 2 || KKstats$nRC <= 2 || (KKstats$nRT + KKstats$nRC <= private$seq_des_obj$.__enclos_env__$private$p + 2)){
				coefs_matched = coef(summary(lm(KKstats$y_matched_diffs ~ KKstats$X_matched_diffs)))
				
				KKstats$beta_hat_T = coefs_matched[1, 1]
				KKstats$s_beta_hat_T = coefs_matched[1, 2]
				KKstats$p_val = coefs_matched[1, 4]
				
			#and sometimes there's no matches	
			} else if (KKstats$m == 0){			
				coefs_reservoir = coef(summary(lm(KKstats$y_reservoir ~ cbind(KKstats$w_reservoir, KKstats$X_reservoir))))		
				KKstats$beta_hat_T = coefs_reservoir[2, 1]
				KKstats$s_beta_hat_T = coefs[2, 2]
				KKstats$p_val = coefs[2, 4]
				
			#but most of the time... we have matches and a nice-sized reservoir
			} else {
				#compute estimator from matched pairs by regression
				coefs_matched = coef(summary(lm(KKstats$y_matched_diffs ~ KKstats$X_matched_diffs)))
				beta_match_regression = coefs_matched[1, 1]
				ssqd_match_regression = coefs_matched[1, 2]^2 #lin mod returns SE not VAR, so square it				
				
				#compute estimator reservoir sample std error
				coefs_reservoir = coef(summary(lm(KKstats$y_reservoir ~ cbind(KKstats$w_reservoir, KKstats$X_reservoir))))
				beta_reservoir_regression = coefs_reservoir[2, 1]
				ssqd_reservoir_regression = coefs_reservoir[2, 2]^2 #lin mod returns SE not VAR, so square it
				
				w_star = ssqd_reservoir_regression / (ssqd_reservoir_regression + ssqd_match_regression) #just a convenience for faster runtime	
				KKstats$beta_hat_T = w_star * beta_match_regression + (1 - w_star) * beta_reservoir_regression #proper weighting
				KKstats$s_beta_hat_T = sqrt(ssqd_match_regression * ssqd_reservoir_regression / (ssqd_match_regression + ssqd_reservoir_regression)) #analagous eq's
				KKstats$p_val = 2 * (pnorm(-abs(KKstats$beta_hat_T / KKstats$s_beta_hat_T))) #approximate by using N(0, 1) distribution
			}
			KKstats
		},
		
		post_matching_data_KK_inference_helper = function(){
			KKstats = private$compute_post_matching_data_KK()
			#sometimes the reservoir just isn't large enough
			if (KKstats$nRT <= 1 || KKstats$nRC <= 1){
				KKstats$beta_hat_T = KKstats$d_bar	
				KKstats$s_beta_hat_T = sqrt(KKstats$ssqD_bar)
			} else if (KKstats$m == 0){ #sometimes there's no matches
				KKstats$beta_hat_T = KKstats$r_bar		
				KKstats$s_beta_hat_T = sqrt(KKstats$ssqR)		
			} else {
				KKstats$beta_hat_T = KKstats$w_star * KKstats$d_bar + (1 - KKstats$w_star) * KKstats$r_bar #proper weighting
				KKstats$s_beta_hat_T = sqrt(KKstats$ssqR * KKstats$ssqD_bar / (KKstats$ssqR + KKstats$ssqD_bar))
			}
			KKstats$p_val = 2 * (pnorm(-abs(KKstats$beta_hat_T / KKstats$s_beta_hat_T))) #approximate by using real Z
			KKstats
		},
		
		compute_post_matching_data_KK = function(){	
			if (is.null(private$KKstats)){
				#get matched data
				match_indic = private$seq_des_obj$.__enclos_env__$private$match_indic
				m = max(match_indic)
				y_matched_diffs = array(NA, m)
				X_matched_diffs = matrix(NA, nrow = m, ncol = ncol(private$seq_des_obj$X))
				if (m > 0){
					for (match_id in 1 : m){ #we want to just calculate the diffs inside matches and ignore the reservoir
						yTs = private$seq_des_obj$y[private$seq_des_obj$w == 1 & match_indic == match_id]
						yCs = private$seq_des_obj$y[private$seq_des_obj$w == 0 & match_indic == match_id]
						y_matched_diffs[match_id] = mean(yTs) - mean(yCs)
						
						xmTs = private$seq_des_obj$X[private$seq_des_obj$w == 1 & match_indic == match_id, ]
						xmCs = private$seq_des_obj$X[private$seq_des_obj$w == 0 & match_indic == match_id, ]
						X_matched_diffs[match_id, ] = mean(xmTs) - mean(xmCs)
					}
				}			
				
				#get reservoir data
				X_reservoir = 	private$seq_des_obj$X[match_indic == 0, ]
				y_reservoir = 	private$seq_des_obj$y[match_indic == 0]
				w_reservoir = 	private$seq_des_obj$w[match_indic == 0]
				y_reservoir_T = y_reservoir[w_reservoir == 1] #get the reservoir responses from the treatment
				y_reservoir_C = y_reservoir[w_reservoir == 0] #get the reservoir responses from the control
				r_bar = mean(y_reservoir_T) - mean(y_reservoir_C) #compute the classic estimator from the reservoir: ybar_T - ybar_C
				
				#get reservoir sample sizes
				nRT = length(y_reservoir_T) #how many treatment observations are there in the reservoir?
				nRC = length(y_reservoir_C) #how many control observations are there in the reservoir?
				nR = nRT + nRC #how many observations are there in the reservoir?
				
				ssqR = (var(y_reservoir_T) * (nRT - 1) + var(y_reservoir_C) * (nRC - 1)) / (nR - 2) * (1 / nRT + 1 / nRC)
				ssqD_bar = var(y_matched_diffs) / m
				
				w_star = ssqR / (ssqR + ssqD_bar)
				
				private$KKstats = list(
					X_matched_diffs = X_matched_diffs,
					y_matched_diffs = y_matched_diffs,
					X_reservoir = X_reservoir,
					y_reservoir = y_reservoir,
					w_reservoir = w_reservoir,
					nRT = nRT,
					nRC = nRC,
					m = m,
					d_bar = mean(y_matched_diffs),
					ssqD_bar = ssqD_bar,
					r_bar = r_bar,
					ssqR = ssqR,
					w_star = w_star
				)
			}
			private$KKstats
		}
	)
)