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
		#' @field estimate_type		The estimate type (see initializer documentation).
		estimate_type = NULL,		
		#' 
		#' @field test_type			The type of test to run (see initializer documentation).
		test_type = NULL,
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' 
		#' 
		#' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param estimate_type		The type of estimate to compute of which there are many and identified by the response type as its first word (except
		#' 							\code{simple_mean_difference} and \code{KK_compound_mean_difference} which are available for all response types except
		#' 							survival when there are any censored observations). If the string "KK" 
		#' 							appears after the first word, then this estimate type is only applicable to KK14, KK21, KK21stepwise designs.
		#' 							  * "simple_mean_difference"
		#' 								assumes the treatment effect parameter is an additive treatment effect 
		#' 							  * "KK_compound_mean_difference" 	
		#' 								assumes the treatment effect parameter is an additive treatment effect	
		#' 								and estimates via combining a simple average difference estimator for both the matches and the reservoir
		#' 								and estimates via the simple average difference
		#' 							  * "continuous_regression_with_covariates"
		#' 								assumes the treatment effect parameter is an additive treatment effect	
		#' 								and the presence of linear additive covariates 
		#' 								and estimates via OLS 
		#' 							  * "continuous_KK_compound_multivariate_regression"	
		#' 								assumes the treatment effect parameter is an additive treatment effect
		#' 								and estimates via combining an OLS estimator for bothe ther matches and the reservoir
		#' 							  * "continuous_KK_regression_with_covariates_with_matching_dummies"
		#' 								assumes the treatment effect parameter is an additive treatment effect 
		#' 								and the presence of linear additive covariates treating the match ID as a factor and estimates via OLS (not recommended)	
		#' 							  * "continuous_KK_regression_with_covariates_with_random_intercepts"	
		#' 								assumes the treatment effect parameter is an additive treatment effect 
		#' 								and the presence of linear additive covariates and random intercepts on the match ID 
		#' 								and estimates via restricted maximum likelihood
		#' 							  * "incidence_simple_log_odds"
		#' 								assumes the treatment effect parameter is additive in the log odds probability of the positive class
		#' 								and estimates via maximum likelihood 	
		#' 							  * "incidence_logistic_regression"
		#' 								assumes the treatment effect parameter is additive in the log odds probability of the positive class
		#' 								and the presence of linear additive covariates also in the log odds probability of the positive class
		#' 								and estimates via maximum likelihood
		#' 							  * "incidence_KK_compound_multivariate_logistic_regression"	
		#' 								assumes the treatment effect parameter is additive in the log odds probability of the positive class
		#' 								and the presence of linear additive covariates treating the match ID as a factor also in the log odds probability of the positive class
		#' 								and estimates via maximum likelihood
		#' 							  * "incidence_KK_multivariate_logistic_regression_with_matching_dummies"	
		#' 								assumes the treatment effect parameter is additive in the log odds probability of the positive class
		#' 								and the presence of linear additive covariates treating the match ID as a factor also in the log odds probability of the positive class
		#' 								and estimates via maximum likelihood
		#' 							  * "incidence_KK_compound_multivariate_logistic_regression_with_random_intercepts_for_matches"
		#' 								assumes the treatment effect parameter is additive in the log odds probability of the positive class
		#' 								and the presence of linear additive covariates 
		#' 								and random intercepts on the match ID also in units of log odds probability of the positive class
		#' 								and estimates via restricted maximum likelihood
		#' 							  * "proportion_simple_logodds_regression"	
		#' 								assumes the treatment effect parameter is additive in the log odds proportion
		#' 								and estimates via beta regression
		#' 							  * "proportion_beta_regression"
		#' 								assumes the treatment effect parameter is additive in the log odds proportion
		#' 								and the presence of linear additive covariates 
		#' 								and estimates via beta regression 
		#'							  * "proportion_KK_compound_univariate_beta_regression"
		#' 								assumes the treatment effect parameter is an additive treatment effect in log odds of proportion
		#' 								and the presence of linear additive covariates also in the log odds of proportion
		#' 								and estimates via combining a simple average difference estimator for both the matches and the reservoir
		#'							  *	"proportion_KK_compound_multivariate_beta_regression"
		#' 								assumes the treatment effect parameter is an additive treatment effect in log odds	
		#' 								and estimates via combining a simple average difference estimator for both the matches and the reservoir
		#' 							  * "proportion_KK_multivariate_beta_regression_with_matching_dummies"
		#' 								assumes the treatment effect parameter is additive in the log odds proportion
		#' 								and the presence of linear additive covariates 
		#' 								and estimates via beta regression
		#' 							  * "count_univariate_negative_binomial_regression"		
		#' 								assumes the treatment effect parameter is additive in the log count
		#' 								and estimates via negative binomial regression
		#' 							  * "count_multivariate_negative_binomial_regression"	
		#' 								assumes the treatment effect parameter is additive in the log count
		#' 								and the presence of linear additive covariates 
		#' 								and estimates via negative binomial regression
		#' 							  * "count_KK_compound_univariate_negative_binomial_regression"	
		#' 								assumes the treatment effect parameter is additive in the log count
		#' 								and treating the match ID as a factor
		#' 								and estimates via maximum likelihood	
		#' 							  * "count_KK_multivariate_negative_binomial_regression_with_matching_dummies"	
		#' 								assumes the treatment effect parameter is additive in the log count
		#' 								and the presence of linear additive covariates
		#' 								and treating the match ID as a factor
		#' 								and estimates via maximum likelihood
		#' 							  * "count_KK_multivariate_negative_binomial_regression_with_random_intercepts_for_matches"
		#' 								assumes the treatment effect parameter is additive in the log count
		#' 								and the presence of linear additive covariates in units of log count
		#' 								and random intercepts on the match ID in the log count
		#' 								and estimates via maximum likelihood							
		#' 							  * "survival_simple_median_difference"	
		#' 								assumes the treatment effect parameter is the difference in survival medians
		#' 								and estimates via Kaplan-Meier
		#' 							  * "survival_simple_restricted_mean_difference"
		#' 								assumes the treatment effect parameter is the difference in survival means
		#' 								and estimates via restricted means (assuming the largest survival time is the absolute limit)
		#' 							  * "survival_univariate_weibull_regression"	
		#' 								assumes the treatment effect parameter is the additive mean survival difference
		#' 								and estimates via Weibull regression
		#' 							  * "survival_multivariate_weibull_regression"
		#' 								assumes the treatment effect parameter is the additive mean survival difference
		#' 								and the presence of linear additive covariates
		#' 								and estimates via Weibull regression	
		#' 							  * "survival_KK_multivariate_weibull_regression_with_matching_dummies"
		#' 								assumes the treatment effect parameter is the additive mean survival difference
		#' 								and the presence of linear additive covariates
		#' 								and treating the match ID as a factor
		#' 								and estimates via Weibull regression	
		#'           				  * "survival_univariate_coxph_regression"
		#' 								assumes the treatment effect is a log difference in hazard which is constant conditional on covariate values
		#' 								and estimates via maximum likelihood
		#'  						  * "survival_multivariate_coxph_regression"
		#' 								assumes the treatment effect is a log difference in hazard which is constant conditional on covariate values
		#' 								and the presence of linear additive covariates in log hazard
		#' 								and estimates via maximum likelihood		
		#'          				  * "survival_KK_multivariate_coxph_regression_with_matching_dummies"
		#' 								assumes the treatment effect is a log difference in hazard which is constant conditional on covariate values
		#' 								and the presence of linear additive covariates in log hazard
		#' 								and treating the match ID as a factor
		#' 								and estimates via maximum likelihood
		#'          				  * "survival_KK_multivariate_coxph_regression_with_random_intercepts_for_matches"
		#' 								assumes the treatment effect is a log difference in hazard which is constant conditional on covariate values
		#' 								and the presence of linear additive covariates in log hazard
		#' 								and random intercepts on the match ID in units of log hazard
		#' 								and estimates via maximum likelihood	
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
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInference$new(seq_des)
		#'  
		#' @return A new `SeqDesignTest` object.
		initialize = function(seq_des_obj, estimate_type, test_type = "randomization-exact", num_cores = 1, verbose = TRUE){
			assertClass(seq_des_obj, "SeqDesign")
			seq_des_obj$assert_experiment_completed()
			assertChoice(estimate_type, c(
				########################################### ALL RESPONSE TYPES
				"simple_mean_difference",				
				"KK_compound_mean_difference",  	
				########################################### CONTINUOUS
				"continuous_regression_with_covariates",
				"continuous_KK_compound_multivariate_regression",
				"continuous_KK_regression_with_covariates_with_matching_dummies",
				"continuous_KK_regression_with_covariates_with_random_intercepts",
				########################################### INCIDENCE
				"incidence_simple_log_odds",	
				"incidence_logistic_regression",
				"incidence_KK_compound_univariate_logistic_regression",
				"incidence_KK_compound_multivariate_logistic_regression",	
				"incidence_KK_multivariate_logistic_regression_with_matching_dummies",	
				"incidence_KK_multivariate_logistic_regression_with_random_intercepts_for_matches",
				########################################### PROPORTION
				"proportion_simple_logodds_regression",
				"proportion_beta_regression",
				"proportion_KK_compound_univariate_beta_regression",
				"proportion_KK_compound_multivariate_beta_regression",
				"proportion_KK_multivariate_beta_regression_with_matching_dummies",
				########################################### COUNT
				"count_univariate_negative_binomial_regression",
				"count_multivariate_negative_binomial_regression",
				"count_KK_compound_univariate_negative_binomial_regression",	
				"count_KK_compound_multivariate_negative_binomial_regression",
				"count_KK_multivariate_negative_binomial_regression_with_matching_dummies",
				"count_KK_multivariate_negative_binomial_regression_with_random_intercepts_for_matches",
				########################################### SURVIVAL
				"survival_simple_median_difference",	
				"survival_simple_restricted_mean_difference",
				"survival_univariate_weibull_regression",	
				"survival_multivariate_weibull_regression",
				"survival_KK_compound_univariate_weibull_regression",	
				"survival_KK_compound_multivariate_weibull_regression",
				"survival_KK_multivariate_weibull_regression_with_matching_dummies",
				"survival_univariate_coxph_regression",	
				"survival_multivariate_coxph_regression",		
				"survival_KK_multivariate_coxph_regression_with_matching_dummies",		
				"survival_KK_multivariate_coxph_regression_with_random_intercepts_for_matches"	
				
			))
			assertChoice(test_type, c("MLE-or-KM-based", "randomization-exact"))
			assertCount(num_cores, positive = TRUE)
			assertFlag(verbose)
			
			if (seq_des_obj$design %in% c("CRD", "iBCRD", "Efron")){ #imputations were never done yet
				seq_des_obj$.__enclos_env__$private$covariate_impute_if_necessary_and_then_create_model_matrix()
			}
			self$estimate_type = estimate_type
			self$test_type = test_type
			private$seq_des_obj = seq_des_obj
			private$isKK = seq_des_obj$.__enclos_env__$private$isKK
			private$match_indic = private$seq_des_obj$.__enclos_env__$private$match_indic
			if (!is.null(private$match_indic)){
				mm = model.matrix(~ 0 + factor(private$match_indic)) 
				mm = mm[, 2 : (ncol(mm) - 1)]
				private$match_indic_model_matrix = mm
			}
			private$n = private$seq_des_obj$.__enclos_env__$private$n
			private$X = private$seq_des_obj$.__enclos_env__$private$compute_all_subject_data()$X_all
			private$yTs = private$seq_des_obj$y[private$seq_des_obj$w == 1]
			private$yCs = private$seq_des_obj$y[private$seq_des_obj$w == 0]
			private$deadTs = private$seq_des_obj$dead[private$seq_des_obj$w == 1]
			private$deadCs = private$seq_des_obj$dead[private$seq_des_obj$w == 0]
			private$num_cores = num_cores
			private$verbose = verbose
			private$prob_T = private$seq_des_obj$prob_T
			private$prob_T_fifty_fifty = private$prob_T == 0.5 
			
			response_type = private$seq_des_obj$response_type
			private$is_continuous = response_type == "continuous"
			private$is_incidence = response_type == "incidence"
			private$is_count = response_type == "count"
			private$is_proportion = response_type == "proportion"
			private$is_survival = response_type == "survival"
			
			#now run more checks
			if (!(estimate_type %in% c("simple_mean_difference", "KK_compound_mean_difference")) & strsplit(estimate_type, "_")[[1]][1] != response_type){
				stop(paste(estimate_type, "cannot be estimated for response type:", response_type))
			}
			if (estimate_type == "simple_mean_difference" & private$is_survival & any(private$seq_des_obj$dead == 0)){
				stop("estimate type \"simple_mean_difference\" cannot be estimated when there is censoring")
			}
			if (estimate_type == "KK_compound_mean_difference" & private$is_survival & any(private$seq_des_obj$dead == 0)){
				stop("estimate type \"KK_compound_mean_difference\" cannot be estimated when there is censoring")
			}
			if (grepl("KK", estimate_type) & !private$isKK){
				stop(paste("estimate type", estimate_type, "cannot be estimated for design type", seq_des_obj$design_type), "(as it is not a KK-type design)")
			}			
			if (estimate_type == "median_difference" & !private$is_survival){
				stop("estimate type \"median_difference\" is only available for response_type \"survival\"")
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
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInference$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' 		
		compute_treatment_estimate = function(){
			private$fetch_and_or_cache_inference_model()$beta_hat_T
		},
		
		#' @description
		#' Computes a 1-alpha level frequentist confidence interval differently for all response types, estimate types and test types.
		#' 
		#' For "mean_difference" it computes
		#' (1a) for incidence outcomes (ignoring the KK design structure), 
		#' the p-value for the test of the additive log odds treatment effect being zero using logistic regression's MLE normal approximation
		#' (1b) for survival outcomes (ignoring the KK design structure), the median difference for survival using the Kaplan-Meier estimates for both arms 
		#' (1c) for count, proportion and continous outcomes (all ignoring the KK design structure), 
		#' the classic mean_difference estimate of the additive treatment effect, 
		#' (1d) for continuous outcome, equal allocation to arms and KK designs, there's a special match-reservoir weighted 
		#' classic mean_difference estimate
		#' 
		#' For "medial_difference" it computes only
		#' (2) for survival outcomes (ignoring the KK design structure), the difference of medians of the two arms
		#' 
		#' Computes for estimte type "default_regression" either
		#' (3a) for incidence outcomes, the additive log odds treatment effect using logistic regression controlled for all other covariates
		#' (3b) for survival outcomes, the additive treatment effect on log suvival using Weibull regression controlled for all other covariates
		#' (3c) for count outcomes, the additive treatment effect on log count using negative binomial regression controlled for all other covariates
		#' (3d) for proportion outcome, the additive treatment effect on proportion using beta regression controlled for all other covariates
		#' (3e) for continous outcomes but not under an equal allocation KK design, the additive treatment effect using OLS regression controlled for all other covariates
		#' (3f) for continuous outcome, equal allocation to arms and KK designs, there's a special match-reservoir weighted 
		#' OLS regression controlled for all other covariates
		#' 
		#' The confidence interval is computed differently for 
		#' [I] test type "MLE-or-KM-based"
		#' Here we use the theory that MLE's computed for GLM's are asymptotically normal (except in the case 
		#' of estimat_type "median difference" where a nonparametric bootstrap confidence interval (see the \code{controlTest::quantileControlTest} method)
		#' is employed. Hence these confidence intervals are asymptotically valid and thus approximate for any sample size.
		#' 
		#' [II] test type "randomization-exact"
		#' Here we invert the randomization test that tests the strong null H_0: y_T_i - y_C_i = delta <=> (y_T_i - delta) - y_C_i = 0 so 
		#' we adjust the treatment responses downward by delta. We then find the set of all delta values that is above 1 - alpha/2 (i.e. two-sided)
		#' This is accomplished via a bisection algorithm (algorithm 1 of Glazer and Stark, 2025 available at
		#' https://arxiv.org/abs/2405.05238). These confidence intervals are exact to within tolerance \code{pval_epsilon}.
		#' 
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#' @param nsim_exact_test		The number of randomization vectors (applicable for test type "randomization-exact" only). 
		#' 								The default is 1000 providing good resolutions to confidence intervals.
		#' @param B						Number of bootstrap samples for the survival response where \code{estimate_type} is "median_difference"
		#' 								(see the \code{controlTest::quantileControlTest} method). The default is NULL which corresponds to B=501 
		#' 								providing pvalue resolution to a fifth of a percent.
		#' @param pval_epsilon			The bisection algorithm tolerance for the test inversion (applicable for test type "randomization-exact" only). 
		#' 								The default is to find a CI accurate to within a tenth of a percent.
		#' 
		#' @return 	A 1 - alpha sized frequentist confidence interval for the treatment effect
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInference$new(seq_des, test_type = "MLE-or-KM-based")
		#' seq_des_inf$compute_confidence_interval()
		#' 		
		compute_confidence_interval = function(alpha = 0.05, nsim_exact_test = 501, pval_epsilon = 0.001, B = NULL){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			
			if (self$test_type == "MLE-or-KM-based"){
				private$handle_B_for_survival_median_inference(B)
				one_minus_alpha_over_two = 1 - alpha / 2
				inference_mod = private$fetch_and_or_cache_inference_model()
				z_or_t_val = 	if (inference_mod$is_z){
									qnorm(one_minus_alpha_over_two)
								} else {
									qt(one_minus_alpha_over_two, inference_mod$df)
								}
				moe = z_or_t_val * inference_mod$s_beta_hat_T
				#return the CI
				inference_mod$beta_hat_T + c(-moe, moe)
			} else { #randomization
				assertCount(nsim_exact_test, positive = TRUE)
				assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
				if (private$is_continuous){
					#to get bounds, temporarily look at the MLE CI
					#a CI with two orders of magnitude more in wideness should provide for us decent bounds to run the bisection search
					lower_upper_ci_bounds = private$compute_mle_or_km_based_confidence_interval(alpha / 100)
					stop("boom")
					c(
						private$compute_ci_lower_by_inverting_the_randomization_test(nsim_exact_test, private$num_cores, self$compute_treatment_estimate(), lower_upper_ci_bounds[2], alpha / 2, pval_epsilon),
						private$compute_ci_upper_by_inverting_the_randomization_test(nsim_exact_test, private$num_cores, lower_upper_ci_bounds[1], self$compute_treatment_estimate(), alpha / 2, pval_epsilon)
					)
				} else if (private$is_incidence){
					stop("Confidence intervals are not supported for randomization tests for mean difference in incidence outomes")
				} else if (private$is_count){
					stop("Confidence intervals are not supported for randomization tests for mean difference in count outomes")
				} else if (private$is_proportion){
					stop("Confidence intervals are not supported for randomization tests for mean difference in proportion outomes")
				} else if (private$is_survival){
					stop("Confidence intervals are not supported for randomization tests for mean difference in survival outomes")
				}
			}
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
		#' 								(see the \code{controlTest::quantileControlTest} method). The default is 501 providing pvalue resolution 
		#' 								to a fifth of a percent.
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#' 
		#' @return 	The approximate frequentist p-value
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInference$new(seq_des)
		#' seq_des_inf$compute_two_sided_pval_for_treatment_effect()
		#' 		
		compute_two_sided_pval_for_treatment_effect = function(nsim_exact_test = 501, B = NULL, delta = 0){
			assertNumeric(delta)
			private$handle_B_for_survival_median_inference(B)
			
			if (self$test_type == "MLE-or-KM-based"){
				if (delta != 0){
					stop("nonzero treatment effect tests not yet supported for MLE or KM based tests")
				}
				private$fetch_and_or_cache_inference_model()$p_val
				
			} else { #randomization
				assertCount(nsim_exact_test, positive = TRUE)
				private$compute_randomization_test_p_val(nsim_exact_test, private$num_cores, delta)				
			}
		}
	),
	
	private = list(
		seq_des_obj = NULL,
		isKK = FALSE,
		match_indic = NULL,
		match_indic_model_matrix = NULL,
		rand_inf_b_T_sims = NULL,
		num_cores = NULL,
		verbose = FALSE,
		n = NULL,
		p = NULL,
		X = NULL,
		yTs = NULL,
		yCs = NULL,
		deadTs = NULL,
		deadCs = NULL,
		beta_T_hat = NULL,
		prob_T = NULL,
		prob_T_fifty_fifty = NULL,
		KKstats = NULL,
		
		inference_model = NULL,
		default_B_for_median_inference = 501,
		is_continuous = FALSE,
		is_incidence = FALSE,
		is_count = FALSE,
		is_proportion = FALSE,
		is_survival = FALSE,
		
		handle_B_for_survival_median_inference = function(B){
			if (is.null(B)){
				B = private$default_B_for_median_inference
			}
			assertCount(B, positive = TRUE)
			if (self$estimate_type == "survival_simple_median_difference" & B != private$default_B_for_median_inference){
				private$inference_model = private$compute_survival_simple_median_inference(B) #recompute it since we now need to pass in the new B value
			}
		},
		
		############# RANDOMIZATION TEST HELPER FUNCTIONS
		
		compute_ci_lower_by_inverting_the_randomization_test = function(nsim_exact_test, num_cores, l, u, pval_th, tol){
			return(NA)
			pval_l = private$compute_randomization_test_p_val(nsim_exact_test, num_cores, l)
			pval_u = private$compute_randomization_test_p_val(nsim_exact_test, num_cores, u)
			repeat {
				if (pval_u - pval_l <= tol){
					return(l)
				}
				m = (l + u) / 2
				pval_m = private$compute_randomization_test_p_val(nsim_exact_test, num_cores, m)
				if (pval_m >= pval_th){
					u = m
					pval_u = pval_m
				} else {
					l = m
					pval_l = pval_m
				}
			}
		},
		compute_ci_upper_by_inverting_the_randomization_test = function(nsim_exact_test, num_cores, l, u, pval_th, tol){	
			return(NA)
			pval_l = private$compute_randomization_test_p_val(nsim_exact_test, num_cores, l)
			pval_u = private$compute_randomization_test_p_val(nsim_exact_test, num_cores, u)
			repeat {
				if (pval_l - pval_u <= tol){
					return(u)
				}
				m = (l + u) / 2
				pval_m = private$compute_randomization_test_p_val(nsim_exact_test, num_cores, m)
				if (pval_m >= pval_th){
					l = m
					pval_l = pval_m
				} else {
					u = m
					pval_u = pval_m
				}
			}
		},
		
		compute_randomization_test_p_val = function(nsim_exact_test, num_cores, delta){			
			y = private$seq_des_obj$y
			if (delta != 0){
				if (!private$is_continuous){
					stop("randomization tests with delta nonzero only works for continuous type!!!!")
				}
				#we are testing against H_0: y_T_i - y_C_i = delta <=> (y_T_i - delta) - y_C_i = 0 so adjust the treatment responses downward by delta
				y[private$seq_des_obj$w == 1] = y[private$seq_des_obj$w == 1] - delta
			}
			
			estimate_type = self$estimate_type
			seq_des_obj = private$seq_des_obj
			if (private$num_cores == 1){
				#easier on the OS I think...
				b_T_sims = array(NA, nsim_exact_test)
				for (r in 1 : nsim_exact_test){
					#cat("		r =", r, "/", nsim_exact_test, "\n")
					#make a copy of the object and then permute the allocation vector according to the design
					seq_des_r = seq_des_obj$.__enclos_env__$private$duplicate()
					seq_des_r$y = y #set the new responses
					seq_des_r$.__enclos_env__$private$redraw_w_according_to_design()
					b_T_sims[r] = SeqDesignInference$new(seq_des_r, estimate_type = estimate_type, verbose = FALSE)$compute_treatment_estimate()
				}
				#print(ggplot2::ggplot(data.frame(sims = b_T_sims)) + ggplot2::geom_histogram(ggplot2::aes(x = sims), bins = 50))
			} else {	
				cl = makeCluster(private$num_cores)
				registerDoParallel(cl)			
				#now copy them to each core's memory
				clusterExport(cl, list("seq_des_obj", "y", "estimate_type"), envir = environment())
				#now do the parallelization
				b_T_sims = foreach(r = 1 : nsim_exact_test, .inorder = FALSE, .combine = c) %dopar% {
					#make a copy of the object and then permute the allocation vector according to the design
					seq_des_r = seq_des_obj$.__enclos_env__$private$duplicate()
					seq_des_r$y = y #set the new responses
					seq_des_r$.__enclos_env__$private$redraw_w_according_to_design()				
					SeqDesignInference$new(seq_des_r, estimate_type = estimate_type, verbose = FALSE)$compute_treatment_estimate()				
				}
				stopCluster(cl)
				rm(cl); gc()
			}
			#this calculates the two-sided pval
			beta_hat_T = self$compute_treatment_estimate()
			#sum(abs(beta_hat_T ) < abs(private$rand_inf_b_T_sims)) / nsim_exact_test ####check this
			2 * min(
				sum(b_T_sims > beta_hat_T) / nsim_exact_test, 
				sum(b_T_sims < beta_hat_T) / nsim_exact_test
			)
		},
		
		######### MODEL INFERENCE FUNCTIONS	
		
		fetch_and_or_cache_inference_model = function(){
			if (is.null(private$inference_model)){
				private$inference_model = switch(self$estimate_type,
					########################################### ALL RESPONSE TYPES
					"simple_mean_difference" =
							private$compute_simple_mean_difference_inference(),
					"KK_compound_mean_difference" =
							private$compute_KK_compound_mean_difference_inference(),  
					########################################### CONTINUOUS
					"continuous_regression_with_covariates" =
							private$compute_continuous_multivariate_ols_inference(),	
					"continuous_KK_compound_multivariate_regression" =
							private$compute_continuous_KK_compound_multivariate_ols_inference(),
					"continuous_KK_regression_with_covariates_with_matching_dummies" = 
							private$compute_continuous_KK_multivariate_with_matching_dummies_ols_inference(),
					"continuous_KK_regression_with_covariates_with_random_intercepts" = 
							private$compute_continuous_KK_multivariate_and_matching_random_intercepts_regression_inference(),
					########################################### INCIDENCE
					"incidence_simple_log_odds" =
							private$compute_incidence_univariate_logistic_regression_inference(),	
					"incidence_logistic_regression" =
							private$compute_incidence_multivariate_logistic_regression_inference(),
					"incidence_KK_compound_univariate_logistic_regression" =
							private$compute_incidence_KK_compound_univariate_logistic_regression_inference(),
					"incidence_KK_compound_multivariate_logistic_regression" =
							private$compute_incidence_KK_compound_multivariate_logistic_regression_inference(),	
					"incidence_KK_multivariate_logistic_regression_with_matching_dummies" =
							private$compute_incidence_KK_multivariate_logistic_regression_with_matching_dummies_inference(),	
					"incidence_KK_multivariate_logistic_regression_with_random_intercepts_for_matches" =
							private$compute_incidence_KK_multivariate_logistic_regression_with_random_intercepts_for_matches_inference(),
					########################################### PROPORTION
					"proportion_simple_logodds_regression" =
							private$compute_proportion_univariate_beta_regression_inference(),
					"proportion_beta_regression" =
							private$compute_proportion_multivariate_beta_regression_inference(),
					"proportion_KK_compound_univariate_beta_regression" =
							private$compute_proportion_KK_compound_univariate_beta_regression_inference(),
					"proportion_KK_compound_multivariate_beta_regression" =
							private$compute_proportion_KK_compound_multivariate_beta_regression_inference(),
					"proportion_KK_multivariate_beta_regression_with_matching_dummies" =
							private$compute_proportion_KK_multivariate_beta_regression_with_matching_dummies_inference(),
					########################################### COUNT
					"count_univariate_negative_binomial_regression" =
							private$compute_count_univariate_negative_binomial_inference(),
					"count_multivariate_negative_binomial_regression" =
							private$compute_count_multivariate_negative_binomial_inference(),
					"count_KK_compound_univariate_negative_binomial_regression" =
							private$compute_count_KK_compound_univariate_negative_binomial_inference(),	
					"count_KK_compound_multivariate_negative_binomial_regression" =
							private$compute_count_KK_compound_multivariate_negative_binomial_inference(),
					"count_KK_multivariate_negative_binomial_regression_with_matching_dummies" =
							private$compute_count_KK_multivariate_with_matching_dummies_negative_binomial_inference(),
					"count_KK_multivariate_negative_binomial_regression_with_random_intercepts_for_matches" =
							private$compute_count_KK_multivariate_negative_binomial_with_random_intercepts_for_matches_inference(),
					########################################### SURVIVAL
					"survival_simple_median_difference" =
							private$compute_survival_simple_median_inference(),	
					"survival_simple_restricted_mean_difference" = 
							private$compute_survival_simple_restricted_mean_inference(),
					"survival_univariate_weibull_regression" =
							private$compute_survival_univariate_weibull_regression_inference(),	
					"survival_multivariate_weibull_regression" =
							private$compute_survival_multivariate_weibull_regression_inference(),
					"survival_KK_compound_univariate_weibull_regression" =
							private$compute_survival_KK_compound_univariate_weibull_inference(),	
					"survival_KK_compound_multivariate_weibull_regression" =
							private$compute_survival_KK_compound_multivariate_weibull_inference(),
					"survival_KK_multivariate_weibull_regression_with_matching_dummies" =
							private$compute_survival_multivariate_with_matching_dummies_weibull_regression_inference(),					
					"survival_univariate_coxph_regression" =
							private$compute_survival_univariate_coxph_regression_inference(),	
					"survival_multivariate_coxph_regression" =
							private$compute_survival_multivariate_coxph_regression_inference(),		
					"survival_KK_multivariate_coxph_regression_with_matching_dummies" =
							private$compute_survival_multivariate_with_matching_dummies_coxph_regression_inference(),		
					"survival_KK_multivariate_coxph_regression_with_random_intercepts_for_matches" =
							private$compute_survival_multivariate_with_random_intercepts_coxph_regression_inference()	
				)
			}
			private$inference_model
		},
		
		compute_simple_mean_difference_inference = function(){
			nT = length(private$yTs)
			nC = length(private$yCs)
			beta_hat_T = mean(private$yTs) - mean(private$yCs)
			s_1_sq = var(private$yTs) / nT 
			s_2_sq = var(private$yCs) / nC
			s_beta_hat_T = sqrt(s_1_sq + s_2_sq)
			df = (s_1_sq + s_2_sq)^2 / 
					(s_1_sq^2 / (nT - 1) + s_2_sq^2 / (nC - 1)) #Welch-Satterthwaite formula
			list(
				beta_hat_T = beta_hat_T,
				s_beta_hat_T = s_beta_hat_T, 
				is_z = FALSE,
				df = df,
				p_val = 2 * pt(-abs(beta_hat_T) / s_beta_hat_T, df)
			)		
		},	
		
		compute_KK_compound_mean_difference_inference = function(){
			KKstats = private$compute_continuous_post_matching_data_KK()
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
			KKstats$p_val = 2 * (pnorm(-abs(KKstats$beta_hat_T / KKstats$s_beta_hat_T))) #approximate by using N(0, 1) distribution			
			KKstats$is_z = TRUE #see KK14 paper for details about how assuming the normal distr here may be too liberal
			KKstats
		},
		
		compute_continuous_multivariate_ols_inference = function(){
			tryCatch({
				ols_regr_mod = lm(private$seq_des_obj$y ~ ., data = cbind(data.frame(w = private$seq_des_obj$w), private$X))
				summary_table = coef(summary(ols_regr_mod))
				list(
					mod = ols_regr_mod,
					summary_table = summary_table,	
					beta_hat_T = summary_table[2, 1],
					s_beta_hat_T = summary_table[2, 2],
					is_z = FALSE,
					df = private$n - nrow(summary_table),
					p_val = summary_table[2, 4]
				)
			}, error = function(e){ #very difficult to get rid of errors here due to Error in vcov.merMod(object, use.hessian = use.hessian)... tried to write a robust function but it didn't work
				list(
					mod = NA,
					summary_table = NA,	
					beta_hat_T = NA,
					s_beta_hat_T = NA,
					is_z = TRUE,
					p_val = NA
				)				
			})				
		},
		
		compute_continuous_KK_compound_multivariate_ols_inference = function(){
			KKstats = private$compute_continuous_post_matching_data_KK()
			
			if (KKstats$nRT <= 2 || KKstats$nRC <= 2 || (KKstats$nRT + KKstats$nRC <= ncol(private$X) + 2)){
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
			KKstats$is_z = TRUE #see KK14 paper for details about how assuming the normal distr here may be too liberal
			KKstats
		},
		
		compute_continuous_KK_multivariate_with_matching_dummies_ols_inference = function(){
			tryCatch({
				ols_regr_mod = lm(private$seq_des_obj$y ~ ., 
						data = private$generate_data_frame_with_matching_dummies())
				summary_table = suppressWarnings(coef(summary(ols_regr_mod)))
				list(
					mod = ols_regr_mod,
					summary_table = summary_table,	
					beta_hat_T = summary_table[2, 1],
					s_beta_hat_T = summary_table[2, 2],
					is_z = FALSE,
					df = private$n - nrow(summary_table),
					p_val = summary_table[2, 4]
				)
			}, error = function(e){ #very difficult to get rid of errors here due to Error in vcov.merMod(object, use.hessian = use.hessian)... tried to write a robust function but it didn't work
				list(
					mod = NA,
					summary_table = NA,	
					beta_hat_T = NA,
					s_beta_hat_T = NA,
					is_z = TRUE,
					p_val = NA
				)				
			})				
		},
		
		generate_data_frame_with_matching_dummies = function(){
			cbind(data.frame(w = private$seq_des_obj$w), private$X, private$match_indic_model_matrix)
		},
		
		compute_continuous_KK_multivariate_and_matching_random_intercepts_regression_inference = function(){
			tryCatch({
				mixed_regr_mod = suppressWarnings(lmerTest::lmer(y ~ . - match_indic + (1 | match_indic), 
						data = cbind(data.frame(y = private$seq_des_obj$y, w = private$seq_des_obj$w, match_indic = factor(private$match_indic)), private$X)))
				summary_table = coef(summary(mixed_regr_mod))
				list(
					mod = mixed_regr_mod,
					summary_table = summary_table,	
					beta_hat_T = summary_table[2, 1],
					s_beta_hat_T = summary_table[2, 2],
					is_z = FALSE,
					df = summary_table[2, 3],
					p_val = summary_table[2, 5]
				)
			}, error = function(e){ #very difficult to get rid of errors here due to Error in vcov.merMod(object, use.hessian = use.hessian)... tried to write a robust function but it didn't work
				list(
					mod = NA,
					summary_table = NA,	
					beta_hat_T = NA,
					s_beta_hat_T = NA,
					is_z = TRUE,
					p_val = NA
				)				
			})			
		},
		
		compute_incidence_univariate_logistic_regression_inference = function(){
			tryCatch({
				logistic_regr_mod = suppressWarnings(glm(private$seq_des_obj$y ~ ., 
						data = data.frame(w = private$seq_des_obj$w), family = "binomial"))
				summary_table = coef(summary_glm_lean(logistic_regr_mod))
				list(
					mod = logistic_regr_mod,
					summary_table = summary_table,	
					beta_hat_T = summary_table[2, 1],
					s_beta_hat_T = summary_table[2, 2],
					is_z = TRUE,
					p_val = summary_table[2, 4]
				)
			}, error = function(e){ #very difficult to get rid of errors here due to Error in vcov.merMod(object, use.hessian = use.hessian)... tried to write a robust function but it didn't work
				list(
					mod = NA,
					summary_table = NA,	
					beta_hat_T = NA,
					s_beta_hat_T = NA,
					is_z = TRUE,
					p_val = NA
				)				
			})			
		},
		
		compute_incidence_multivariate_logistic_regression_inference = function(){
			tryCatch({
				logistic_regr_mod = suppressWarnings(glm(private$seq_des_obj$y ~ ., 
						data = cbind(data.frame(w = private$seq_des_obj$w), private$X), family = "binomial"))
				summary_table = coef(summary_glm_lean(logistic_regr_mod))
				list(
					mod = logistic_regr_mod,
					summary_table = summary_table,	
					beta_hat_T = summary_table[2, 1],
					s_beta_hat_T = summary_table[2, 2],
					is_z = TRUE,
					p_val = summary_table[2, 4]
				)
			}, error = function(e){ #very difficult to get rid of errors here due to Error in vcov.merMod(object, use.hessian = use.hessian)... tried to write a robust function but it didn't work
				list(
					mod = NA,
					summary_table = NA,	
					beta_hat_T = NA,
					s_beta_hat_T = NA,
					is_z = TRUE,
					p_val = NA
				)				
			})
		},
		
		compute_incidence_KK_compound_univariate_logistic_regression_inference = function(){
			stop("not implemented yet")
		},
		
		compute_incidence_KK_compound_multivariate_logistic_regression_inference = function(){
			stop("not implemented yet")
		},
		
		compute_incidence_KK_multivariate_logistic_regression_with_matching_dummies_inference = function(){
			tryCatch({
				logistic_regr_mod = suppressWarnings(glm(private$seq_des_obj$y ~ ., 
						data = private$generate_data_frame_with_matching_dummies(), family = "binomial"))
				summary_table = coef(summary_glm_lean(logistic_regr_mod))
				list(
					mod = logistic_regr_mod,
					summary_table = summary_table,	
					beta_hat_T = summary_table[2, 1],
					s_beta_hat_T = summary_table[2, 2],
					is_z = TRUE,
					p_val = summary_table[2, 4]
				)						
			}, error = function(e){ #very difficult to get rid of errors here due to Error in vcov.merMod(object, use.hessian = use.hessian)... tried to write a robust function but it didn't work
				list(
					mod = NA,
					summary_table = NA,	
					beta_hat_T = NA,
					s_beta_hat_T = NA,
					is_z = TRUE,
					p_val = NA
				)				
			})
		},
		
		compute_incidence_KK_multivariate_logistic_regression_with_random_intercepts_for_matches_inference = function(){
			tryCatch({
				mixed_logistic_regr_mod = suppressWarnings(lme4::glmer(y ~ . - match_indic + (1 | match_indic), 
						data = cbind(data.frame(y = private$seq_des_obj$y, w = private$seq_des_obj$w, match_indic = factor(private$match_indic)), private$X),
						family = "binomial"))
				summary_table = coef(summary(mixed_logistic_regr_mod))
				list(
					mod = mixed_logistic_regr_mod,
					summary_table = summary_table,	
					beta_hat_T = summary_table[2, 1],
					s_beta_hat_T = summary_table[2, 2],
					is_z = TRUE,
					p_val = summary_table[2, 4]
				)
			}, error = function(e){ #very difficult to get rid of errors here due to Error in vcov.merMod(object, use.hessian = use.hessian)... tried to write a robust function but it didn't work
				list(
					mod = NA,
					summary_table = NA,	
					beta_hat_T = NA,
					s_beta_hat_T = NA,
					is_z = TRUE,
					p_val = NA
				)				
			})
		},
		
		compute_proportion_univariate_beta_regression_inference = function(){
			private$shared_beta_regression_inference(data.frame(y = private$seq_des_obj$y, w = private$seq_des_obj$w))
		},
		
		compute_proportion_multivariate_beta_regression_inference = function(){
			private$shared_beta_regression_inference(cbind(data.frame(y = private$seq_des_obj$y, w = private$seq_des_obj$w), private$X))
		},
		
		compute_proportion_KK_multivariate_beta_regression_with_matching_dummies_inference = function(){
			private$shared_beta_regression_inference(cbind(data.frame(y = private$seq_des_obj$y), private$generate_data_frame_with_matching_dummies()))
		},
		
		compute_proportion_KK_compound_univariate_beta_regression_inference = function(){
			stop("not implemented yet")
		}, 
		
		compute_proportion_KK_compound_multivariate_beta_regression_inference = function(){
			stop("not implemented yet")
		},
		
		shared_beta_regression_inference = function(data_obj){
			tryCatch({
				beta_regr_mod = robust_betareg(y ~ ., data_obj)
				beat_regr_mod_sub = coef(summary(beta_regr_mod))
				summary_table = if (!is.null(beat_regr_mod_sub$mean)){ #beta model
							beat_regr_mod_sub$mean 
						} else if (!is.null(beat_regr_mod_sub$mu)){ #extended-support xbetax model
							beat_regr_mod_sub$mu
						}
				list(
					mod = beta_regr_mod,
					summary_table = summary_table,	
					beta_hat_T = summary_table[2, 1],
					s_beta_hat_T = summary_table[2, 2],
					is_z = TRUE,
					p_val = summary_table[2, 4]
				)			
			}, error = function(e){
				list(
					mod = NA,
					summary_table = NA,	
					beta_hat_T = NA,
					s_beta_hat_T = NA,
					is_z = TRUE,
					p_val = NA
				)				
			})
		},
		
		compute_count_univariate_negative_binomial_inference = function(){
			private$compute_count_negative_binomial_regression(data.frame(y = private$seq_des_obj$y, w = private$seq_des_obj$w))
		},
		
		compute_count_multivariate_negative_binomial_inference = function(){
			private$compute_count_negative_binomial_regression(cbind(data.frame(y = private$seq_des_obj$y, w = private$seq_des_obj$w), private$X))
		},
		
		compute_count_KK_compound_univariate_negative_binomial_inference = function(){
			stop("not implemented yet")
		},
		
		compute_count_KK_compound_multivariate_negative_binomial_inference = function(){
			stop("not implemented yet")
		},
		
		compute_count_KK_multivariate_negative_binomial_with_random_intercepts_for_matches_inference = function(){			
			tryCatch({
				mixed_neg_bin_regr_mod = suppressWarnings(lme4::glmer.nb(y ~ . - match_indic + (1 | match_indic), 
								data = cbind(data.frame(y = private$seq_des_obj$y, w = private$seq_des_obj$w, match_indic = factor(private$match_indic)), private$X)))
				summary_table = coef(summary(mixed_neg_bin_regr_mod))
				list(
					mod = mixed_neg_bin_regr_mod,
					summary_table = summary_table,	
					beta_hat_T = summary_table[2, 1],
					s_beta_hat_T = summary_table[2, 2],
					is_z = TRUE,
					p_val = summary_table[2, 4]
				)						
			}, error = function(e){ #very difficult to get rid of errors here due to VTV not positive definite... tried to write a robust function but it didn't work
				list(
					mod = NA,
					summary_table = NA,	
					beta_hat_T = NA,
					s_beta_hat_T = NA,
					is_z = TRUE,
					p_val = NA
				)				
			})
		},
		
		compute_count_KK_multivariate_with_matching_dummies_negative_binomial_inference = function(){
			private$compute_count_negative_binomial_regression(cbind(data.frame(y = private$seq_des_obj$y), private$generate_data_frame_with_matching_dummies()))
		},
		
		compute_count_negative_binomial_regression = function(data_obj){
			tryCatch({
				negbin_regr_mod = suppressWarnings(MASS::glm.nb(y ~ ., data = data_obj))
				summary_table = coef(summary_glm_lean(negbin_regr_mod))
				list(
					mod = negbin_regr_mod,
					summary_table = summary_table,	
					beta_hat_T = summary_table[2, 1],
					s_beta_hat_T = summary_table[2, 2],
					is_z = TRUE,
					p_val = summary_table[2, 4]
				)			
			}, error = function(e){ #very difficult to get rid of errors here due to VTV not positive definite... tried to write a robust function but it didn't work
				list(
					mod = NA,
					summary_table = NA,	
					beta_hat_T = NA,
					s_beta_hat_T = NA,
					is_z = TRUE,
					p_val = NA
				)				
			})
		},
		
		compute_survival_simple_median_inference = function(B = NULL){
			if (is.null(B)){
				B = private$default_B_for_median_inference
			}
			survival_obj = survival::Surv(private$seq_des_obj$y, private$seq_des_obj$dead)
			survival_fit_obj = survival::survfit(survival_obj ~ private$seq_des_obj$w)
			survival_fit_res = summary(survival_fit_obj)$table
			beta_hat_T = survival_fit_res[2, 7] - survival_fit_res[1, 7]
			if (is.na(beta_hat_T)){
				list(
					mod = survival_fit_obj,
					summary_table = survival_fit_res,
					beta_hat_T = beta_hat_T,
					s_beta_hat_T = NA,
					is_z = TRUE,
					p_val = NA
				)
			} else {
				test_obj = suppressWarnings(controlTest::quantileControlTest(private$yTs, private$deadTs, private$yCs, private$deadCs, B = B))
				list(
					mod = survival_fit_obj,
					summary_table = survival_fit_res,
					beta_hat_T = beta_hat_T,
					s_beta_hat_T = abs(beta_hat_T / test_obj$Z), #wtf is quantileControlTest's se field if not this??????????
					is_z = TRUE,
					p_val = test_obj$pval
				)
			}
		},
		
		compute_survival_simple_restricted_mean_inference = function(){
			survival_obj = survival::Surv(private$seq_des_obj$y, private$seq_des_obj$dead)
			survival_fit_obj = survival::survfit(survival_obj ~ private$seq_des_obj$w)
			survival_fit_res = summary(survival_fit_obj)$table
			beta_hat_T = survival_fit_res[2, 5] - survival_fit_res[1, 5]
			s_beta_hat_T = sqrt(survival_fit_res[2, 6]^2 + survival_fit_res[1, 6]^2)
			z_beta_hat_T = beta_hat_T / s_beta_hat_T
			p_val = 2 * min(pnorm(z_beta_hat_T), 1 - pnorm(z_beta_hat_T))
			list(
				mod = survival_fit_obj,
				summary_table = survival_fit_res,
				beta_hat_T = beta_hat_T,
				s_beta_hat_T = s_beta_hat_T,
				is_z = TRUE,
				p_val = p_val
			)
		},
				
		compute_survival_univariate_weibull_regression_inference = function(){
			private$compute_survival_weibull_regression(private$seq_des_obj$w)
		},
		
		compute_survival_multivariate_weibull_regression_inference = function(){
			private$compute_survival_weibull_regression(cbind(private$seq_des_obj$w, private$X))
		},
		
		compute_survival_KK_compound_univariate_weibull_inference = function(){
			stop("not implemented yet")
		},
		
		compute_survival_KK_compound_multivariate_weibull_inference = function(){
			stop("not implemented yet")
		},
		
		compute_survival_multivariate_with_matching_dummies_weibull_regression_inference = function(){
			private$compute_survival_weibull_regression(private$generate_data_frame_with_matching_dummies())
		},
		
		compute_survival_weibull_regression = function(data_obj){
			surv_regr_mod = robust_survreg(private$seq_des_obj$y, private$seq_des_obj$dead, data_obj)
			if (is.null(surv_regr_mod)){
				list(
					mod = surv_regr_mod,
					summary_table = NA,	
					beta_hat_T = NA,
					s_beta_hat_T = NA,
					is_z = TRUE,
					p_val = NA
				)
			} else {
				summary_table = summary(surv_regr_mod)$table
				if (is.na(summary_table[2, 1])){
					list(
						mod = surv_regr_mod,
						summary_table = NA,	
						beta_hat_T = NA,
						s_beta_hat_T = NA,
						is_z = TRUE,
						p_val = NA
					)
				} else {
					list(
						mod = surv_regr_mod,
						summary_table = summary_table,	
						beta_hat_T = summary_table[2, 1],
						s_beta_hat_T = summary_table[2, 2],
						is_z = TRUE,
						p_val = summary_table[2, 4]
					)
				}				
			}

		},
		
		compute_survival_univariate_coxph_regression_inference = function(){
			private$compute_cox_regression(data.frame(w = private$seq_des_obj$w))
		},
		
		compute_survival_multivariate_coxph_regression_inference = function(){
			private$compute_cox_regression(cbind(data.frame(w = private$seq_des_obj$w), private$X))
		},
		
		compute_survival_multivariate_with_matching_dummies_coxph_regression_inference = function(){
			private$compute_cox_regression(private$generate_data_frame_with_matching_dummies())
		},
		
		compute_cox_regression = function(data_obj){
			y = private$seq_des_obj$y
			dead = private$seq_des_obj$dead
			
			stop("boom")
			tryCatch({
				surv_obj = survival::Surv(y, dead)
				coxph_mod = suppressWarnings(coxph(surv_obj ~ ., data = data_obj))
				summary_table = coef(summary(coxph_mod))
				list(
					mod = coxph_mod,
					summary_table = summary_table,	
					beta_hat_T = summary_table[1, 1],
					s_beta_hat_T = summary_table[1, 3],
					is_z = TRUE,
					p_val = summary_table[1, 5]
				)
			}, error = function(e){ #very difficult to get rid of errors here due to VTV not positive definite... tried to write a robust function but it didn't work
				list(
					mod = NA,
					summary_table = NA,	
					beta_hat_T = NA,
					s_beta_hat_T = NA,
					is_z = TRUE,
					p_val = NA
				)				
			})			
		},
		
		compute_survival_multivariate_with_random_intercepts_coxph_regression_inference = function(){
			y = private$seq_des_obj$y
			dead = private$seq_des_obj$dead
			surv_obj = survival::Surv(y, dead)
			X = private$X
			colnames(X) = paste0("x", 1 : ncol(X)) #as there may be spaces in the colnames which blows the formula up... so damn annoying
			data = cbind(data.frame(w = private$seq_des_obj$w, match_indic = factor(private$match_indic)), X)
			f = as.formula(paste("surv_obj ~ (1 | match_indic) + w +", paste(colnames(X), collapse = "+")))
			tryCatch({
				coxph_mod = suppressWarnings(coxme::coxme(f, data = data)) #https://stats.oarc.ucla.edu/r/dae/mixed-effects-cox-regression/
				summary_table = coef(summary(coxph_mod))
				list(
					mod = coxph_mod,
					summary_table = summary_table,	
					beta_hat_T = summary_table[1, 2],
					s_beta_hat_T = summary_table[1, 3],
					is_z = TRUE,
					p_val = summary_table[1, 5]
				)					
			}, error = function(e){
				list(
					mod = NA,
					summary_table = NA,	
					beta_hat_T = NA,
					s_beta_hat_T = NA,
					is_z = TRUE,
					p_val = NA
				)
			})
		},
				
		compute_continuous_post_matching_data_KK = function(){	
			if (is.null(private$KKstats)){
				#get matched data
				
				m = max(private$match_indic)
				y_matched_diffs = array(NA, m)
				X_matched_diffs = matrix(NA, nrow = m, ncol = ncol(private$X))
				if (m > 0){
					for (match_id in 1 : m){ #we want to just calculate the diffs inside matches and ignore the reservoir
						yTs = private$seq_des_obj$y[private$seq_des_obj$w == 1 & private$match_indic == match_id]
						yCs = private$seq_des_obj$y[private$seq_des_obj$w == 0 & private$match_indic == match_id]
						y_matched_diffs[match_id] = mean(yTs) - mean(yCs)
						
						xmTs = private$X[private$seq_des_obj$w == 1 & private$match_indic == match_id, ]
						xmCs = private$X[private$seq_des_obj$w == 0 & private$match_indic == match_id, ]
						X_matched_diffs[match_id, ] = mean(xmTs) - mean(xmCs)
					}
				}			
				
				#get reservoir data
				X_reservoir = 	private$X[private$match_indic == 0, ]
				y_reservoir = 	private$seq_des_obj$y[private$match_indic == 0]
				w_reservoir = 	private$seq_des_obj$w[private$match_indic == 0]
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