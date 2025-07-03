		#' @param estimate_type		The type of estimate to compute of which there are many and identified by the response type as its first word (except
		#' 							\code{simple_mean_difference} and \code{KK_compound_mean_difference} which are available for all response types except
		#' 							survival when there are any censored observations). If the string "KK" 
		#' 							appears after the first word, then this estimate type is only applicable to KK14, KK21, KK21stepwise designs.
		#' 							  * "simple_mean_difference" (or "*_simple_mean_difference" where * is any of the five supported response types)
		#' 								assumes the treatment effect parameter is an additive treatment effect 
		#' 							  * "KK_compound_mean_difference" (or "*_KK_compound_mean_difference" where * is any of the five supported response types) 	
		#' 								assumes the treatment effect parameter is an additive treatment effect	
		#' 								and estimates via combining a simple average difference estimator for both the matches and the reservoir
		#' 								and estimates via the simple average difference
		#' 							  * "continuous_multivariate_regression"
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
		#' 							  * "incidence_multivariate_logistic_regression"
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
		#' 							  * "proportion_multivariate_beta_regression"
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
		
		
		
		
					assertChoice(estimate_type, c(
				########################################### ALL RESPONSE TYPES
				"continuous_simple_mean_difference", "incidence_simple_mean_difference", "proportion_simple_mean_difference", "count_simple_mean_difference", "survival_simple_mean_difference",
				"continuous_KK_compound_mean_difference", "incidence_KK_compound_mean_difference", "proportion_KK_compound_mean_difference", "count_KK_compound_mean_difference", "survival_KK_compound_mean_difference",
				########################################### CONTINUOUS
				"continuous_multivariate_regression",
				"continuous_KK_compound_multivariate_regression",
				"continuous_KK_regression_with_covariates_with_matching_dummies",
				"continuous_KK_regression_with_covariates_with_random_intercepts",
				########################################### INCIDENCE
				"incidence_simple_log_odds",	
				"incidence_multivariate_logistic_regression",
				"incidence_KK_compound_univariate_logistic_regression",
				"incidence_KK_compound_multivariate_logistic_regression",	
				"incidence_KK_multivariate_logistic_regression_with_matching_dummies",	
				"incidence_KK_multivariate_logistic_regression_with_random_intercepts_for_matches",
				########################################### PROPORTION
				"proportion_simple_logodds_regression",
				"proportion_multivariate_beta_regression",
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
	
	
	
	
		compute_continuous_KK_multivariate_with_matching_dummies_ols_inference = function(){
			tryCatch({
				ols_regr_mod = lm(private$seq_des_obj_priv_int$y ~ ., 
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
		
		
		
		
		
		
				compute_continuous_KK_multivariate_and_matching_random_intercepts_regression_inference = function(){
			tryCatch({
				mixed_regr_mod = suppressWarnings(lmerTest::lmer(y ~ . - match_indic + (1 | match_indic), 
						data = cbind(data.frame(y = private$seq_des_obj_priv_int$y, w = private$seq_des_obj_priv_int$w, match_indic = factor(private$match_indic)), private$get_X())))
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