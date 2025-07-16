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
		
		
		
			fetch_and_or_cache_inference_model = function(compute_ci = FALSE, compute_pval = FALSE){
			if (is.null(private$inference_model)){
				private$inference_model = switch(self$estimate_type,
					########################################### APPLICABLE TO ALL RESPONSE TYPES
#					"continuous_simple_mean_difference" =, "incidence_simple_mean_difference" =, "proportion_simple_mean_difference" =, "count_simple_mean_difference" =, "survival_simple_mean_difference" =
#							private$compute_simple_mean_difference_inference(compute_ci, compute_pval),
#					"continuous_KK_compound_mean_difference" =, "incidence_KK_compound_mean_difference" =, "proportion_KK_compound_mean_difference" =, "count_KK_compound_mean_difference" =, "survival_KK_compound_mean_difference" =
#							private$compute_KK_compound_mean_difference_inference(compute_ci, compute_pval),  
					########################################### CONTINUOUS
#					"continuous_multivariate_regression" =
#							private$compute_continuous_multivariate_ols_inference(compute_ci, compute_pval),	
#					"continuous_KK_compound_multivariate_regression" =
#							private$compute_continuous_KK_compound_multivariate_ols_inference(compute_ci, compute_pval),
					"continuous_KK_regression_with_covariates_with_matching_dummies" = 
							private$compute_continuous_KK_multivariate_with_matching_dummies_ols_inference(compute_ci, compute_pval),
					"continuous_KK_regression_with_covariates_with_random_intercepts" = 
							private$compute_continuous_KK_multivariate_and_matching_random_intercepts_regression_inference(compute_ci, compute_pval),
					########################################### INCIDENCE
#					"incidence_simple_log_odds" =
#							private$compute_incidence_univariate_logistic_regression_inference(compute_ci, compute_pval),	
#					"incidence_multivariate_logistic_regression" =
#							private$compute_incidence_multivariate_logistic_regression_inference(compute_ci, compute_pval),
#					"incidence_KK_compound_univariate_logistic_regression" =
#							private$compute_incidence_KK_compound_univariate_logistic_regression_inference(compute_ci, compute_pval),
					"incidence_KK_compound_multivariate_logistic_regression" =
							private$compute_incidence_KK_compound_multivariate_logistic_regression_inference(compute_ci, compute_pval),	
					"incidence_KK_multivariate_logistic_regression_with_matching_dummies" =
							private$compute_incidence_KK_multivariate_logistic_regression_with_matching_dummies_inference(compute_ci, compute_pval),	
					"incidence_KK_multivariate_logistic_regression_with_random_intercepts_for_matches" =
							private$compute_incidence_KK_multivariate_logistic_regression_with_random_intercepts_for_matches_inference(compute_ci, compute_pval),
					########################################### PROPORTION
#					"proportion_simple_logodds_regression" =
#							private$compute_proportion_univariate_beta_regression_inference(compute_ci, compute_pval),
#					"proportion_multivariate_beta_regression" =
#							private$compute_proportion_multivariate_beta_regression_inference(compute_ci, compute_pval),
#					"proportion_KK_compound_univariate_beta_regression" =
							private$compute_proportion_KK_compound_univariate_beta_regression_inference(compute_ci, compute_pval),
					"proportion_KK_compound_multivariate_beta_regression" =
							private$compute_proportion_KK_compound_multivariate_beta_regression_inference(compute_ci, compute_pval),
					"proportion_KK_multivariate_beta_regression_with_matching_dummies" =
							private$compute_proportion_KK_multivariate_beta_regression_with_matching_dummies_inference(compute_ci, compute_pval),
					########################################### COUNT
#					"count_univariate_negative_binomial_regression" =
#							private$compute_count_univariate_negative_binomial_inference(compute_ci, compute_pval),
#					"count_multivariate_negative_binomial_regression" =
#							private$compute_count_multivariate_negative_binomial_inference(compute_ci, compute_pval),
#					"count_KK_compound_univariate_negative_binomial_regression" =
							private$compute_count_KK_compound_univariate_negative_binomial_inference(compute_ci, compute_pval),	
					"count_KK_compound_multivariate_negative_binomial_regression" =
							private$compute_count_KK_compound_multivariate_negative_binomial_inference(compute_ci, compute_pval),
					"count_KK_multivariate_negative_binomial_regression_with_matching_dummies" =
							private$compute_count_KK_multivariate_with_matching_dummies_negative_binomial_inference(compute_ci, compute_pval),
					"count_KK_multivariate_negative_binomial_regression_with_random_intercepts_for_matches" =
							private$compute_count_KK_multivariate_negative_binomial_with_random_intercepts_for_matches_inference(compute_ci, compute_pval),
					########################################### SURVIVAL
					"survival_simple_median_difference" =
							private$compute_survival_simple_median_inference(compute_ci, compute_pval),	
					"survival_simple_restricted_mean_difference" = 
							private$compute_survival_simple_restricted_mean_inference(compute_ci, compute_pval),
					"survival_univariate_weibull_regression" =
							private$compute_survival_univariate_weibull_regression_inference(compute_ci, compute_pval),	
					"survival_multivariate_weibull_regression" =
							private$compute_survival_multivariate_weibull_regression_inference(compute_ci, compute_pval),
					"survival_KK_compound_univariate_weibull_regression" =
							private$compute_survival_KK_compound_univariate_weibull_inference(compute_ci, compute_pval),	
					"survival_KK_compound_multivariate_weibull_regression" =
							private$compute_survival_KK_compound_multivariate_weibull_inference(compute_ci, compute_pval),
					"survival_KK_multivariate_weibull_regression_with_matching_dummies" =
							private$compute_survival_multivariate_with_matching_dummies_weibull_regression_inference(compute_ci, compute_pval),					
					"survival_univariate_coxph_regression" =
							private$compute_survival_univariate_coxph_regression_inference(compute_ci, compute_pval),	
					"survival_multivariate_coxph_regression" =
							private$compute_survival_multivariate_coxph_regression_inference(compute_ci, compute_pval),		
					"survival_KK_multivariate_coxph_regression_with_matching_dummies" =
							private$compute_survival_multivariate_with_matching_dummies_coxph_regression_inference(compute_ci, compute_pval),		
					"survival_KK_multivariate_coxph_regression_with_random_intercepts_for_matches" =
							private$compute_survival_multivariate_with_random_intercepts_coxph_regression_inference(compute_ci, compute_pval)	
				)
			}
			private$inference_model
		},
		
		
		
		
		#takes care of aliases without duplication
		simplified_estimate_type = function(){
		  	if (grepl("simple_mean_difference", self$estimate_type)){
		    	"simple_mean_difference"
		  	} else if (grepl("KK_compound_mean_difference", self$estimate_type)){
		    	"KK_compound_mean_difference"
		  	} else {
				self$estimate_type
			}
		},
		
		
		
		
		compute_incidence_KK_multivariate_logistic_regression_with_matching_dummies_inference = function(){
			tryCatch({
				logistic_regr_mod = suppressWarnings(glm(private$seq_des_obj_priv_int$y ~ ., 
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
						data = cbind(data.frame(y = private$seq_des_obj_priv_int$y, w = private$seq_des_obj_priv_int$w, match_indic = factor(private$match_indic)), private$get_X()),
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
		
		
						compute_count_KK_multivariate_negative_binomial_with_random_intercepts_for_matches_inference = function(){			
			tryCatch({
				mixed_neg_bin_regr_mod = suppressWarnings(lme4::glmer.nb(y ~ . - match_indic + (1 | match_indic), 
								data = cbind(data.frame(y = private$seq_des_obj_priv_int$y, w = private$seq_des_obj_priv_int$w, match_indic = factor(private$match_indic)), private$get_X())))
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
			private$compute_count_negative_binomial_regression(cbind(data.frame(y = private$seq_des_obj_priv_int$y), private$generate_data_frame_with_matching_dummies()))
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
		
		
		compute_survival_multivariate_with_matching_dummies_weibull_regression_inference = function(){
			private$compute_survival_weibull_regression(private$generate_data_frame_with_matching_dummies())
		},
		
		
		compute_survival_multivariate_with_matching_dummies_coxph_regression_inference = function(){
			private$compute_cox_regression(private$generate_data_frame_with_matching_dummies())
		},	
		
		
		compute_survival_multivariate_with_random_intercepts_coxph_regression_inference = function(){
			y = private$seq_des_obj_priv_int$y
			dead = private$seq_des_obj_priv_int$dead
			surv_obj = survival::Surv(y, dead)
			X = private$get_X()
			colnames(X) = paste0("x", 1 : ncol(X)) #as there may be spaces in the colnames which blows the formula up... so damn annoying
			data = cbind(data.frame(w = private$seq_des_obj_priv_int$w, match_indic = factor(private$match_indic)), X)
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
		}