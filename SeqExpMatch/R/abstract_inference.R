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
		initialize = function(seq_des_obj, num_cores = 1, verbose = TRUE){
			assertClass(seq_des_obj, "SeqDesign")
			assertCount(num_cores, positive = TRUE)
			assertFlag(verbose)
			seq_des_obj$assert_experiment_completed()
			
			private$any_censoring = seq_des_obj$any_censoring()
			private$seq_des_obj_priv_int = seq_des_obj$.__enclos_env__$private
			private$n = private$seq_des_obj_priv_int$n
			private$yTs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 1]
			private$yCs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 0]
			private$deadTs = private$seq_des_obj_priv_int$dead[private$seq_des_obj_priv_int$w == 1]
			private$deadCs = private$seq_des_obj_priv_int$dead[private$seq_des_obj_priv_int$w == 0]
			private$num_cores = num_cores
			private$verbose = verbose			
#			if (private$verbose){
#				cat(paste0("Intialized inference methods for a ", seq_des_obj$design, " design, response type ", response_type, ", estimation type ", estimate_type, " and test type: ", test_type, ".\n"))
#			}	
		},		
		
		#' @description
		#' Computes a 1-alpha level frequentist confidence interval for the randomization test
		#' 
		#' Here we invert the randomization test that tests the strong null H_0: y_T_i - y_C_i = delta <=> (y_T_i - delta) - y_C_i = 0 so 
		#' we adjust the treatment responses downward by delta. We then find the set of all delta values that is above 1 - alpha/2 (i.e. two-sided)
		#' This is accomplished via a bisection algorithm (algorithm 1 of Glazer and Stark, 2025 available at
		#' https://arxiv.org/abs/2405.05238). These confidence intervals are exact to within tolerance \code{pval_epsilon}.
		#' As far as we know, this only works for response type continuous.
		#' 
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#' @param nsim_exact_test		The number of randomization vectors (applicable for test type "randomization-exact" only). 
		#' 								The default is 1000 providing good resolutions to confidence intervals.
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
		compute_confidence_interval_rand = function(alpha = 0.05, nsim_exact_test = 501, pval_epsilon = 0.001){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertCount(nsim_exact_test, positive = TRUE)
			assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
			
			switch(private$seq_des_obj_priv_int$response_type,
				continuous = {
					#####TO-DO
					lower_upper_ci_bounds = private$compute_mle_or_km_based_confidence_interval(alpha / 100)
					c(
						private$compute_ci_lower_by_inverting_the_randomization_test(nsim_exact_test, private$num_cores, self$compute_treatment_estimate(), lower_upper_ci_bounds[2], alpha / 2, pval_epsilon),
						private$compute_ci_upper_by_inverting_the_randomization_test(nsim_exact_test, private$num_cores, lower_upper_ci_bounds[1], self$compute_treatment_estimate(), alpha / 2, pval_epsilon)
					)					
				},
				incidence =  stop("Confidence intervals are not supported for randomization tests for mean difference in incidence outomes"),
				count =      stop("Confidence intervals are not supported for randomization tests for mean difference in count outomes"),
				proportion = stop("Confidence intervals are not supported for randomization tests for mean difference in proportion outomes"),
				survival =   stop("Confidence intervals are not supported for randomization tests for mean difference in survival outomes")
			
			)
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
		compute_two_sided_pval_for_treatment_effect_rand = function(nsim_exact_test = 501, delta = 0){
			assertNumeric(delta)
			assertCount(nsim_exact_test, positive = TRUE)

			seq_inf_class_constructor = get(class(self)[1])$new 		
			y = private$seq_des_obj_priv_int$y
			X = private$X	
			if (delta != 0){
				if (private$seq_des_obj_priv_int$response_type != "continous"){
					stop("randomization tests with delta nonzero only works for continuous type!!!!")
				}
				#we are testing against H_0: y_T_i - y_C_i = delta <=> (y_T_i - delta) - y_C_i = 0 so adjust the treatment responses downward by delta
				y[private$seq_des_obj_priv_int$w == 1] = y[private$seq_des_obj_priv_int$w == 1] - delta
			}
			
			seq_des_obj_priv_int = private$seq_des_obj_priv_int
			if (private$num_cores == 1){
				#easier on the OS I think...
				b_T_sims = array(NA, nsim_exact_test)
				for (r in 1 : nsim_exact_test){
					#cat("		r =", r, "/", nsim_exact_test, "\n")
					#make a copy of the object and then permute the allocation vector according to the design
					seq_des_r = seq_des_obj_priv_int$duplicate()
					seq_des_r$.__enclos_env__$private$y = y #set the new responses
					seq_des_r$.__enclos_env__$private$X = X	
					seq_des_r$.__enclos_env__$private$redraw_w_according_to_design()
					#now initialize a new inference object and compute beta_T_hat
					seq_inf_r = do.call(seq_inf_class_constructor, args = list(seq_des = seq_des_r, verbose = FALSE))
					seq_inf_r$.__enclos_env__$private$X = X					
					b_T_sims[r] = seq_inf_r$compute_treatment_estimate()
				}
				#print(ggplot2::ggplot(data.frame(sims = b_T_sims)) + ggplot2::geom_histogram(ggplot2::aes(x = sims), bins = 50))
			} else {	
				cl = doParallel::makeCluster(private$num_cores)
				doParallel::registerDoParallel(cl)	
				#now copy them to each core's memory
				doParallel::clusterExport(cl, list("seq_des_obj_priv_int", "y", "seq_inf_class_constructor", "X"), envir = environment())
				#now do the parallelization
				b_T_sims = doParallel::foreach(r = 1 : nsim_exact_test, .inorder = FALSE, .combine = c) %dopar% {
					#make a copy of the object and then permute the allocation vector according to the design
					seq_des_r = seq_des_obj_priv_int$duplicate()
					seq_des_r$.__enclos_env__$private$y = y #set the new responses
					seq_des_r$.__enclos_env__$private$X = X
					seq_des_r$.__enclos_env__$private$redraw_w_according_to_design()				
					seq_inf_r = do.call(seq_inf_class_constructor, args = list(seq_des = seq_des_r, verbose = FALSE))
					seq_inf_r$.__enclos_env__$private$X = X
					seq_inf_r$compute_treatment_estimate()			
				}
				doParallel::stopCluster(cl)
				rm(cl); gc()
			}
			#this calculates the two-sided pval
			beta_hat_T = self$compute_treatment_estimate()
			#finally compute the p-value
			2 * min(
				sum(b_T_sims > beta_hat_T) / nsim_exact_test, 
				sum(b_T_sims < beta_hat_T) / nsim_exact_test
			)
		}		

	),
	
	private = list(
		seq_des_obj_priv_int = NULL, 			seq_des_obj_priv_int_obj = function(){private$seq_des_obj_priv_int},
		any_censoring = NULL,					get_any_censoring = function(){private$any_censoring},
		num_cores = NULL,		 				get_num_cores = function(){private$num_cores},
		verbose = FALSE,		 				get_verbose = function(){private$verbose},
		n = NULL,		 						get_n = function(){private$n},
		p = NULL,		 						get_p = function(){private$p},
		X = NULL,								#get_X is defined later as it needs some logic dependent on the design type
		yTs = NULL,		 						get_yTs = function(){private$yTs},
		yCs = NULL,		 						get_yCs = function(){private$yCs},
		deadTs = NULL,							get_deadTs = function(){private$deadTs},
		deadCs = NULL,							get_deadCs = function(){private$deadCs},
		cached_values = list(),					get_cached_values = function(){private$cached_values},
		
		get_X = function(){
			if (is.null(private$X)){
				if (is.null(private$seq_des_obj_priv_int$X)){
					private$seq_des_obj_priv_int$covariate_impute_if_necessary_and_then_create_model_matrix()
				}
	#			if (seq_des_obj$design %in% c("CRD", "iBCRD", "Efron")){ #imputations were never done yet
	#				seq_des_obj$.__enclos_env__$private$covariate_impute_if_necessary_and_then_create_model_matrix()
	#			}
	#			if (private$isKK & uniqueN(private$match_indic) == 1){ #imputations were never done yet either
	#				seq_des_obj$.__enclos_env__$private$covariate_impute_if_necessary_and_then_create_model_matrix()
	#			}	
				private$X = private$seq_des_obj_priv_int$compute_all_subject_data()$X_all
			}
			private$X
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
		
		######### MODEL INFERENCE FUNCTIONS	
		
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
					"proportion_simple_logodds_regression" =
							private$compute_proportion_univariate_beta_regression_inference(compute_ci, compute_pval),
					"proportion_multivariate_beta_regression" =
							private$compute_proportion_multivariate_beta_regression_inference(compute_ci, compute_pval),
					"proportion_KK_compound_univariate_beta_regression" =
							private$compute_proportion_KK_compound_univariate_beta_regression_inference(compute_ci, compute_pval),
					"proportion_KK_compound_multivariate_beta_regression" =
							private$compute_proportion_KK_compound_multivariate_beta_regression_inference(compute_ci, compute_pval),
					"proportion_KK_multivariate_beta_regression_with_matching_dummies" =
							private$compute_proportion_KK_multivariate_beta_regression_with_matching_dummies_inference(compute_ci, compute_pval),
					########################################### COUNT
					"count_univariate_negative_binomial_regression" =
							private$compute_count_univariate_negative_binomial_inference(compute_ci, compute_pval),
					"count_multivariate_negative_binomial_regression" =
							private$compute_count_multivariate_negative_binomial_inference(compute_ci, compute_pval),
					"count_KK_compound_univariate_negative_binomial_regression" =
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
		
		
				
		compute_ci_lower_by_inverting_the_randomization_test = function(nsim_exact_test, l, u, pval_th, tol){
			return(NA) #not implemented yet!!!!
			################################################################################################
			pval_l = private$compute_randomization_test_p_val(nsim_exact_test, l)
			pval_u = private$compute_randomization_test_p_val(nsim_exact_test, u)
			repeat {
				if (pval_u - pval_l <= tol){
					return(l)
				}
				m = (l + u) / 2
				pval_m = private$compute_randomization_test_p_val(nsim_exact_test, m)
				if (pval_m >= pval_th){
					u = m
					pval_u = pval_m
				} else {
					l = m
					pval_l = pval_m
				}
			}
		},
		
		compute_ci_upper_by_inverting_the_randomization_test = function(nsim_exact_test, l, u, pval_th, tol){
			return(NA) #not implemented yet!!!!
			################################################################################################
			pval_l = private$compute_randomization_test_p_val(nsim_exact_test, l)
			pval_u = private$compute_randomization_test_p_val(nsim_exact_test, u)
			repeat {
				if (pval_l - pval_u <= tol){
					return(u)
				}
				m = (l + u) / 2
				pval_m = private$compute_randomization_test_p_val(nsim_exact_test, m)
				if (pval_m >= pval_th){
					l = m
					pval_l = pval_m
				} else {
					u = m
					pval_u = pval_m
				}
			}
		},

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		

		
		

		

		

		

		

		
		compute_incidence_KK_compound_univariate_logistic_regression_inference = function(){
			stop("not implemented yet")
		},
		
		compute_incidence_KK_compound_multivariate_logistic_regression_inference = function(){
			stop("not implemented yet")
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
		
		compute_proportion_KK_multivariate_beta_regression_with_matching_dummies_inference = function(){
			private$shared_beta_regression_inference(cbind(data.frame(y = private$seq_des_obj_priv_int$y), private$generate_data_frame_with_matching_dummies()))
		},


		
		compute_count_univariate_negative_binomial_inference = function(){
			private$compute_count_negative_binomial_regression(data.frame(y = private$seq_des_obj_priv_int$y, w = private$seq_des_obj_priv_int$w))
		},
		
		compute_count_multivariate_negative_binomial_inference = function(){
			private$compute_count_negative_binomial_regression(cbind(data.frame(y = private$seq_des_obj_priv_int$y, w = private$seq_des_obj_priv_int$w), private$get_X()))
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
		
		compute_survival_simple_median_inference = function(compute_ci, compute_pval, B = NULL){			
			if (is.null(private$beta_hat_T)){
				survival_obj = survival::Surv(private$seq_des_obj_priv_int$y, private$seq_des_obj_priv_int$dead)
				survival_fit_obj = survival::survfit(survival_obj ~ private$seq_des_obj_priv_int$w)
				survival_fit_res = summary(survival_fit_obj)$table
				private$beta_hat_T = survival_fit_res[2, 7] - survival_fit_res[1, 7]
			}
			
			if (compute_ci){
				if (is.null(B)){
					B = private$default_B_for_median_inference
				}
				assertCount(B, positive = TRUE)
				if (B != private$previous_B_value){
					private$ci_mle = NULL
				}
				
				if (is.null(private$ci_mle)){
					test_obj = suppressWarnings(controlTest::quantileControlTest(private$yTs, private$deadTs, private$yCs, private$deadCs, B = B))
					list(
						mod = survival_fit_obj,
						summary_table = survival_fit_res,
						beta_hat_T = beta_hat_T,
						s_beta_hat_T = abs(beta_hat_T / test_obj$Z), #wtf is quantileControlTest's se field if not this??????????
						is_z = TRUE
					)
				}
			}
			
			if (compute_pval & is.null(private$pval_mle)){
				survival_diff_obj = survdiff(Surv(y, 1 - c_vec) ~ is_male)
				private$pval_mle = survival_diff_obj$pvalue
			}
		},
		
		compute_survival_simple_restricted_mean_inference = function(){
			survival_obj = survival::Surv(private$seq_des_obj_priv_int$y, private$seq_des_obj_priv_int$dead)
			survival_fit_obj = survival::survfit(survival_obj ~ private$seq_des_obj_priv_int$w)
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
			private$compute_survival_weibull_regression(private$seq_des_obj_priv_int$w)
		},
		
		compute_survival_multivariate_weibull_regression_inference = function(){
			private$compute_survival_weibull_regression(cbind(private$seq_des_obj_priv_int$w, private$get_X()))
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
			surv_regr_mod = robust_survreg(private$seq_des_obj_priv_int$y, private$seq_des_obj_priv_int$dead, data_obj)
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
			private$compute_cox_regression(data.frame(w = private$seq_des_obj_priv_int$w))
		},
		
		compute_survival_multivariate_coxph_regression_inference = function(){
			private$compute_cox_regression(cbind(data.frame(w = private$seq_des_obj_priv_int$w), private$get_X()))
		},
		
		compute_survival_multivariate_with_matching_dummies_coxph_regression_inference = function(){
			private$compute_cox_regression(private$generate_data_frame_with_matching_dummies())
		},
		
		compute_cox_regression = function(data_obj){
			y = private$seq_des_obj_priv_int$y
			dead = private$seq_des_obj_priv_int$dead
			
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
	)
)
