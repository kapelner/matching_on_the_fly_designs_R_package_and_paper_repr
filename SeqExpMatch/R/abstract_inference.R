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
			private$n = seq_des_obj$get_n()
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
		#' Under the sharp null of 
		#' forall i H_0: y_i_T - y_i_C = delta 
		#' there will be a distribution of the estimates of the treatment effect (over many realizations of assignments)
		#'
		#' @param nsim_exact_test		The number of randomization vectors to use. The default is 501.
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#' 
		#' @return 	A vector of size \code{nsim_exact_test} that has the values of beta_hat_T over many w draws.
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
		#' beta_hat_T_diff_ws = seq_des_inf$compute_beta_hat_T_randomization_distr_under_sharp_null()
		#' ggplot(data.frame(beta_hat_T_diff_ws = beta_hat_T_diff_ws)) + geom_histogram(aes(x = beta_hat_T_diff_ws))
		compute_beta_hat_T_randomization_distr_under_sharp_null = function(nsim_exact_test = 501, delta = 0){
			assertNumeric(delta)
			assertCount(nsim_exact_test, positive = TRUE)

			seq_inf_class_constructor = get(class(self)[1])$new 		
			X = private$X	
			if (delta != 0){
				if (private$seq_des_obj_priv_int$response_type != "continous"){
					stop("randomization tests with delta nonzero only works for continuous type!!!!")
				}
				y = copy(private$seq_des_obj_priv_int$y) #copy to ensure we don't edit it
				#we are testing against H_0: y_T_i - y_C_i = delta <=> (y_T_i - delta) - y_C_i = 0 
				#so adjust the treatment responses downward by delta
				y[private$seq_des_obj_priv_int$w == 1] = y[private$seq_des_obj_priv_int$w == 1] - delta
			} else {
				y = private$seq_des_obj_priv_int$y #no need for copy as we are not mutating
			}
			
			seq_des_obj_priv_int = private$seq_des_obj_priv_int
			if (private$num_cores == 1){
				#easier on the OS I think...
				beta_hat_T_diff_ws = array(NA, nsim_exact_test)
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
					beta_hat_T_diff_ws[r] = seq_inf_r$compute_treatment_estimate()
				}
				#print(ggplot2::ggplot(data.frame(sims = beta_hat_T_diff_ws)) + ggplot2::geom_histogram(ggplot2::aes(x = sims), bins = 50))
			} else {	
				cl = doParallel::makeCluster(private$num_cores)
				doParallel::registerDoParallel(cl)	
				#now copy them to each core's memory
				doParallel::clusterExport(cl, list("seq_des_obj_priv_int", "y", "seq_inf_class_constructor", "X"), envir = environment())
				#now do the parallelization
				beta_hat_T_diff_ws = doParallel::foreach(r = 1 : nsim_exact_test, .inorder = FALSE, .combine = c) %dopar% {
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
		#' @param na.rm 				Should we remove beta_hat_T's that are NA's? Default is \code{FALSE}.
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
		compute_two_sided_pval_for_treatment_effect_rand = function(nsim_exact_test = 501, delta = 0, na.rm = FALSE){
			assertLogical(na.rm)
			
			beta_hat_T_diff_ws = self$compute_beta_hat_T_randomization_distr_under_sharp_null(nsim_exact_test, delta)
			#this calculates the two-sided pval
			beta_hat_T = self$compute_treatment_estimate()
			#finally compute the p-value
			if (na.rm){
				nsim_exact_test = sum(!is.na(beta_hat_T))
			}
			2 * min(
				sum(beta_hat_T_diff_ws > beta_hat_T, na.rm = na.rm) / nsim_exact_test, 
				sum(beta_hat_T_diff_ws < beta_hat_T, na.rm = na.rm) / nsim_exact_test
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
		}
	)
)
