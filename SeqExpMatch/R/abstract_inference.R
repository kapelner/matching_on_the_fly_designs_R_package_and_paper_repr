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
		#' @param thin		For internal use only. Do not specify. You can thank R6's single constructor-only for this coding noise.
		#' @return A new `SeqDesignTest` object.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE, thin = FALSE){	
			if (!thin){
				assertClass(seq_des_obj, "SeqDesign")
				assertCount(num_cores, positive = TRUE)
				assertFlag(verbose)
				seq_des_obj$assert_experiment_completed()
				
				private$any_censoring = seq_des_obj$any_censoring()
				private$seq_des_obj_priv_int = seq_des_obj$.__enclos_env__$private
				private$n = seq_des_obj$get_n()
				private$num_cores = num_cores
				private$verbose = verbose
				if (private$verbose){
					cat(paste0("Intialized inference methods for a ", class(seq_des_obj), " design and response type ", response_type, ".\n"))
				}
			}
		},		
		
		#' @description
		#' Computes a 1-alpha level frequentist confidence interval for the randomization test
		#' 
		#' Here we invert the randomization test that tests the strong null H_0: y_T_i - y_C_i = delta <=> (y_T_i - delta) - y_C_i = 0 so 
		#' we adjust the treatment responses downward by delta. We then find the set of all delta values that is above 1 - alpha/2 (i.e. two-sided)
		#' This is accomplished via a bisection algorithm (algorithm 1 of Glazer and Stark, 2025 available at
		#' https://arxiv.org/abs/2405.05238). These confidence intervals are exact to within tolerance \code{pval_epsilon}.
		#' As far as we know, this only works for response type continuous and survival uncensored (where the CI is on log mean).
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
					lower_upper_ci_bounds = private$compute_mle_or_km_based_confidence_interval(alpha / 100) #ensure a wider CI to be the starting position then pare down
					c(
						private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test, 
							l = lower_upper_ci_bounds[1], 
							u = self$compute_treatment_estimate(), 
							pval_th = alpha / 2, 
							tol = pval_epsilon, 
							log_responses = FALSE, 
							lower = TRUE),
						private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test, 
							l = self$compute_treatment_estimate(), 
							u = lower_upper_ci_bounds[2],  
							pval_th = alpha / 2, 
							tol = pval_epsilon, 
							log_responses = FALSE, 
							lower = FALSE)
					)					
				},
				incidence =  stop("Confidence intervals are not supported for randomization tests for mean difference in incidence outomes"),
				count =      stop("Confidence intervals are not supported for randomization tests for mean difference in count outomes"),
				proportion = stop("Confidence intervals are not supported for randomization tests for mean difference in proportion outomes"),
				survival =   {
					assertNoCensoring(private$any_censoring)
					lower_upper_ci_bounds = private$compute_mle_or_km_based_confidence_interval(alpha / 100) #ensure a wider CI to be the starting position then pare down
					c(
						private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test, 
							l = lower_upper_ci_bounds[1], 
							u = self$compute_treatment_estimate(), 
							pval_th = alpha / 2, 
							tol = pval_epsilon,
							log_responses = TRUE, 
							lower = TRUE),
						private$compute_ci_by_inverting_the_randomization_test_iteratively(nsim_exact_test, 
							l = self$compute_treatment_estimate(), 
							u = lower_upper_ci_bounds[2], 
							pval_th = alpha / 2, 
							tol = pval_epsilon, 
							log_responses = TRUE, 
							lower = FALSE)
					)						
				} #compute_ci_by_inverting_the_randomization_test_iteratively = function(nsim_exact_test, num_cores, l, u, pval_th, tol, log_responses, lower)			
			)
		},		
		
		#' @description
		#' Under the sharp null of 
		#' forall i H_0: y_i_T - y_i_C = delta 
		#' there will be a distribution of the estimates of the treatment effect (over many realizations of assignments)
		#'
		#' @param nsim_exact_test		The number of randomization vectors to use. The default is 501.
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#' @param log_responses			Work in the log space of the responses. This is mostly an internal parameter set to TRUE only if inferring a survival response.
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
		compute_beta_hat_T_randomization_distr_under_sharp_null = function(nsim_exact_test = 501, delta = 0, log_responses = FALSE){
			assertNumeric(delta)
			assertCount(nsim_exact_test, positive = TRUE)
			
			is_KK =          is(self, "SeqDesignInferenceMLEorKMKK")
			is_KK_compound = is(self, "SeqDesignInferenceAllKKCompoundMeanDiff")
			is_Efron = !is.null(private$seq_des_obj_priv_int$weighted_coin_prob)

			X = private$X	
			if (delta != 0){
				if (private$seq_des_obj_priv_int$response_type != "continous"){
					stop("randomization tests with delta nonzero only works for continuous type!!!!")
				}
				y = if (log_responses){
						log(copy(private$seq_des_obj_priv_int$y)) #copy to ensure we don't edit it
					} else {
						copy(private$seq_des_obj_priv_int$y) #copy to ensure we don't edit it
					}
				#we are testing against H_0: y_T_i - y_C_i = delta <=> (y_T_i - delta) - y_C_i = 0 
				#so adjust the treatment responses downward by delta
				y[private$seq_des_obj_priv_int$w == 1] = y[private$seq_des_obj_priv_int$w == 1] - delta
			} else {
				y = if (log_responses){
						log(private$seq_des_obj_priv_int$y) #no need for copy as we are not mutating
					} else {
						private$seq_des_obj_priv_int$y #no need for copy as we are not mutating
					}				
			}
			dead = private$seq_des_obj_priv_int$dead
			w = private$seq_des_obj_priv_int$w
			t = private$seq_des_obj_priv_int$t
			prob_T = private$seq_des_obj_priv_int$prob_T
			p_raw_t = private$seq_des_obj_priv_int$p_raw_t
			if (is_Efron){
				weighted_coin_prob = private$seq_des_obj_priv_int$weighted_coin_prob
			}
			
			seq_des_obj_priv_int = private$seq_des_obj_priv_int	
			#make a copy of the object and then permute the allocation vector according to the design
			seq_des_r = seq_des_obj_priv_int$thin_duplicate()
			#set only the data we need for all designs
			seq_des_r$.__enclos_env__$private$y = y #set the new responses
			seq_des_r$.__enclos_env__$private$X = X	
			seq_des_r$.__enclos_env__$private$dead = dead
			seq_des_r$.__enclos_env__$private$w = w
			seq_des_r$.__enclos_env__$private$t = t
			seq_des_r$.__enclos_env__$private$prob_T = prob_T
			seq_des_r$.__enclos_env__$private$p_raw_t = p_raw_t
			if (is_Efron){ #for Efron, one more piece of data is needed
				seq_des_r$.__enclos_env__$private$weighted_coin_prob = weighted_coin_prob
			}
			#now initialize a new inference object and compute beta_T_hat
			seq_inf_r = private$thin_duplicate()	
			seq_inf_r$.__enclos_env__$private$X = X			
										
			if (private$num_cores == 1){ #easier on the OS I think...
				beta_hat_T_diff_ws = array(NA, nsim_exact_test)
				
				for (r in 1 : nsim_exact_test){
					#cat("		r =", r, "/", nsim_exact_test, "\n")	
					#scramble the allocation vector based on the design algorithm
					seq_des_r$.__enclos_env__$private$redraw_w_according_to_design()
					#set the internals of the design data to the inference object
					seq_inf_r$.__enclos_env__$private$seq_des_obj_priv_int = seq_des_r$.__enclos_env__$private
					seq_inf_r$.__enclos_env__$private$cached_values = list() #ensure nothing is kept between iterations
					if (is_KK){
						seq_inf_r$.__enclos_env__$private$compute_basic_match_data()
					}	
					if (is_KK_compound){
						seq_inf_r$.__enclos_env__$private$compute_reservoir_and_match_statistics()
					}					
					beta_hat_T_diff_ws[r] = seq_inf_r$compute_treatment_estimate()
				}
				#print(ggplot2::ggplot(data.frame(sims = beta_hat_T_diff_ws)) + ggplot2::geom_histogram(ggplot2::aes(x = sims), bins = 50))
			} else {	
				cl = doParallel::makeCluster(private$num_cores)
				doParallel::registerDoParallel(cl)	
				#now copy them to each core's memory
				doParallel::clusterExport(cl, list("seq_des_r", "seq_inf_r", "is_KK", "is_KK_compound", "is_Efron", "seq_des_obj_priv_int", "y", "X", "dead", "w", "t", "prob_T", "p_raw_t", "weighted_coin_prob"), envir = environment())
				#now do the parallelization
				beta_hat_T_diff_ws = doParallel::foreach(r = 1 : nsim_exact_test, .inorder = FALSE, .combine = c) %dopar% {
					seq_des_r$.__enclos_env__$private$redraw_w_according_to_design()
					seq_inf_r$.__enclos_env__$private$seq_des_obj_priv_int = seq_des_r$.__enclos_env__$private
					seq_inf_r$.__enclos_env__$private$cached_values = list() #ensure nothing is kept between iterations
					if (is_KK){ #for matching-on-the-fly there is some more data required for inference
						seq_inf_r$.__enclos_env__$private$compute_basic_match_data()
					}		
					if (is_KK_compound){
						seq_inf_r$.__enclos_env__$private$compute_reservoir_and_match_statistics()
					}					
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
		#' @param log_responses			Work in the log space of the responses. This is mostly an internal parameter set to TRUE only if inferring a survival response.
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
		compute_two_sided_pval_for_treatment_effect_rand = function(nsim_exact_test = 501, delta = 0, log_responses = FALSE, na.rm = FALSE){
			assertLogical(na.rm)
			#approximate the null distribution by computing estimates on many draws of w
			beta_hat_T_diff_ws = self$compute_beta_hat_T_randomization_distr_under_sharp_null(nsim_exact_test, delta ,log_responses)
			#this calculates the actual estimate to compare against the null distribution
			beta_hat_T = self$compute_treatment_estimate()
			#finally compute the p-value
			if (na.rm){
				nsim_exact_test = sum(!is.na(beta_hat_T))
			}
			2 * min(
				sum(beta_hat_T_diff_ws > beta_hat_T, na.rm = na.rm) / nsim_exact_test, #na.rm is because some runs produce NA... TO-DO is to trace these down
				sum(beta_hat_T_diff_ws < beta_hat_T, na.rm = na.rm) / nsim_exact_test
			)
		}
	),
	
	private = list(
		seq_des_obj_priv_int = NULL, 			#seq_des_obj_priv_int_obj = function(){private$seq_des_obj_priv_int},
		any_censoring = NULL,					#get_any_censoring = function(){private$any_censoring},
		num_cores = NULL,		 				#get_num_cores = function(){private$num_cores},
		verbose = FALSE,		 				#get_verbose = function(){private$verbose},
		n = NULL,		 						#get_n = function(){private$n},
		p = NULL,		 						#get_p = function(){private$p},
		X = NULL,								#get_X is defined later as it needs some logic dependent on the design type
		cached_values = list(),					#get_cached_values = function(){private$cached_values},
				
		thin_duplicate = function(){
			do.call(get(class(self)[1])$new, args = list(thin = TRUE))		
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
		
		compute_ci_by_inverting_the_randomization_test_iteratively = function(nsim_exact_test, l, u, pval_th, tol, log_responses, lower){
			pval_l = private$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test, delta = l, log_responses)
			pval_u = private$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test, delta = u, log_responses)
			repeat {
				if (pval_u - pval_l <= tol & lower){
					return(l)
				}
				if (pval_u - pval_l <= tol & !lower){
					return(u)
				}				
				m = (l + u) / 2
				pval_m = private$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test, delta = m, log_responses)
				if (pval_m >= pval_th & lower){
					u = m
					pval_u = pval_m
				} else if (pval_m >= pval_th & !lower){
					l = m
					pval_l = pval_m
				} else if (lower){
					l = m
					pval_l = pval_m
				} else if (upper){
					u = m
					pval_u = pval_m
				}
			}
		}
	)
)
