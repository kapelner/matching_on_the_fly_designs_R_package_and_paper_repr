#' Inference for A Sequential Design
#'
#' @description
#' An abstract R6 Class that provides MLE-based tests and intervals for a treatment effect in a sequential design.
#' 
SeqDesignInferenceMLEorKM = R6::R6Class("SeqDesignInferenceMLEorKM",
	inherit = SeqDesignInference,
	public = list(
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
        #' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference 
		#' 							(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 							for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		#'				
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		},
		
		#' @description
		#' Creates the boostrap distribution of the estimate for the treatment effect
		#' 
		#' @param B						Number of bootstrap samples. The default is 501.
		#' 
		#' @return 	A vector of length \code{B} with the bootstrap values of the estimates of the treatment effect
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
		#' beta_hat_T_bs = seq_des_inf$approximate_boostrap_distribution_beta_hat_T()
		#' ggplot(data.frame(beta_hat_T_bs = beta_hat_T_bs)) + geom_histogram(aes(x = beta_hat_T_bs))
		#' 			
		approximate_boostrap_distribution_beta_hat_T = function(B = 501){
			assertCount(B, positive = TRUE)
									
			n = private$seq_des_obj_priv_int$n	
			y = private$seq_des_obj_priv_int$y
			dead = private$seq_des_obj_priv_int$dead
			w = private$seq_des_obj_priv_int$w
			X = private$get_X()
			#now duplicate the design and the inference objects so we can set new data within them for each iteration
			seq_des_r = private$seq_des_obj_priv_int$duplicate()				
			seq_inf_r = private$duplicate()
			
			if (private$num_cores == 1){ #easier on the OS I think...
				beta_hat_T_bs = array(NA, B)
				for (r in 1 : B){
					#draw a bootstrap sample
					i_b = sample_int_replace_cpp(n, n)
					seq_des_r$.__enclos_env__$private$y = y[i_b]
					seq_des_r$.__enclos_env__$private$dead = dead[i_b]
					seq_des_r$.__enclos_env__$private$X = X[i_b, ]
					seq_des_r$.__enclos_env__$private$w = w[i_b]
					#compute beta_T_hat
					seq_inf_r$.__enclos_env__$private$seq_des_obj_priv_int = seq_des_r$.__enclos_env__$private
					seq_inf_r$.__enclos_env__$private$cached_values = list() #ensure nothing is kept between iterations
					beta_hat_T_bs[r] = seq_inf_r$compute_treatment_estimate()
				}
				#print(ggplot2::ggplot(data.frame(sims = beta_hat_T_bs)) + ggplot2::geom_histogram(ggplot2::aes(x = sims), bins = 50))
			} else {	
				cl = doParallel::makeCluster(private$num_cores)
				doParallel::registerDoParallel(cl)
				#now copy them to each core's memory
				doParallel::clusterExport(cl, list("seq_des_r", "seq_inf_r", "n", "y", "dead", "X", "w"), envir = environment())
				#now do the parallelization
				beta_hat_T_bs = doParallel::foreach(r = 1 : B, .inorder = FALSE, .combine = c) %dopar% {
					#draw a bootstrap sample
					i_b = sample_int_replace_cpp(n, n)
					seq_des_r$.__enclos_env__$private$y = y[i_b]
					seq_des_r$.__enclos_env__$private$dead = dead[i_b]
					seq_des_r$.__enclos_env__$private$X = X[i_b, ]
					seq_des_r$.__enclos_env__$private$w = w[i_b]
					#compute beta_T_hat
					seq_inf_r$.__enclos_env__$private$seq_des_obj_priv_int = seq_des_r$.__enclos_env__$private
					seq_inf_r$.__enclos_env__$private$cached_values = list() #ensure nothing is kept between iterations
					seq_inf_r$compute_treatment_estimate()			
				}
				doParallel::stopCluster(cl)
				rm(cl); gc()
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
		#' seq_des_inf$compute_bootstrap_two_sided_pval()
		#' 					
		compute_bootstrap_two_sided_pval = function(delta = 0, B = 501, na.rm = FALSE){
			assertNumeric(delta)
			assertLogical(na.rm)

			beta_hat_T_bs = private$get_or_cache_bootstrap_samples(B)
			
			2 * min(
				mean(delta < beta_hat_T_bs, na.rm = na.rm), 
				mean(delta > beta_hat_T_bs, na.rm = na.rm)
			)
		}		
	),
	private = list(		
		get_or_cache_bootstrap_samples = function(B){			
			if (is.null(private$cached_values$beta_hat_T_bs)){
				private$cached_values$beta_hat_T_bs = self$approximate_boostrap_distribution_beta_hat_T(B)
			} else {
				B_0 = length(private$cached_values$beta_hat_T_bs)
				if (B_0 > B){
					return (private$cached_values$beta_hat_T_bs[1 : B]) #send back what we need but don't reduce the cache
				} else if (B_0 < B){ 
					private$cached_values$beta_hat_T_bs = c(private$cached_values$beta_hat_T_bs, #go get more and add them to the cache in case we need them later
						self$approximate_boostrap_distribution_beta_hat_T(B - B_0))
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
						pnorm(z_or_t_stats)
					} else {
						pt(z_or_t_stats, private$cached_values$df)
					}
			2 * min(probs)
		}	
	)
)
		
		