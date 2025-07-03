#' Inference based on Maximum Likelihood for KK designs  
#'
#' @description
#' An abstract class
#' 
#'
SeqDesignInferenceMLEorKMKK = R6::R6Class("SeqDesignInferenceMLEorKMKK",
	inherit = SeqDesignInferenceMLEorKM,
	public = list(
		
        #' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference 
		#' 							(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 							for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		#'
		initialize = function(seq_des_obj, num_cores = 1, verbose = TRUE){
			if (!seq_des_obj$.__enclos_env__$private$isKK){
				stop("This type of inference is only available for KK designs.")
			}
			super$initialize(seq_des_obj, num_cores, verbose)	
			private$helper = SeqDesignInferenceHelperKK$new(seq_des_obj, private$get_X())
			private$KKstats = private$helper$get_post_matching_data()
		},
		
			
		#' @description
		#' Computes a 1-alpha level frequentist bootstrap confidence interval for KK tests
		#' 
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#' @param B						Number of bootstrap samples. The default is NA which corresponds to B=501.
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
		compute_bootstrap_confidence_interval = function(alpha = 0.05, B = 501){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertCount(B, positive = TRUE)

			n = private$seq_des_obj_priv_int$n	
			y = private$seq_des_obj_priv_int$y
			dead = private$seq_des_obj_priv_int$dead
			w = private$seq_des_obj_priv_int$w
			X = private$get_X()
			match_indic = private$KKstats$match_indic
			i_reservoir = which(match_indic == 0)
			n_reservoir = length(i_reservoir)
			m = private$KKstats$m
			seq_inf_class_constructor = get(class(self)[1])$new 
			match_indic_b = c(rep(0, n_reservoir), rep(1 : m, each = 2))
			
			if (super$get_num_cores() == 1){ #easier on the OS I think...
				b_T_sims = array(NA, B)
				for (r in 1 : B){
					#draw a bootstrap sample of the reservoir
					i_reservoir_b = sample(i_reservoir, n_reservoir, replace = TRUE)
					ms_b = sample_int_replace_cpp(m, m)
					i_b = c(i_reservoir_b, array(NA, m * 2))
					i_b_idx = n_reservoir
					for (m0 in 1 : m){
						i_b[(i_b_idx + 1) : (i_b_idx + 2)] = which(match_indic == ms_b[m0])
						i_b_idx = i_b_idx + 2
					}
					
					seq_des_r = private$seq_des_obj_priv_int$duplicate()
					seq_des_r$.__enclos_env__$private$y = y[i_b]
					seq_des_r$.__enclos_env__$private$dead = dead[i_b]
					seq_des_r$.__enclos_env__$private$X = X[i_b, ]
					seq_des_r$.__enclos_env__$private$w = w[i_b]
					seq_des_r$.__enclos_env__$private$match_indic = match_indic_b
					#compute beta_T_hat					
					seq_inf_r = do.call(seq_inf_class_constructor, args = list(seq_des = seq_des_r, verbose = FALSE))
					seq_inf_r$.__enclos_env__$private$X = seq_des_r$.__enclos_env__$private$X		
					b_T_sims[r] = seq_inf_r$compute_treatment_estimate()
				}
				#print(ggplot2::ggplot(data.frame(sims = b_T_sims)) + ggplot2::geom_histogram(ggplot2::aes(x = sims), bins = 50))
			} else {	
				cl = doParallel::makeCluster(super$get_num_cores())
				doParallel::registerDoParallel(cl)			
				#now copy them to each core's memory
				doParallel::clusterExport(cl, list("seq_des_obj", "n", "match_indic", "i_reservoir", "n_reservoir", "match_indic_b", "m", "y", "dead", "X", "w", "seq_inf_class_constructor"), envir = environment())
				#now do the parallelization
				b_T_sims = doParallel::foreach(r = 1 : B, .inorder = FALSE, .combine = c) %dopar% {
					#draw a bootstrap sample of the reservoir
					i_reservoir_b = sample(i_reservoir, n_reservoir, replace = TRUE)
					ms_b = sample_int_replace_cpp(m, m)
					i_b = c(i_reservoir_b, array(NA, m * 2))
					i_b_idx = n_reservoir
					for (m0 in 1 : m){
						i_b[(i_b_idx + 1) : (i_b_idx + 2)] = which(match_indic == ms_b[m0])
						i_b_idx = i_b_idx + 2
					}
					
					seq_des_r = private$seq_des_obj_priv_int$duplicate()
					seq_des_r$.__enclos_env__$private$y = y[i_b]
					seq_des_r$.__enclos_env__$private$dead = dead[i_b]
					seq_des_r$.__enclos_env__$private$X = X[i_b, ]
					seq_des_r$.__enclos_env__$private$w = w[i_b]
					seq_des_r$.__enclos_env__$private$match_indic = match_indic_b
					#compute beta_T_hat					
					seq_inf_r = do.call(seq_inf_class_constructor, args = list(seq_des = seq_des_r, verbose = FALSE))
					seq_inf_r$.__enclos_env__$private$X = seq_des_r$.__enclos_env__$private$X		
					seq_inf_r$compute_treatment_estimate()			
				}
				doParallel::stopCluster(cl)
				rm(cl); gc()
			}
			
			#this finally computes the ci
			quantile(b_T_sims, c(alpha / 2, 1 - alpha / 2))
		}			
		
		
	),
	private = list(
		helper = NULL,
		KKstats = NULL,
		
		compute_model_matrix_with_matching_dummies = function(){
			private$helper$data_frame_with_matching_dummies()
		}	
	)
)
