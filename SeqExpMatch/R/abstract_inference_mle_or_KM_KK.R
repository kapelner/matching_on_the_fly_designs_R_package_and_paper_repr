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
		#' @param thin		For internal use only. Do not specify. You can thank R6's single constructor-only for this coding noise.
		#'
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE, thin = FALSE){
			if (!thin){
				assertClass(seq_des_obj, "SeqDesignKK14")
				super$initialize(seq_des_obj, num_cores, verbose)
				private$setup_matching_data(seq_des_obj, private$get_X())				
			}
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
			match_indic = private$KKstats$match_indic
			i_reservoir = which(match_indic == 0)
			n_reservoir = length(i_reservoir)
			m = private$KKstats$m
			match_indic_b = c(rep(0, n_reservoir), rep(1 : m, each = 2))
			
			if (private$num_cores == 1){ #easier on the OS I think...
				beta_hat_T_bs = array(NA, B)
				for (r in 1 : B){
					#draw a bootstrap of both the reservoir and matches - first create index vector					
					i_b = array(NA, n_reservoir + 2 * m)
					#draw a bootstrap sample of the reservoir
					i_reservoir_b_i = sample_int_replace_cpp(n_reservoir, n_reservoir)
					#load it into the index vector
					i_b[1 : n_reservoir] = i_reservoir[i_reservoir_b_i]
					#draw a bootstrap sample of the matches
					ms_b = sample_int_replace_cpp(m, m)
					#load it into the index vector 2 by 2
#					i_b_idx = n_reservoir
#					for (m0 in 1 : m){
#						i_b[(i_b_idx + 1) : (i_b_idx + 2)] = which(match_indic == ms_b[m0])
#						i_b_idx = i_b_idx + 2
#					}
					fill_i_b_with_matches_loop_cpp(i_b, match_indic, ms_b, n_reservoir)
					
					seq_des_r = private$seq_des_obj_priv_int$thin_duplicate()
					seq_des_r$.__enclos_env__$private$y = y[i_b]
					seq_des_r$.__enclos_env__$private$dead = dead[i_b]
					seq_des_r$.__enclos_env__$private$X = X[i_b, ]
					seq_des_r$.__enclos_env__$private$w = w[i_b]
					seq_des_r$.__enclos_env__$private$match_indic = match_indic_b
					#compute beta_T_hat					
					seq_inf_r = private$thin_duplicate()
					seq_inf_r$.__enclos_env__$private$X = seq_des_r$.__enclos_env__$private$X
					seq_inf_r$.__enclos_env__$private$setup_matching_data(seq_des_r, seq_des_r$.__enclos_env__$private$X)					
					beta_hat_T_bs[r] = seq_inf_r$compute_treatment_estimate()
				}
			} else {	
				cl = doParallel::makeCluster(private$num_cores)
				doParallel::registerDoParallel(cl)			
				#now copy them to each core's memory
				doParallel::clusterExport(cl, list("seq_des_obj", "n", "match_indic", "i_reservoir", "n_reservoir", "match_indic_b", "m", "y", "dead", "X", "w"), envir = environment())
				#now do the parallelization
				beta_hat_T_bs = doParallel::foreach(r = 1 : B, .inorder = FALSE, .combine = c) %dopar% {
					#draw a bootstrap of both the reservoir and matches - first create index vector					
					i_b = array(NA, n_reservoir + 2 * m)
					#draw a bootstrap sample of the reservoir
					i_reservoir_b_i = sample_int_replace_cpp(n_reservoir, n_reservoir)
					#load it into the index vector
					i_b[1 : n_reservoir] = i_reservoir[i_reservoir_b_i]
					#draw a bootstrap sample of the matches
					ms_b = sample_int_replace_cpp(m, m)
					#load it into the index vector 2 by 2
#					i_b_idx = n_reservoir
#					for (m0 in 1 : m){
#						i_b[(i_b_idx + 1) : (i_b_idx + 2)] = which(match_indic == ms_b[m0])
#						i_b_idx = i_b_idx + 2
#					}
					fill_i_b_with_matches_loop_cpp(i_b, match_indic, ms_b, n_reservoir)
					
					seq_des_r = private$seq_des_obj_priv_int$thin_duplicate()
					seq_des_r$.__enclos_env__$private$y = y[i_b]
					seq_des_r$.__enclos_env__$private$dead = dead[i_b]
					seq_des_r$.__enclos_env__$private$X = X[i_b, ]
					seq_des_r$.__enclos_env__$private$w = w[i_b]
					seq_des_r$.__enclos_env__$private$match_indic = match_indic_b
					#compute beta_T_hat					
					seq_inf_r = private$thin_duplicate()
					seq_inf_r$.__enclos_env__$private$seq_des_obj_priv_int = seq_des_r$.__enclos_env__$private
					seq_inf_r$.__enclos_env__$private$X = seq_des_r$.__enclos_env__$private$X
					seq_inf_r$.__enclos_env__$private$setup_matching_data(seq_des_r, seq_des_r$.__enclos_env__$private$X)
					seq_inf_r$compute_treatment_estimate()			
				}
				doParallel::stopCluster(cl)
				rm(cl); gc()
			}
			
			beta_hat_T_bs
		}
	),
	private = list(
		helper = NULL,
		KKstats = NULL,
		
		setup_matching_data = function(d, X){
			private$helper = SeqDesignInferenceHelperKK$new(d, X)
			private$KKstats = private$helper$get_post_matching_data()	
		},
		
		thin_duplicate = function(){
			i = super$thin_duplicate()
			h = SeqDesignInferenceHelperKK$new(seq_des_obj, private$get_X()) 
			i$.__enclos_env__$private$KKstats = h$get_post_matching_data()
			i
		},
		
		compute_model_matrix_with_matching_dummies = function(){
			private$helper$data_frame_with_matching_dummies()
		}	
	)
)
