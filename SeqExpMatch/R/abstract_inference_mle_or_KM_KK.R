#' Inference based on Maximum Likelihood for KK designs  
#'
#' @description
#' An abstract class
#' 
#'
SeqDesignInferenceMLEorKMKK = R6::R6Class("SeqDesignInferenceMLEorKMKK",
	inherit = SeqDesignInference,
	public = list(
		
        #' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference 
		#' 							(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 							for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		#'
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertClass(seq_des_obj, "SeqDesignKK14")
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
			
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}		

			n = private$seq_des_obj_priv_int$n	
			y = private$seq_des_obj_priv_int$y
			dead = private$seq_des_obj_priv_int$dead
			w = private$seq_des_obj_priv_int$w
			X = private$get_X()
			match_indic = private$seq_des_obj_priv_int$match_indic
			i_reservoir = which(match_indic == 0)
			n_reservoir = length(i_reservoir)
			m = private$cached_values$KKstats$m
			match_indic_b = c(rep(0, n_reservoir), rep(1 : m, each = 2))
			
			seq_des_r = private$seq_des_obj_priv_int$duplicate()
			seq_inf_r = private$duplicate()
			seq_inf_r$.__enclos_env__$private$X = seq_des_r$.__enclos_env__$private$X
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
					
					seq_des_r$.__enclos_env__$private$y = y[i_b]
					seq_des_r$.__enclos_env__$private$dead = dead[i_b]
					seq_des_r$.__enclos_env__$private$X = X[i_b, ]
					seq_des_r$.__enclos_env__$private$w = w[i_b]
					seq_des_r$.__enclos_env__$private$match_indic = match_indic_b
					#compute beta_T_hat					
					seq_inf_r$.__enclos_env__$private$seq_des_obj_priv_int = seq_des_r$.__enclos_env__$private
					seq_inf_r$.__enclos_env__$private$cached_values = list() #ensure nothing is kept between iterations		
					beta_hat_T_bs[r] = seq_inf_r$compute_treatment_estimate()
#					if (is.nan(beta_hat_T_bs[r])){
#						stop("boom")
#					}
				}
			} else {	
				cl = doParallel::makeCluster(private$num_cores)
				doParallel::registerDoParallel(cl)			
				#now copy them to each core's memory
				doParallel::clusterExport(cl, 
					list("seq_des_r", "seq_inf_r", "match_indic", "i_reservoir", "n_reservoir", "match_indic_b", "m", "y", "dead", "X", "w"), 
					envir = environment()
				)
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
					
					seq_des_r$.__enclos_env__$private$y = y[i_b]
					seq_des_r$.__enclos_env__$private$dead = dead[i_b]
					seq_des_r$.__enclos_env__$private$X = X[i_b, ]
					seq_des_r$.__enclos_env__$private$w = w[i_b]
					seq_des_r$.__enclos_env__$private$match_indic = match_indic_b
					#compute beta_T_hat
					seq_inf_r$.__enclos_env__$private$seq_des_obj_priv_int = seq_des_r$.__enclos_env__$private
					seq_inf_r$.__enclos_env__$private$cached_values = list() #ensure nothing is kept between iterations	
					seq_inf_r$compute_treatment_estimate()			
				}
				doParallel::stopCluster(cl)
				rm(cl); gc()
			}
			
			beta_hat_T_bs
		}
	),
	private = list(		
		compute_basic_match_data = function(){
			if (is.null(private$X)){
				private$X = private$get_X()
			}
			#cache data for speed	
			match_indic = private$seq_des_obj_priv_int$match_indic			
			m = max(match_indic)
			y = private$seq_des_obj_priv_int$y
			w = private$seq_des_obj_priv_int$w
			
			yTs_matched = array(NA, m)
			yCs_matched = array(NA, m)
			y_matched_diffs = array(NA, m)
			X_matched_diffs = matrix(NA, nrow = m, ncol = ncol(private$X))
			if (m > 0){
#				for (match_id in 1 : m){ #we want to just calculate the diffs inside matches and ignore the reservoir
#					yTs_matched[match_id] = y[w == 1 & match_indic == match_id]
#					yCs_matched[match_id] = y[w == 0 & match_indic == match_id]
#					
#					xmTvec = private$X[w == 1 & match_indic == match_id, ]
#					xmCvec = private$X[w == 0 & match_indic == match_id, ]
#					X_matched_diffs[match_id, ] = xmTvec - xmCvec
#				}
				match_data = match_diffs_cpp(w, match_indic, y, private$X, m)
				yTs_matched = match_data$yTs_matched
				yCs_matched = match_data$yCs_matched
				X_matched_diffs = match_data$X_matched_diffs
				y_matched_diffs = yTs_matched - yCs_matched
			}
			w_reservoir = w[match_indic == 0]

			private$cached_values$KKstats = list(
				X_matched_diffs = X_matched_diffs,
				yTs_matched = yTs_matched,
				yCs_matched = yCs_matched,
				y_matched_diffs = y_matched_diffs,
				X_reservoir = private$X[match_indic == 0, ],
				y_reservoir = y[match_indic == 0],
				w_reservoir = w_reservoir,
				nRT = sum(w_reservoir), #how many treatment observations are there in the reservoir?
				nRC = sum(w_reservoir == 0), #how many control observations are there in the reservoir?
				m = m
			)
		},
		
		compute_reservoir_and_match_statistics = function(){	
			nRC = private$cached_values$KKstats$nRC
			nRT = private$cached_values$KKstats$nRT		
			nR = nRT + nRC #how many observations are there in the reservoir?
					
			y_reservoir_T = private$cached_values$KKstats$y_reservoir[private$cached_values$KKstats$w_reservoir == 1] #get the reservoir responses from the treatment
			y_reservoir_C = private$cached_values$KKstats$y_reservoir[private$cached_values$KKstats$w_reservoir == 0] #get the reservoir responses from the control
			
			ssqD_bar = var_cpp(private$cached_values$KKstats$y_matched_diffs) / private$cached_values$KKstats$m
			ssqR = (var_cpp(y_reservoir_T) * (nRT - 1) + var_cpp(y_reservoir_C) * (nRC - 1)) / 
						(nR - 2) * (1 / nRT + 1 / nRC)
			private$cached_values$KKstats$d_bar = mean_cpp(private$cached_values$KKstats$y_matched_diffs)			
			private$cached_values$KKstats$ssqD_bar = ssqD_bar
			private$cached_values$KKstats$r_bar = mean_cpp(y_reservoir_T) - mean_cpp(y_reservoir_C) #compute the classic estimator from the reservoir: ybar_T - ybar_C
			private$cached_values$KKstats$ssqR = ssqR
			private$cached_values$KKstats$w_star = ssqR / (ssqR + ssqD_bar)
		},
		
		#not used now, but could be used for random effects models in the future
		compute_model_matrix_with_matching_dummies = function(){
			if (is.null(private$cached_values$data_frame_with_matching_dummies)){
				if (!is.null(private$seq_des_obj_priv_int$match_indic) & uniqueN(private$seq_des_obj_priv_int$match_indic) > 1){
					mm = model.matrix(~ 0 + factor(private$seq_des_obj_priv_int$match_indic)) 
					mm = mm[, 2 : (ncol(mm) - 1)]
				} else {
					mm = NULL
				}
				private$cached_values$data_frame_with_matching_dummies = cbind(data.frame(w = private$seq_des_obj_priv_int$w), mm)
			}
			private$cached_values$data_frame_with_matching_dummies
		}
	)
)
