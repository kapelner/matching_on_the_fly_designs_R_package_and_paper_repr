#' A class that provides for relevant methods when the designs are KK matching-on-the-fly
#'
#' @description
#' An abstract class
SeqDesignInferenceKKPassThrough = R6::R6Class("SeqDesignInferenceKKPassThrough",
	inherit = SeqDesignInference,
	public = list(
		
		#' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 								(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 								for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
			
			#there is no situation where we don't need the basic match data, so hit it right away
			private$match_indic = seq_des_obj$.__enclos_env__$private$match_indic
			private$compute_basic_match_data()
		},
		
			
		#' @description
		#' Creates the boostrap distribution of the estimate for the treatment effect
		#' 
		#' @param B						Number of bootstrap samples. The default is 501.
		#' 
		#' @return 	A vector of length \code{B} with the bootstrap values of the estimates of the treatment effect
		#' 
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignKK14$new(n = 6, response_type = "continuous")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceContinMultOLSKK$new(seq_des)
		#' beta_hat_T_bs = seq_des_inf$approximate_bootstrap_distribution_beta_hat_T(B = 5)
		#' beta_hat_T_bs
		#' }
		#' 			
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501){
			if (!private$is_KK){
				super$approximate_bootstrap_distribution_beta_hat_T(B)
			} else {
				assertCount(B, positive = TRUE)	
	
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
				indices = bootstrap_match_indices_cpp(match_indic, i_reservoir, n_reservoir, m, B)

				# Create lightweight function that computes estimate from KK stats
				# Avoids R6 object mutation overhead
				# Note: Each thread will create its own copy via duplicate_inference_fn

				duplicate_inference_fn <- function(){
					private$duplicate()
				}

				compute_estimate_from_kk_stats <- function(kk_stats){
					# Use the thread-local inference object passed in kk_stats
					thread_inf_obj = kk_stats$inf_obj

					# Store stats in cached_values (without inf_obj)
					kk_stats$inf_obj = NULL  # Remove to avoid storing SEXP
					thread_inf_obj$.__enclos_env__$private$cached_values$KKstats = kk_stats

					# Compute reservoir/match statistics if needed (for compound classes)
					if ("compute_reservoir_and_match_statistics" %in% names(thread_inf_obj$.__enclos_env__$private)){
						thread_inf_obj$.__enclos_env__$private$compute_reservoir_and_match_statistics()
					}

					# Compute and return the treatment estimate
					thread_inf_obj$compute_treatment_estimate()
				}

				# Use C++ loop with OpenMP parallelization
				# Each thread creates its own inference object copy for thread safety
				beta_hat_T_bs = kk_bootstrap_loop_cpp(
					indices,
					y,
					w,
					X,
					match_indic_b,
					m,
					duplicate_inference_fn,
					compute_estimate_from_kk_stats,
					private$num_cores
				)
				
				beta_hat_T_bs
			}			
		}
	),
	private = list(		
		match_indic = NULL,
		
		compute_basic_match_data = function(){
			if (is.null(private$X)){
				private$X = private$get_X()
			}
			#cache data for speed	
			match_indic = private$seq_des_obj_priv_int$match_indic
			if (is.null(match_indic)){
				match_indic = rep(0, private$n)
			}
			match_indic[is.na(match_indic)] = 0
			m = max(match_indic, na.rm = TRUE)
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
				X_reservoir = private$X[match_indic == 0, , drop = FALSE],
				y_reservoir = y[match_indic == 0],
				w_reservoir = w_reservoir,
				nRT = sum(w_reservoir, na.rm = TRUE), #how many treatment observations are there in the reservoir?
				nRC = sum(w_reservoir == 0, na.rm = TRUE), #how many control observations are there in the reservoir?
				m = m
			)
		},
		
		compute_concordant_and_discordant_match_statistics = function(){
			m = private$cached_values$KKstats$m
			y_matched_diffs = private$cached_values$KKstats$y_matched_diffs
			i_m_conc = which(y_matched_diffs == 0)
			i_m_disc = setdiff(seq_len(m), i_m_conc)
			private$cached_values$KKstats$i_m_conc = i_m_conc
			private$cached_values$KKstats$i_m_disc = i_m_disc
			private$cached_values$KKstats$n_m_conc = length(i_m_conc)
			private$cached_values$KKstats$n_m_disc = length(i_m_disc)
			private$cached_values$KKstats$y_matched_diffs_disc = y_matched_diffs[i_m_disc]
			private$cached_values$KKstats$X_matched_diffs_disc =
				private$cached_values$KKstats$X_matched_diffs[i_m_disc, drop = FALSE]
			i_conc = which(private$seq_des_obj_priv_int$match_indic %in% i_m_conc)
			private$cached_values$KKstats$X_conc = private$X[i_conc, drop = FALSE]
			private$cached_values$KKstats$y_conc = private$seq_des_obj_priv_int$y[i_conc]
			private$cached_values$KKstats$w_conc = private$seq_des_obj_priv_int$w[i_conc]
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
