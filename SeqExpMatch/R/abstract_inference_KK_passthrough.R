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

				n = private$n
				y = private$y
				dead = private$dead
				w = private$w
				X = private$get_X()
				match_indic = private$match_indic
				if (is.null(match_indic)){
					match_indic = rep(0, n)
				}
				match_indic[is.na(match_indic)] = 0
				m = private$cached_values$KKstats$m

				# Identify reservoir and matched observations
				i_reservoir = which(match_indic == 0)
				n_reservoir = length(i_reservoir)

				# Pure R KK bootstrap implementation
				beta_hat_T_bs = rep(NA_real_, B)

				for (b in 1:B) {
					# Bootstrap indices for KK design:
					# - Resample reservoir with replacement
					# - Resample entire matched pairs with replacement
					i_reservoir_b = sample(i_reservoir, n_reservoir, replace = TRUE)

				# For matched pairs, sample which pairs to include (with replacement)
				if (m > 0) {
					pairs_to_include = sample(1:m, m, replace = TRUE)
					i_matched_b = integer(0)
					match_indic_b_matched = integer(0)
					for (new_pair_id in 1:m) {
						# Get both members of the original pair
						original_pair_id = pairs_to_include[new_pair_id]
						pair_indices = which(match_indic == original_pair_id)
						i_matched_b = c(i_matched_b, pair_indices)
						# Assign new pair IDs sequentially
						match_indic_b_matched = c(match_indic_b_matched, new_pair_id, new_pair_id)
					}
				} else {
					i_matched_b = integer(0)
					match_indic_b_matched = integer(0)
				}

					# Combine reservoir and matched indices
					i_b = c(i_reservoir_b, i_matched_b)

					# Create bootstrap sample
					y_b = y[i_b]
					dead_b = dead[i_b]
					w_b = w[i_b]
					X_b = X[i_b, , drop = FALSE]
					match_indic_b = c(rep(0, n_reservoir), match_indic_b_matched)

					# Create duplicate inference object
					boot_inf_obj = self$duplicate()

					# Update with bootstrap data
					boot_inf_obj$.__enclos_env__$private$y = y_b
					boot_inf_obj$.__enclos_env__$private$dead = dead_b
					boot_inf_obj$.__enclos_env__$private$w = w_b
					boot_inf_obj$.__enclos_env__$private$X = X_b
					boot_inf_obj$.__enclos_env__$private$cached_values = list()

					# Compute KK match statistics for bootstrap sample
					tryCatch({
						# Call compute_basic_match_data with the new match_indic
						boot_inf_obj$.__enclos_env__$private$match_indic = match_indic_b
						boot_inf_obj$.__enclos_env__$private$compute_basic_match_data()

						# For compound classes, compute reservoir/match statistics
						if (private$object_has_private_method(boot_inf_obj, "compute_reservoir_and_match_statistics")){
							boot_inf_obj$.__enclos_env__$private$compute_reservoir_and_match_statistics()
						}

						# Compute treatment estimate
						beta_hat_T_bs[b] = boot_inf_obj$compute_treatment_estimate()
					}, error = function(e) {
						if (private$verbose) {
							cat("      KK Bootstrap sample", b, "failed:", e$message, "\n")
						}
						beta_hat_T_bs[b] = NA_real_
					})
				}

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
			match_indic = private$match_indic
			if (is.null(match_indic)){
				match_indic = rep(0, private$n)
			}
			match_indic[is.na(match_indic)] = 0
			m = max(match_indic, na.rm = TRUE)
			y = private$y
			w = private$w

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
			i_conc = which(private$match_indic %in% i_m_conc)
			private$cached_values$KKstats$X_conc = private$X[i_conc, drop = FALSE]
			private$cached_values$KKstats$y_conc = private$y[i_conc]
			private$cached_values$KKstats$w_conc = private$w[i_conc]
		},

		#not used now, but could be used for random effects models in the future
		compute_model_matrix_with_matching_dummies = function(){
			if (is.null(private$cached_values$data_frame_with_matching_dummies)){
				if (!is.null(private$match_indic) & uniqueN(private$match_indic) > 1){
					mm = model.matrix(~ 0 + factor(private$match_indic))
					mm = mm[, 2 : (ncol(mm) - 1)]
				} else {
					mm = NULL
				}
				private$cached_values$data_frame_with_matching_dummies = cbind(data.frame(w = private$w), mm)
			}
			private$cached_values$data_frame_with_matching_dummies
		}
	)
)
