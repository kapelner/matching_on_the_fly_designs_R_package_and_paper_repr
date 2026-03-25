#' A Sequential Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data initialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#'
#' @export
SeqDesignKK14 = R6::R6Class("SeqDesignKK14",
	inherit = SeqDesign,
	public = list(
		#' @description
		#' Initialize a KK14 sequential experimental design
		#'
		#' @param	response_type 	"continuous", "incidence", "proportion", "count", "survival", or "ordinal".
		#' @param	prob_T	Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param num_cores The number of CPU cores to use.
		#' @param verbose A flag for verbosity.
		#' @param lambda The penalty parameter for covariate imbalance.
		#' @param t_0_pct The percentage of subjects to allocate before matching begins.
		#' @param morrison If TRUE, use Morrison's method for matching.
		#' @param p The number of covariates to use for matching.
		#'
		#' @return	A new `SeqDesignKK14` object
		initialize = function(
						response_type = "continuous",
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,
						num_cores = 1,
						verbose = FALSE,
						lambda = NULL,
						t_0_pct = NULL,
						morrison = FALSE,
						p = NULL
					) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			private$uses_covariates = TRUE
		},

		#' @description
		#' Assign the next subject to a treatment group
		#'
		#' @return 	The treatment assignment (0 or 1)
		assign_wt = function(){
			wt = 	if (private$too_early_to_match()){
						#we're early, so randomize
						private$m[private$t] = 0 #zero means "reservoir", >0 means match number
						self$assign_wt_Bernoulli()
					} else {
						all_subject_data = private$compute_all_subject_data()
						#first calculate the threshold we're operating at
						#when inverting, ensure full rank by adding eps * I
						inv_Xt_X_plus_lambda_I = solve(all_subject_data$XtX_all + private$lambda * diag(all_subject_data$p_all))
						#formula 14
						threshold = as.numeric(all_subject_data$xt_all %*% inv_Xt_X_plus_lambda_I %*% all_subject_data$xt_all)
						
						#now find the best match in the reservoir
						reservoir_indices = which(private$m[1 : (private$t - 1)] == 0)
						sqd_mahal_dists = compute_mahal_distances_cpp(all_subject_data$X_prev, all_subject_data$xt_prev, inv_Xt_X_plus_lambda_I)
						min_sqd_dist_index = which.min(sqd_mahal_dists)
						min_sqd_dist = sqd_mahal_dists[min_sqd_dist_index]
						
						if (min_sqd_dist < threshold){
							#we matched!
							match_num = max(private$m) + 1
							#update previous subject's match number
							private$m[reservoir_indices[min_sqd_dist_index]] = match_num
							#update current subject's match number
							private$m[private$t] = match_num
							#assign opposite
							1 - private$w[reservoir_indices[min_sqd_dist_index]]
						} else { #otherwise, randomize and add it to the reservoir
							private$m[private$t] = 0
							self$assign_wt_Bernoulli()
						}
					}
			if (is.na(private$m[private$t])){
				stop("no match data recorded")
			}
			wt
		},

		#' @description
		#' Draw multiple treatment assignment vectors according to KK14 design.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			generate_permutations_kk_cpp(
				as.integer(private$m),
				as.integer(r),
				as.numeric(private$prob_T)
			)$w_mat
		}
	),
	private = list(
		m = NULL,
		lambda = NULL,
		t_0_pct = NULL,
		morrison = NULL,
		p = NULL,

		too_early_to_match = function(){
			private$t <= private$t_0_pct * private$n
		},

		redraw_w_according_to_design = function(){
			private$w = redraw_w_kk14_cpp(private$m, private$w)
		}
	)
)
