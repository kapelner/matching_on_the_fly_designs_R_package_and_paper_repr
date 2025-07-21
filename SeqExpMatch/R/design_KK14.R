#' A Sequential Design
#' 
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data intialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#' 
#' @export
SeqDesignKK14 = R6::R6Class("SeqDesignKK14",
	inherit = SeqDesign,
	public = list(
		#' 				
		#' @description
		#' Initialize a matching-on-the-fly sequential experimental design which matches based on Kapelner and Krieger (2014) or Morrison and Owen (2025)
		#'
		#' @param response_type 	The data type of response values which must be one of the following: 
		#' 							"continuous", 
		#' 							"incidence", 
		#' 							"proportion", 
		#' 							"count", 
		#' 							"survival".
		#' 							This package will enforce that all added responses via the \code{add_subject_response} method will be
		#' 							of the appropriate type.
		#' @param prob_T	The probability of the treatment assignment. This defaults to \code{0.5}.
		#' @param include_is_missing_as_a_new_feature	If missing data is present in a variable, should we include another dummy variable for its
		#' 												missingness in addition to imputing its value? If the feature is type factor, instead of creating
		#' 												a new column, we allow missingness to be its own level. The default is \code{TRUE}.
		#' @param n			The sample size (if fixed). Default is \code{NULL} for not fixed.
		#' @param verbose	A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}.
		#' @param lambda   The quantile cutoff of the subject distance distribution for determining matches. If unspecified and \code{morrison = FALSE}, default is 10\%.
		#' @param t_0_pct  The percentage of total sample size n where matching begins. If unspecified and \code{morrison = FALSE}, default is 35\%.
		#' @param morrison 	Default is \code{FALSE} which implies matching via the KK14 algorithm using \code{lambda} and \code{t_0_pct} matching.
		#'					If \code{TRUE}, we use Morrison and Owen (2025)'s formula for \code{lambda} which differs in the fixed n versus variable n
		#'					settings and matching begins immediately with no wait for a certain reservoir size like in KK14.
		#' @param p			The number of covariate features. Must be specified when \code{morrison = TRUE} otherwise do not specify this argument.
  		#' @return 			A new `SeqDesignKK14` object 
		#' 
		#' @examples
		#' seq_des = SeqDesignKK14$new(response_type = "continuous")
		#'  
		initialize = function(
			response_type, 
			prob_T = 0.5,
			include_is_missing_as_a_new_feature = TRUE, 
			verbose = FALSE,
			n = NULL,
			lambda = NULL,
			t_0_pct = NULL,
			morrison = FALSE,
			p = NULL
		){
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, verbose, n)
			self$assert_even_allocation()
			private$assert_KK_and_morrison_parameters_corrrect(lambda, t_0_pct, morrison, p)
			private$morrison = morrison
			
			if (morrison){
				private$p = p				
				private$compute_lambda = 	if (super$is_fixed_sample()){
												function(){private$n^(-1 / (2 * private$p))}
											} else {
												function(){private$t^(-1 / (2 * private$p))}
											}				
			} else {
				private$match_indic = array(NA, n)		
				private$lambda = ifelse(is.null(lambda), 0.1, lambda) #10% is the default	
				private$compute_lambda = function(){private$lambda}
				private$t_0 = round(ifelse(is.null(t_0_pct), 0.35, t_0_pct) * n) #35% is the default			
			}				
		},
		
		#' @description
		#' For KK designs only, this returns a list with useful matching statistics.
		#' 
		#' @return 	A list with the following data: \code{num_matches}, \code{prop_subjects_matched}, 
		#' 			\code{num_subjects_remaining_in_reservoir}, \code{prop_subjects_remaining_in_reservoir}.
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "KK14", response_type = "continuous")
		#' 
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' 
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des$matching_statistics()
		#'
		matching_statistics = function(){
			if (private$t == 0){
				stop("The experiment has not begun yet")
			}
			num_subjects_matched = sum(private$match_indic != 0, na.rm = TRUE)
			num_subjects_remaining_in_reservoir = private$t - num_subjects_matched
			list(
				num_matches = length(unique(private$match_indic[private$match_indic != 0])) ,
				prop_subjects_matched = num_subjects_matched / private$t,
				num_subjects_remaining_in_reservoir = num_subjects_remaining_in_reservoir,
				prop_subjects_remaining_in_reservoir = num_subjects_remaining_in_reservoir / private$t
			)					
		}
	),
	private = list(
		uses_covariates = TRUE,
		morrison = NULL,
		t_0 = NULL,	
		lambda = NULL,
		p = NULL,
		compute_lambda = NULL,	
		match_indic = NA, #works for Morrison and non-Morrison
		
		duplicate = function(){
			d = super$duplicate()
			d$.__enclos_env__$private$morrison = private$morrison
			d$.__enclos_env__$private$t_0 = private$t_0
			d$.__enclos_env__$private$lambda = private$lambda
			d$.__enclos_env__$private$p = private$p
			d$.__enclos_env__$private$compute_lambda = private$compute_lambda
			d$.__enclos_env__$private$match_indic = private$match_indic
			d
		},
	
		assert_KK_and_morrison_parameters_corrrect = function(lambda, t_0_pct, morrison, p){
			if (morrison){
				assert_count(p)				
			} else {
				self$assert_fixed_sample()
				if (!is.null(lambda)){
					assertNumeric(lambda, lower = 0, upper = 1)
				}
				if (!is.null(t_0_pct)){
					assertNumeric(t_0_pct, lower = .Machine$double.eps, upper = 1)
				}						
			}
		},	
		
		too_early_to_match = function(){
			sum(private$match_indic == 0, na.rm = TRUE) == 0 | private$t <= (ncol(private$Xraw) + 2) | 
							(!private$morrison & private$t <= private$t_0)
		},
		
		assign_wt = function(){
			wt = 	if (private$too_early_to_match()){
						#we're early, so randomize
						private$match_indic[private$t] = 0 #zero means "reservoir", >0 means match number 
						private$assign_wt_CRD()
					} else {
						all_subject_data = private$compute_all_subject_data()
						# cat("else\n")
						#first calculate the threshold we're operating at	
						#when inverting, ensure full rank by adding eps * I			
						S_xs_inv = solve(var(all_subject_data$X_prev) + diag(.Machine$double.eps, all_subject_data$rank_prev), tol = .Machine$double.xmin)
						F_crit =  qf(private$compute_lambda(), all_subject_data$rank_prev, private$t - all_subject_data$rank_prev)
						n = self$get_n()
						T_cutoff_sq = all_subject_data$rank_prev * (n - 1) / (n - all_subject_data$rank_prev) * F_crit
						#now iterate over all items in reservoir and take the minimum distance x
						reservoir_indices = which(private$match_indic == 0)
						sqd_distances_times_two = array(NA, length(reservoir_indices))
						for (r in 1 : length(reservoir_indices)){
							x_r_x_new_delta = all_subject_data$xt_prev - all_subject_data$X_prev[reservoir_indices[r], ]
							sqd_distances_times_two[r] = t(x_r_x_new_delta) %*% S_xs_inv %*% x_r_x_new_delta		
						}					
						#find minimum distance index
						min_sqd_dist_index = which(sqd_distances_times_two == min(sqd_distances_times_two))
						if (length(sqd_distances_times_two[min_sqd_dist_index]) > 1 || length(T_cutoff_sq) > 1){
							min_sqd_dist_index = min_sqd_dist_index[1] #if there's a tie, just take the first one
						}
						#if it's smaller than the threshold, we're in business: match it
						if (sqd_distances_times_two[min_sqd_dist_index] < T_cutoff_sq){
							match_num = max(private$match_indic, na.rm = TRUE) + 1
							private$match_indic[reservoir_indices[min_sqd_dist_index]] = match_num
							private$match_indic[private$t] = match_num
							#assign opposite
							1 - private$w[reservoir_indices[min_sqd_dist_index]]
						} else { #otherwise, randomize and add it to the reservoir
							private$match_indic[private$t] = 0
							private$assign_wt_CRD()
						}
					}
			if (is.na(private$match_indic[private$t])){
				stop("no match data recorded")
			}
			wt
		},
		
		redraw_w_according_to_design = function(){
			#we rearrange within each match set (and the reservoir which is when m = 0)
			for (m in 0 : max(private$match_indic)){
				private$w[private$match_indic == m] = shuffle_cpp(private$w[private$match_indic == m]) 
			}
		}	
	)
)