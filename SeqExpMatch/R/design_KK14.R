#' Kapelner and Krieger's (2014) Covariate-Adjusted Matching on the Fly Sequential Design
#' 
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data intialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#' 
#' @export
SeqDesignKK14 = R6::R6Class("SeqDesignKK14",
	inherit = SeqDesignKK14abstract,
	public = list(
		#' 				
		#' @description
		#' Initialize a sequential experimental design
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
		#' @param lambda   The quantile cutoff of the subject distance distribution for determining matches. If unspecified, default is 10%.
		#' @param t_0_pct  The percentage of total sample size n where matching begins. If unspecified, default is 35%.
		#' @return 			A new `SeqDesign` object of the specific type
		#' 
		#' @examples
		#' seq_des = SeqDesignKK14$new(response_type = "continuous")
		#'  
		initialize = function(
			response_type, 
			prob_T = 0.5,
			include_is_missing_as_a_new_feature = TRUE, 
			verbose = TRUE,
			n,
			lambda = NULL,
			t_0_pct = NULL
		){
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, verbose, n)
			private$match_indic = array(NA, n)
			
			if (is.null(lambda)){
				lambda = 0.1 #default
			} else {
				assertNumeric(lambda, lower = 0, upper = 1)
			}			
			private$lambda = lambda	
			if (is.null(t_0_pct)){
				private$t_0 = round(0.35 * n) #default
			} else {
				assertNumeric(t_0_pct, lower = .Machine$double.eps, upper = 1)
				private$t_0 = round(t_0_pct * n)
			}	
		}		
	),
	private = list(
		lambda = NULL,
		t_0 = NULL,		
		
		assign_wt = function(){
			wt = 	if (sum(private$match_indic == 0, na.rm = TRUE) == 0 | private$t <= (ncol(private$Xraw) + 2) | private$t <= private$t_0){
						#we're early, so randomize
						private$match_indic[private$t] = 0 #zero means "reservoir", >0 means match number 
						private$assign_wt_CRD()
					} else {
						all_subject_data = private$compute_all_subject_data()
						# cat("else\n")
						#first calculate the threshold we're operating at	
						#when inverting, ensure full rank by adding eps * I			
						S_xs_inv = solve(var(all_subject_data$X_prev) + diag(.Machine$double.eps, all_subject_data$rank_prev), tol = .Machine$double.xmin)
						F_crit =  qf(private$other_params$lambda, all_subject_data$rank_prev, private$t - all_subject_data$rank_prev)
						T_cutoff_sq = all_subject_data$rank_prev * (private$n - 1) / (private$n - all_subject_data$rank_prev) * F_crit
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
		}		
	)
)