#' Atkinson's (1982) Covariate-Adjusted Biased Coin Sequential Design
#' 
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data intialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#' 
#' @export
SeqDesignAtkinson = R6::R6Class("SeqDesignAtkinson",
	inherits = SeqDesign,
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
		#' @return 			A new `SeqDesign` object of the specific type
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(response_type = "continuous")
		#'  
		initialize = function(
						response_type, 
						prob_T = 0.5, 
						include_is_missing_as_a_new_feature = TRUE, 
						verbose = TRUE,
						n = NULL
					) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, verbose, n)
		}
	),
	private = list(
		assign_wt = function(){			
			#if it's too early in the trial or if all the assignments are the same, then randomize
			if (private$t <= ncol(private$Xraw) + 2 + 1 | length(unique(private$t)) == 1){
				private$assign_wt_CRD()
			} else {
				all_subject_data = private$compute_all_subject_data()
				#this matrix is [w | 1 | X]
				Xprev_with_w = cbind(private$w[1 : (private$t - 1)], 1, all_subject_data$X_prev)
				XwtXw = t(Xprev_with_w) %*% Xprev_with_w	
				tryCatch({										
					M = (private$t - 1) * solve(XwtXw, tol = .Machine$double.xmin)
					A = M[1, 2 : (all_subject_data$rank_prev + 2)] %*% c(1, all_subject_data$xt_prev) 
					s_over_A_plus_one_sq = (M[1, 1] / A + 1)^2
					#assign via the Atkinson weighted biased coin
					rbinom(1, 1, s_over_A_plus_one_sq / (s_over_A_plus_one_sq + 1))	
				}, error = function(e){ #sometimes XtX is still not invertible... 
					#so in that case... just randomize
					private$assign_wt_CRD()
				})
			}
		},
		
		redraw_w_according_to_design = function(){
			n = private$t
			private$t = 0 #reset
			for (t in 1 : n){
				private$t = private$t + 1
				private$w[t] = private$assign_wt_Atkinson()
			}
		}
	)
)