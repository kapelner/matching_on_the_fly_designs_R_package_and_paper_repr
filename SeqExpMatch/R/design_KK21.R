#' Kapelner and Krieger's (2014) Covariate-Adjusted Matching on the Fly Sequential Design
#' 
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data intialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#' 
#' @export
SeqDesignKK21 = R6::R6Class("SeqDesignKK21",
	inherit = SeqDesignKK21abstract,
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
		#' @param n			The sample size (must be fixed).
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
			self$assert_fixed_sample()
			private$match_indic = array(NA, n)
			
			if (is.null(lambda)){
				lambda = 0.1 #default
			} else {
				assertNumeric(lambda, lower = 0, upper = 1)
			}			
			private$lambda = lambda	
			private$compute_lambda = function(){private$lambda}
			if (is.null(t_0_pct)){
				private$t_0 = round(0.35 * n) #default
			} else {
				assertNumeric(t_0_pct, lower = .Machine$double.eps, upper = 1)
				private$t_0 = round(t_0_pct * n)
			}	
		}
	),
	private = list(
		t_0 = NULL,	
		lambda = NULL,
		
		too_early_to_match = function(){
			sum(private$match_indic == 0, na.rm = TRUE) == 0 | private$t <= (ncol(private$Xraw) + 2) | private$t <= private$t_0
		}
	)
)