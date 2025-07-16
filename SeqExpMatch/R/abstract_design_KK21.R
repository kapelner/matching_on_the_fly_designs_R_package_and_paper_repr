#' A Sequential Design
#' 
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data intialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#' 
SeqDesignKK21abstract = R6::R6Class("SeqDesignKK21abstract",
	inherit = SeqDesignKK14,
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
			n = NULL,
			lambda = NULL,
			t_0_pct = NULL
		){
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, verbose, n, lambda, t_0_pct)	
		}
	),
	private = list(
		covariate_weights = NULL	
	)
)