#' A completely randomized / Bernoulli Sequential Design
#' 
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data initialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#' 
#' @export
SeqDesignCRD = R6::R6Class("SeqDesignCRD",
	inherit = SeqDesign,
	public = list(
		#' 				
		#' @description
		#' Initialize a completely randomized sequential experimental design
		#'
		#' @param response_type 	The data type of response values which must be one of the following:
		#' 								"continuous" (the default),
		#' 								"incidence",
		#' 								"proportion",
		#' 								"count",
		#' 								"survival".
		#' 								This package will enforce that all added responses via the \code{add_subject_response} method will be
		#' 								of the appropriate type.
		#' @param prob_T	The probability of the treatment assignment. This defaults to \code{0.5}.
		#' @param include_is_missing_as_a_new_feature	If missing data is present in a variable, should we include another dummy variable for its
		#' 								missingness in addition to imputing its value? If the feature is type factor, instead of creating
		#' 								a new column, we allow missingness to be its own level. The default is \code{TRUE}.
		#' @param n			The sample size (if fixed). Default is \code{NULL} for not fixed.
		#' @param verbose	A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}.
		#' @return 			A new `SeqDesignCRD` object
		#' 
		#' @examples
		#' seq_des = SeqDesignCRD$new(response_type = "continuous")
		#'  
		initialize = function(
						response_type = "continuous",  
						prob_T = 0.5, 
						include_is_missing_as_a_new_feature = TRUE, 
						n = NULL,
						verbose = FALSE
					) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
		}
	),
	private = list(				
		assign_wt = function(){
			private$assign_wt_CRD()
		},
		
		redraw_w_according_to_design = function(){
			private$w = rbinom(private$t, 1, private$prob_T)
		}
	)
)
