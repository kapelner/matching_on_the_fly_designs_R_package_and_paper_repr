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
		#' @param p			The number of covariate measurements. Must be specified.
		#' @param verbose	A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}.
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
			n = NULL,
			p
		){
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, verbose, n)
			assert_count(p)
			private$p = p
			
			if (super$is_fixed_sample()){
				private$compute_lambda = function(){private$n^(-1 / (2 * private$p))}
			} else {
				private$compute_lambda = function(){private$t^(-1 / (2 * private$p))}
			}	
	),
	private = list(
		p = NULL,
				
		too_early_to_match = function(){
			sum(private$match_indic == 0, na.rm = TRUE) == 0 | private$t <= (ncol(private$Xraw) + 2)
		}
	)
)