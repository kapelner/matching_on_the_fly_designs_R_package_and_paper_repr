#' An incomplete / balanaced completely randomized Sequential Design
#' 
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data intialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#' for the  with appropriate permuted blocks based on \code{prob_T}
#' 					(e.g., if \code{prob_T = 2}, then this design would enforce n/2 T's and n/2 C's),
#' 
#' @export
SeqDesigniBCRD = R6::R6Class("SeqDesigniBCRD",
	inherit = SeqDesign,
	public = list(
		#' 				
		#' @description
		#' Initialize a balanced completely randomized sequential experimental design
		#'
		#' @param response_type 	The data type of response values which must be one of the following: 
		#' 							"continuous"(the default), 
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
		#' @return 			A new `SeqDesigniBCRD` object
		#' 
		#' @examples
		#' seq_des = SeqDesigniBCRD$new(response_type = "continuous")
		#'  
		initialize = function(
						response_type = "continuous",  
						prob_T = 0.5, 
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL, 
						verbose = FALSE
					) {					
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
			self$assert_fixed_sample()
			if (!all.equal(n * prob_T, as.integer(n * prob_T), check.attributes = FALSE)){
				stop("Design iBCRD requires that the fraction of treatments of the total sample size must be a natural number.")
			}
		}
	),
	private = list(		
		assign_wt = function(){
			n_T_total = round(private$n * private$prob_T) #this quantity should never be a fraction anyway as it was checked during initialization
			nT = sum(private$w == 1, na.rm = TRUE)
			nC = sum(private$w == 0, na.rm = TRUE)
			sample(c(rep(1, n_T_total - nT), rep(0, n_T_total - nC)), 1)
		},
		
		redraw_w_according_to_design = function(){
			n_T_total = round(private$t * private$prob_T) #this quantity should never be a fraction anyway as it was checked during initialization
			private$w = shuffle_cpp(c(rep(1, n_T_total), rep(0, private$t - n_T_total)))
		}
	)
)