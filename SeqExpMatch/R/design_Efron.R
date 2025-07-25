#' Efron's (1971) Biased Coin Sequential Design
#' 
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data intialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#' 
#' @export
SeqDesignEfron = R6::R6Class("SeqDesignEfron",
	inherit = SeqDesign,
	public = list(
		#' 				
		#' @description
		#' Initialize an Efron's weighted coin sequential experimental design
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
		#' @param weighted_coin_prob The probability of Efron's weighted coin. Defualt is \code{NULL} for Efron's default of 2/3.
		#' @param thin		For internal use only. Do not specify. You can thank R6's single constructor-only for this coding noise.
		#'
		#' @return 			A new `SeqDesignEfron` object
		#' 
		#' @examples
		#' seq_des = SeqDesignEfron$new(response_type = "continuous")
		#'  
		initialize = function(
			response_type = "continuous",  
			prob_T = 0.5,
			include_is_missing_as_a_new_feature = TRUE, 
			n = NULL,
			verbose = FALSE,
			weighted_coin_prob = NULL,
			thin = FALSE
		){
			if (!thin){
				super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
				self$assert_fixed_sample()
	
				if (is.null(weighted_coin_prob)){
					weighted_coin_prob = 2 / 3 #default Efron coin
				} else {
					assertNumeric(weighted_coin_prob, lower = 0, upper = 1)
				}
				private$weighted_coin_prob = weighted_coin_prob					
			}
		}
		
	),
	private = list(
		uses_covariates = FALSE,
		weighted_coin_prob = NULL,
		
		duplicate = function(){
			d = super$duplicate()
			d$.__enclos_env__$private$weighted_coin_prob = private$weighted_coin_prob
			d
		},
		
		assign_wt = function(){
			n_T = sum(private$w, na.rm = TRUE)
			n_C = private$t - n_T
			if (n_T * private$prob_T > n_C * (1 - private$prob_T)){
				rbinom(1, 1, 1 - private$weighted_coin_prob)
			} else if (n_T * private$prob_T < n_C * (1 - private$prob_T)){
				rbinom(1, 1, private$weighted_coin_prob)
			} else {
				private$assign_wt_CRD()
			}				
		}
	)
)