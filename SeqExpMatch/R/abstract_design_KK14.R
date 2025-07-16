#' A Sequential Design
#' 
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data intialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#' 
SeqDesignKK14abstract = R6::R6Class("SeqDesignKK14abstract",
	inherit = SeqDesign,
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
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, verbose, n)
			self$assert_even_allocation()
			self$assert_fixed_sample()		
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
		t_0 = NULL,		
		match_indic = NA #works for Morrison and non-Morrison		
	)
)