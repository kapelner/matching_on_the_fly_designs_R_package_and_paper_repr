#' Exact incidence inference via Zhang's combined reservoir test
#'
#' @description
#' Provides incidence-specific exact inference based on the Zhang (2026) combined
#' method. For CRD designs this uses the reservoir exact test only.
#'
#' @details
#' This class provides the \code{compute_exact_two_sided_pval_for_treatment_effect} and
#' \code{compute_exact_confidence_interval} methods. Standard methods like
#' \code{compute_treatment_estimate} and \code{compute_mle_confidence_interval}
#' are not implemented for this class.
#'
#' @export
SeqDesignInferenceIncidExactZhang = R6::R6Class("SeqDesignInferenceIncidExactZhang",
	inherit = SeqDesignInferenceIncidExactZhangAbstract,
	public = list(

		#' @description
		#' Initialize the exact incidence inference object.
		#' @param	seq_des_obj		A SeqDesign object with an incidence response.
		#' @param	num_cores			Number of CPU cores for parallel processing.
		#' @param	verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "incidence")
			private$assert_supported_design(seq_des_obj)
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		}
	),

	private = list(

		assert_supported_design = function(seq_des_obj){
			if (!is(seq_des_obj, "SeqDesignCRD")){
				stop(class(self)[1], " requires a completely randomized design (SeqDesignCRD). Use SeqDesignInferenceIncidKKExactZhang for KK matching-on-the-fly designs.")
			}
		}
	)
)
