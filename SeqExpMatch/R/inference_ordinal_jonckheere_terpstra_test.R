#' Jonckheere-Terpstra (JT) Test for Ordinal Responses
#'
#' @description
#' Exact Jonckheere-Terpstra (JT) rank test for a two-arm ordered alternative with an
#' ordinal response. For treatment versus control, the test statistic is the
#' sum of Mann-Whitney U counts across groups. This class provides the exact
#' distribution-based p-value.
#'
#' @inherit InferenceRand methods
#' @inherit InferenceBoot methods
#' @inherit InferenceAsymp methods
#' @inherit InferenceRandCI methods
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneBernoulli$new(n = nrow(x_dat), response_type = "ordinal", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalJonckheereTerpstraTest$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
InferenceOrdinalJonckheereTerpstraTest = R6::R6Class(
	"InferenceOrdinalJonckheereTerpstraTest",
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize the JT test object.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object.
		#' @param num_cores Number of CPU cores.
		#' @param verbose Whether to print progress.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			super$initialize(des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Returns the estimated treatment effect (JT superiority measure).
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Returns the exact two-sided p-value.
		compute_exact_two_sided_pval_for_treatment_effect = function(){
			private$shared()
			private$cached_values$p_exact
		}
	),

	private = list(
		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			
			res = exact_jonckheere_terpstra_pval_cpp(as.integer(private$y), as.integer(private$w))
			
			private$cached_values$superiority = res$superiority
			private$cached_values$beta_hat_T = res$superiority - 0.5
			private$cached_values$p_exact = res$p_exact
			private$cached_values$p_lower = res$p_lower
			private$cached_values$p_upper = res$p_upper
		}
	)
)
