#' Exact incidence inference for KK designs via Zhang's combined matched-pair / reservoir tests
#'
#' @description
#' Extends \code{SeqDesignInferenceIncidExactZhang} with the KK matched-pair exact test,
#' producing Zhang's combined exact p-values and confidence intervals for KK
#' matching-on-the-fly designs with incidence outcomes.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignKK14$new(n = nrow(x_dat), response_type = "incidence", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(c(0, 1, 0, 1, 0, 1, 1, 0))
#' infer <- SeqDesignInferenceIncidKKExactZhang$new(seq_des, verbose = FALSE)
#' infer
#'
SeqDesignInferenceIncidKKExactZhang = R6::R6Class("SeqDesignInferenceIncidKKExactZhang",
	inherit = SeqDesignInferenceIncidExactZhang,
	public = list(

		#' @description
		#' Initialize the exact incidence inference object for a KK design.
		#' @param	seq_des_obj		A SeqDesignKK14 object (or subclass) with an incidence response.
		#' @param	num_cores			Number of CPU cores for parallel processing.
		#' @param	verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),
	private = list(

		assert_supported_design = function(seq_des_obj){
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
		},

		# Exact McNemar test under H0: OR_pairs = exp(delta_0).
		# Given k = d_plus + d_minus discordant pairs,
		# d_plus | k ~ Binomial(k, expit(delta_0)) under H0.
		# Concordant pairs contribute no information and are ignored.
		compute_exact_pval_matched_pairs = function(delta_0){
			exact_stats = private$get_exact_zhang_stats()
			if (exact_stats$m == 0L) return(NA_real_)
			if (exact_stats$d_plus + exact_stats$d_minus == 0L) return(NA_real_)

			zhang_exact_binom_pval_cpp(exact_stats$d_plus, exact_stats$d_minus, delta_0)
		}
	)
)
