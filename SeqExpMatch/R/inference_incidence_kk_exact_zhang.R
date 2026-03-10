#' Exact incidence inference for KK designs via Zhang's combined matched-pair / reservoir tests
#'
#' @description
#' Extends \code{SeqDesignInferenceIncidExactZhang} with the KK matched-pair exact test,
#' producing Zhang's combined exact p-values and confidence intervals for KK
#' matching-on-the-fly designs with incidence outcomes.
#'
#' @export
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
		compute_rand_pval_matched_pairs = function(delta_0){
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats) || KKstats$m == 0L) return(NA_real_)

			yTs = KKstats$yTs_matched
			yCs = KKstats$yCs_matched
			d_plus  = sum(yTs == 1L & yCs == 0L)
			d_minus = sum(yTs == 0L & yCs == 1L)
			k = d_plus + d_minus
			if (k == 0L) return(NA_real_)

			p0 = exp(delta_0) / (1 + exp(delta_0))
			stats::binom.test(d_plus, k, p = p0, alternative = "two.sided")$p.value
		}
	)
)
