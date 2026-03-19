#' Paired Sign Test Inference for KK Designs with Ordinal Response
#'
#' @description
#' Fits a paired sign test for ordinal responses under a KK matching-on-the-fly design.
#' For matched pairs, it considers the sign of the within-pair differences.
#' Reservoir subjects are not included in this simple paired test.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignKK14$new(n = nrow(x_dat), response_type = "ordinal", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- SeqDesignInferenceOrdinalPairedSignTest$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
SeqDesignInferenceOrdinalPairedSignTest = R6::R6Class("SeqDesignInferenceOrdinalPairedSignTest",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param	seq_des_obj		A SeqDesign object (must be a KK design).
		#' @param	num_cores			Number of CPU cores.
		#' @param	verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "ordinal")
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design.")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		#' @description
		#' Returns the estimated treatment effect (proportion of pairs where T > C).
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the confidence interval for the probability P(T > C).
		#' @param	alpha					The significance level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the p-value for the sign test.
		#' @param	delta					The null difference (must be 0 for sign test).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			if (delta != 0) stop("Sign test only supports testing against delta = 0.")
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			
			diffs = private$cached_values$KKstats$y_matched_diffs
			# Sign test on matched pairs: ignore ties (diff == 0)
			pos = sum(diffs > 0)
			neg = sum(diffs < 0)
			n_eff = pos + neg
			
			if (n_eff == 0){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}
			
			# Estimate: proportion of non-tied pairs favoring treatment
			p_hat = pos / n_eff
			# Beta is usually defined as p_hat - 0.5 for centered tests
			private$cached_values$beta_hat_T = p_hat - 0.5
			# Standard error for proportion
			se = sqrt(p_hat * (1 - p_hat) / n_eff)
			
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z = TRUE
		}
	)
)
