#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring)
#' sequential experimental design estimation and test object
#' after the sequential design is completed.
#'
#'
#' @export
InferenceSurvivalMultiCoxPHRegr = R6::R6Class("InferenceSurvivalMultiCoxPHRegr",
	inherit = InferenceSurvivalUniCoxPHRegr,
	public = list(


		#' @description
		#' Computes the appropriate estimate
		#'
		#' @return	The setting-appropriate (see description) numeric estimate of the treatment effect
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = DesignSeqOneByOneBernoulli$new(n = 6, response_type = "survival")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(
		#'   ys = c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43),
		#'   deads = c(1L, 0L, 1L, 1L, 0L, 1L)
		#' )
		#'
		#' seq_des_inf = InferenceSurvivalMultiCoxPHRegr$new(seq_des)
		#   seq_des_inf$compute_treatment_estimate()
		# }
		#
		),

		private = list(
		generate_mod = function(estimate_only = FALSE){
			surv_obj = survival::Surv(private$y, private$dead)
			tryCatch({
				coxph_mod = suppressWarnings(survival::coxph(surv_obj ~ cbind(private$w, private$get_X())))

				if (estimate_only) {
					list(
						b = c(0, stats::coef(coxph_mod)),
						ssq_b_2 = NA_real_
					)
				} else {
					list(
						b = c(0, stats::coef(coxph_mod)),
						ssq_b_2 = as.numeric(stats::vcov(coxph_mod))[1]
					)
				}
			}, error = function(e){
				list(
					b = c(NA, NA),
					ssq_b_2 = NA
				)
			})
		}
		)
		)
