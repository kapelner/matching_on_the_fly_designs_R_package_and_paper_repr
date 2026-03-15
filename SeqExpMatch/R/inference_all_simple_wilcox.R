#' Wilcoxon Rank-Sum Inference for Any Sequential Design
#'
#' @description
#' Non-parametric inference based on the two-sample Wilcoxon rank-sum test
#' (Mann-Whitney U) for any sequential experimental design (CRD, Efron, iBCRD,
#' or KK). The treatment effect estimate is the Hodges-Lehmann pseudo-median of
#' all pairwise treatment-minus-control response differences. Inference uses the
#' asymptotic normal approximation for the HL estimator.
#'
#' This class supports \code{continuous}, \code{count}, \code{proportion}, and
#' uncensored \code{survival} response types. It throws an informative error for
#' \code{incidence} (binary) responses and for censored survival data.
#'
#' The HL estimator satisfies the exact linearity property
#' \eqn{\hat{\Delta}(\mathbf{y}_T - \delta) = \hat{\Delta}(\mathbf{y}_T) - \delta},
#' so the randomization confidence interval inversion is exact (not approximate)
#' on the original response scale for continuous responses.
#'
#' @export
SeqDesignInferenceAllSimpleWilcox = R6::R6Class("SeqDesignInferenceAllSimpleWilcox",
	inherit = SeqDesignInference,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param seq_des_obj  A completed \code{SeqDesign} object.
		#' @param num_cores    Number of CPU cores for parallel bootstrap/randomization. Default 1.
		#' @param verbose      Whether to print progress messages. Default \code{FALSE}.
		#' @examples
		#' set.seed(1)
		#' x_dat <- data.frame(
		#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
		#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
		#' )
		#' seq_des <- SeqDesignCRD$new(n = nrow(x_dat), response_type = "continuous", verbose = FALSE)
		#' for (i in seq_len(nrow(x_dat))) {
		#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
		#' }
		#' seq_des$add_all_subject_responses(c(1.2, 0.9, 1.5, 1.8, 2.1, 1.7, 2.6, 2.2))
		#' infer <- SeqDesignInferenceAllSimpleWilcox$new(seq_des, verbose = FALSE)
		#' infer
		#'
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			res_type = seq_des_obj$get_response_type()
			if (res_type == "incidence"){
				stop(
					"Wilcoxon rank-sum inference is not implemented for incidence (binary) ",
					"responses: the Hodges-Lehmann estimator is degenerate (almost always 0) ",
					"on 0/1 data. Use SeqDesignInferenceAllSimpleMeanDiff or a clogit estimator instead."
				)
			}
			assertResponseType(res_type, c("continuous", "count", "proportion", "survival", "ordinal"))
			super$initialize(seq_des_obj, num_cores, verbose)
			if (private$any_censoring){
				stop(
					"Wilcoxon rank-sum inference does not support censored survival data. ",
					"Use SeqDesignInferenceSurvivalGehanWilcox for censored survival outcomes."
				)
			}
		},

		#' @description
		#' Returns the Hodges-Lehmann pseudo-median of all pairwise treatment-minus-control
		#' differences.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a \eqn{1 - \alpha} asymptotic confidence interval based on the normal
		#' approximation for the HL estimator.
		#' @param alpha Significance level; default 0.05 gives a 95\% CI.
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Returns a two-sided p-value for \eqn{H_0: \Delta = \delta} using the asymptotic
		#' normal approximation for the HL estimator.
		#' @param delta Null value; default 0.
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			yT = private$y[private$w == 1]
			yC = private$y[private$w == 0]

			if (length(yT) == 0L || length(yC) == 0L){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			mod = tryCatch(
				stats::wilcox.test(yT, yC, conf.int = TRUE),
				error = function(e) NULL
			)

			if (is.null(mod)){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			beta = as.numeric(mod$estimate)
			ci   = mod$conf.int  # 95% CI by default (conf.level = 0.95)
			se   = if (length(ci) == 2L) (ci[2] - ci[1]) / (2 * 1.96) else NA_real_

			private$cached_values$beta_hat_T   = if (length(beta) == 1L && is.finite(beta)) beta else NA_real_
			private$cached_values$s_beta_hat_T = if (length(se)   == 1L && is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z         = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop("Wilcoxon rank-sum: could not compute a finite standard error.")
			}
		}
	)
)
