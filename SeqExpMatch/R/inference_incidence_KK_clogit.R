#' Conditional Logistic Regression Inference for KK Designs with Binary Response
#'
#' @description
#' Fits a conditional logistic regression model (via \code{bclogit::clogit}) using
#' the matched pairs from a KK matching-on-the-fly design. Only matched subjects
#' (i.e. those with a non-zero match indicator) are included in the model; reservoir
#' subjects are excluded because conditional logistic regression requires paired data.
#'
#' @export
SeqDesignInferenceIncidKKClogit = R6::R6Class("SeqDesignInferenceIncidKKClogit",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize a conditional logistic regression inference object for a completed KK design
		#' with a binary (incidence) response.
		#' @param seq_des_obj		A SeqDesign object (must be a KK design) whose entire n subjects
		#' 							are assigned and whose binary response y is recorded.
		#' @param num_cores			The number of CPU cores to use for randomization-based inference.
		#' 							Ignored for MLE-based inference. Default is 1.
		#' @param verbose			Whether to print progress messages. Default is \code{FALSE}.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "incidence")
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop("SeqDesignInferenceIncidKKClogit requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Computes the conditional log-odds ratio estimate for the treatment effect.
		#'
		#' @return The treatment effect estimate (log-odds ratio) from the conditional logistic regression.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignKK14$new(n = 20, response_type = "incidence")
		#' for (i in 1 : 20){
		#'   seq_des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
		#' }
		#' seq_des$add_all_subject_responses(rbinom(20, 1, 0.5))
		#'
		#' seq_des_inf = SeqDesignInferenceIncidKKClogit$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a (1 - alpha)-level Wald confidence interval for the treatment log-odds ratio
		#' using the asymptotic normality of the conditional MLE.
		#'
		#' @param alpha		The significance level. Default is 0.05, yielding a 95\% CI.
		#'
		#' @return A named numeric vector of length 2 with the CI bounds.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des_inf$compute_mle_confidence_interval()
		#' }
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the two-sided Wald p-value for the treatment log-odds ratio.
		#'
		#' @param delta		The null value of the treatment effect to test against. Default is 0.
		#'
		#' @return The two-sided p-value.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
		#' }
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				stop("TO-DO")
			}
		}
	),

	private = list(

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)){
				return(invisible(NULL))
			}
			mod = private$fit_clogit()
			private$cached_values$beta_hat_T   = as.numeric(mod$coefficients["w"])
			private$cached_values$s_beta_hat_T = as.numeric(mod$se[1])
			private$cached_values$is_z = TRUE
		},

		fit_clogit = function(){
			match_indic = private$match_indic
			if (is.null(match_indic)){
				match_indic = rep(0L, private$n)
			}
			match_indic[is.na(match_indic)] = 0L

			# Conditional logistic regression only applies to matched pairs
			i_matched = which(match_indic > 0)
			if (length(i_matched) == 0){
				stop("No matched pairs found; bclogit::clogit requires at least one matched pair.")
			}

			y_m      = private$y[i_matched]
			w_m      = private$w[i_matched]
			strata_m = match_indic[i_matched]
			X_m      = as.data.frame(private$X[i_matched, , drop = FALSE])

			# formula = NULL dispatches to clogit.default; treatment and strata are
			# passed as separate vectors per the bclogit API (not inside the formula)
			bclogit::clogit(
				formula        = NULL,
				y              = y_m,
				X              = X_m,
				treatment      = w_m,
				strata         = strata_m,
				treatment_name = "w"
			)
		}
	)
)
