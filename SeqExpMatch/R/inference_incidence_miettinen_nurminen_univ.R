#' Univariate Miettinen-Nurminen Risk-Difference Inference for Binary Responses
#'
#' @description
#' Fits the classical Miettinen-Nurminen score method for the risk difference in a
#' two-arm binary trial. The point estimate is the observed risk difference, while
#' confidence intervals and p-values are obtained by inverting the constrained
#' score test under the null \eqn{p_T - p_C = \delta}.
#'
#' @details
#' This class is intentionally unadjusted. It operates on the \eqn{2 \times 2}
#' table induced by treatment assignment and incidence response, and is therefore
#' the natural classical binary-endpoint complement to the regression-based
#' incidence methods already in the package.
#'
#' @export
InferenceIncidUnivMiettinenNurminenRiskDiff = R6::R6Class("InferenceIncidUnivMiettinenNurminenRiskDiff",
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize a Miettinen-Nurminen risk-difference inference object for a
		#' completed design with an incidence response.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object with an incidence response.
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = DesignSeqOneByOneBernoulli$new(n = 20, response_type = "incidence")
		#' for (i in 1:20) {
		#' 	x_i = data.frame(x1 = rnorm(1), x2 = rnorm(1))
		#' 	w_i = seq_des$add_subject_to_experiment_and_assign(x_i)
		#' 	p_i = plogis(-0.8 + 0.5 * w_i)
		#' 	seq_des$add_subject_response(i, rbinom(1, 1, p_i))
		#' }
		#' seq_des_inf = InferenceIncidUnivMiettinenNurminenRiskDiff$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "incidence")
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Computes the observed risk-difference estimate.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a 1 - \code{alpha} Miettinen-Nurminen score confidence interval.
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		#' @param pval_epsilon Bisection tolerance for CI bounds.
		compute_asymp_confidence_interval = function(alpha = 0.05, pval_epsilon = 1e-7){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			counts = private$cached_values$mn_counts
			if (is.null(counts)) return(c(NA_real_, NA_real_))
			
			ci = mn_ci_cpp(counts$x_t, counts$n_t, counts$x_c, counts$n_c, counts$p_t, counts$p_c, alpha, pval_epsilon)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},

		#' @description
		#' Computes a two-sided Miettinen-Nurminen score p-value.
		#' @param delta The null treatment effect on the risk-difference scale.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta, len = 1)
			private$shared()
			counts = private$cached_values$mn_counts
			if (is.null(counts)) return(NA_real_)
			
			res = mn_pvalue_cpp(counts$x_t, counts$n_t, counts$x_c, counts$n_c, delta, counts$p_t, counts$p_c)
			if (!is.finite(res)) stop("Miettinen-Nurminen risk-difference: could not compute a finite score statistic.")
			res
		}
	),

	private = list(
		mn_eps = function(){
			sqrt(.Machine$double.eps)
		},

		get_counts = function(){
			i_t = private$w == 1
			i_c = private$w == 0
			n_t = sum(i_t)
			n_c = sum(i_c)
			x_t = sum(private$y[i_t])
			x_c = sum(private$y[i_c])

			list(
				n_t = n_t,
				n_c = n_c,
				x_t = x_t,
				x_c = x_c,
				p_t = if (n_t > 0) x_t / n_t else NA_real_,
				p_c = if (n_c > 0) x_c / n_c else NA_real_
			)
		},

		set_failed_fit_cache = function(){
			private$cached_values$mn_counts = NULL
			private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$is_z = TRUE
			private$cached_values$df = NA_real_
			private$cached_values$full_coefficients = NULL
			private$cached_values$full_vcov = NULL
			private$cached_values$summary_table = NULL
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			counts = private$get_counts()
			if (counts$n_t == 0 || counts$n_c == 0 || !is.finite(counts$p_t) || !is.finite(counts$p_c)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			est = counts$p_t - counts$p_c
			correction = (counts$n_t + counts$n_c) / max(counts$n_t + counts$n_c - 1, 1)
			se = sqrt(correction * (
				counts$p_t * (1 - counts$p_t) / counts$n_t +
				counts$p_c * (1 - counts$p_c) / counts$n_c
			))

			private$cached_values$mn_counts = counts
			private$cached_values$beta_hat_T = est
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z = TRUE
			private$cached_values$df = counts$n_t + counts$n_c - 2
			private$cached_values$full_coefficients = c("(Intercept)" = counts$p_c, "treatment" = est)
			private$cached_values$full_vcov = matrix(
				c(
					counts$p_c * (1 - counts$p_c) / counts$n_c,
					0,
					0,
					se^2
				),
				nrow = 2,
				dimnames = list(c("(Intercept)", "treatment"), c("(Intercept)", "treatment"))
			)
			z_val = if (is.finite(se) && se > 0) est / se else NA_real_
			private$cached_values$summary_table = rbind(
				treatment = c(
					Value = est,
					`Std. Error` = se,
					`z value` = z_val,
					`Pr(>|z|)` = if (is.finite(z_val)) 2 * stats::pnorm(-abs(z_val)) else NA_real_
				)
			)
		}
	)
)
