#' Stereotype Logistic Inference for Ordinal Responses
#'
#' @description
#' Stereotype logistic inference for ordinal responses. The model allows
#' category-specific treatment effects through latent score weights while keeping
#' a common treatment coefficient.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignBernoulli$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "ordinal",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- DesignInferenceOrdinalUniStereotypeLogitRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
DesignInferenceOrdinalUniStereotypeLogitRegr = R6::R6Class("DesignInferenceOrdinalUniStereotypeLogitRegr",
	inherit = DesignInferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj A SeqDesign object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param num_cores The number of CPU cores to use to parallelize
		#'   the sampling during randomization-based inference and
		#'   bootstrap resampling.
		#'   The default is 1 for serial computation. For simple
		#'   estimators (e.g. mean difference and KK compound),
		#'   parallelization is achieved with zero-overhead C++ OpenMP.
		#'   For complex models (e.g. GLMs),
		#'   parallelization falls back to R's
		#'   \code{parallel::mclapply}, which incurs
		#'   session-forking overhead.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{TRUE}.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			super$initialize(des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Computes the stereotype logistic treatment effect estimate.
		#'
		#' @return	The estimated treatment coefficient.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a profile-likelihood confidence interval for the stereotype
		#' treatment effect.
		#'
		#' @param alpha	The confidence level in the computed confidence interval is
		#' 				\eqn{1 - \alpha}. The default is \code{0.05}.
		#'
		#' @return	A profile-likelihood confidence interval for the treatment effect.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			beta_hat = self$compute_treatment_estimate()
			ll_hat = private$profile_loglik(beta_hat)
			cutoff = stats::qchisq(1 - alpha, df = 1)
			target = ll_hat - cutoff / 2

			find_bound = function(direction){
				step = 0.25
				edge = beta_hat + direction * step
				for (iter in 1:40){
					ll_edge = private$profile_loglik(edge)
					if (is.finite(ll_edge) && ll_edge <= target){
						f = function(b) private$profile_loglik(b) - target
						return(stats::uniroot(f, lower = min(beta_hat, edge), upper = max(beta_hat, edge))$root)
					}
					step = step * 2
					edge = beta_hat + direction * step
				}
				NA_real_
			}

			stats::setNames(c(find_bound(-1), find_bound(1)), c("2.5%", "97.5%"))
		},

		#' @description
		#' Computes a profile-likelihood ratio p-value for the stereotype treatment
		#' effect.
		#'
		#' @param delta	The null treatment effect. Only \code{0} is supported.
		#'
		#' @return	The approximate frequentist p-value.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			if (!identical(delta, 0)) stop("Only delta = 0 is currently supported for stereotype logistic inference.")
			beta_hat = self$compute_treatment_estimate()
			ll_hat = private$profile_loglik(beta_hat)
			ll_null = private$profile_loglik(0)
			lr_stat = max(0, 2 * (ll_hat - ll_null))
			stats::pchisq(lr_stat, df = 1, lower.tail = FALSE)
		}
	),

	private = list(
		stereotype_design_matrix = function(){
			Xmm = matrix(private$w, ncol = 1)
			colnames(Xmm) = "treatment"
			Xmm
		},

		generate_mod = function(){
			res = fast_stereotype_logit_with_var_cpp(
				X = private$stereotype_design_matrix(),
				y = as.numeric(private$y)
			)
			if (isTRUE(res$converged) && is.finite(res$ssq_b_1)){
				return(list(
					b = c(NA, res$b[1]),
					ssq_b_2 = res$ssq_b_1
				))
			}
			fallback = private$stereotype_polr_fallback()
			if (!is.null(fallback)){
				return(fallback)
			}
			list(
				b = c(NA, res$b[1]),
				ssq_b_2 = res$ssq_b_1
			)
		},

		profile_loglik = function(beta){
			fast_stereotype_profile_loglik_cpp(
				X = private$stereotype_design_matrix(),
				y = as.numeric(private$y),
				beta_fixed = beta
			)
		},

		stereotype_polr_fallback = function(link = "logistic"){
			if (!requireNamespace("MASS", quietly = TRUE)) return(NULL)
			y_fac = factor(private$y, levels = sort(unique(private$y)))
			if (length(levels(y_fac)) < 2) return(NULL)
			dat = data.frame(y = y_fac, w = private$w)
			mod = tryCatch(
				MASS::polr(y ~ w, data = dat, method = link, Hess = TRUE),
				error = function(e) NULL
			)
			if (is.null(mod) || !"w" %in% names(stats::coef(mod))) return(NULL)
			coef_w = as.numeric(stats::coef(mod)["w"])
			var_w = tryCatch(vcov(mod)["w", "w"], error = function(e) NA_real_)
			ssq = if (is.finite(var_w) && var_w > 0) var_w else NA_real_
			list(b = c(NA, coef_w), ssq_b_2 = ssq)
		}
	)
)
