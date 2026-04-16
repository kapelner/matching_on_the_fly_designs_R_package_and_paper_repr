#' Cumulative Cauchit Inference for Ordinal Responses
#'
#' Cumulative Cauchit model inference for ordinal responses.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneBernoulli$new(n = nrow(x_dat), response_type = "ordinal",
#'   verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalUniCauchitRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
InferenceOrdinalUniCauchitRegr = R6::R6Class("InferenceOrdinalUniCauchitRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(
		#' @description
		#' Initialize a sequential experimental design estimation and test object.
		#' @param des_obj A DesignSeqOneByOne object whose entire n subjects are assigned and
		#'   response y is recorded within.
		#' @param verbose A flag indicating whether messages should be displayed.
		#' @param harden Whether to apply robustness measures.
		initialize = function(des_obj,  verbose = FALSE, harden = TRUE){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			super$initialize(des_obj, verbose, harden)
			assertNoCensoring(private$any_censoring)
		}
	),

	private = list(
		cauchit_polr_fallback = function(){
			if (!check_package_installed("MASS")) return(NULL)
			y_fac = factor(private$y, levels = sort(unique(private$y)))
			if (length(levels(y_fac)) < 2) return(NULL)
			dat = data.frame(y = y_fac, w = private$w)
			mod = tryCatch(
				MASS::polr(y ~ w, data = dat, method = "cauchit", Hess = TRUE),
				error = function(e) NULL
			)
			if (is.null(mod) || !"w" %in% names(stats::coef(mod))) return(NULL)
			coef_w = as.numeric(stats::coef(mod)["w"])
			var_w = tryCatch(vcov(mod)["w", "w"], error = function(e) NA_real_)
			ssq = if (is.finite(var_w) && var_w > 0) var_w else NA_real_
			list(b = c(NA, coef_w), ssq_b_2 = ssq)
		},

		generate_mod = function(estimate_only = FALSE){
			Xmm = matrix(private$w, ncol = 1)
			full_names = c("treatment")
			colnames(Xmm) = full_names[seq_len(ncol(Xmm))]

			res = fast_ordinal_cauchit_regression_with_var_cpp(X = Xmm, y = as.numeric(private$y))
			b1 = tryCatch(res$b[1], error = function(e) NA_real_)
			if (is.finite(b1)){
				if (estimate_only) return(list(b = c(NA, b1), ssq_b_2 = NA_real_))
				if (is.finite(res$ssq_b_2) && res$ssq_b_2 > 0) return(list(b = c(NA, b1), ssq_b_2 = res$ssq_b_2))
			}
			fallback = private$cauchit_polr_fallback()
			if (!is.null(fallback)){
				if (estimate_only) return(list(b = fallback$b, ssq_b_2 = NA_real_))
				return(fallback)
			}
			list(b = c(NA, b1), ssq_b_2 = if (estimate_only) NA_real_ else res$ssq_b_2)
		}
	)
)
#' Multivariate Cumulative Cauchit Inference for Ordinal Responses
#'
#' Multivariate cumulative cauchit model inference for ordinal responses.
#'
#' @export
InferenceOrdinalMultiCauchitRegr = R6::R6Class("InferenceOrdinalMultiCauchitRegr",
	lock_objects = FALSE,
	inherit = InferenceOrdinalUniCauchitRegr,
	public = list(
	),

	private = list(
		generate_mod = function(estimate_only = FALSE){
			X_full = cbind(treatment = private$w, private$get_X())
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					res = fast_ordinal_cauchit_regression_with_var_cpp(X = X_fit, y = as.numeric(private$y))
					b1 = tryCatch(res$b[1], error = function(e) NA_real_)
					if (is.finite(b1)){
						if (estimate_only) return(list(b = c(NA, res$b), ssq_b_2 = NA_real_, converged = res$converged))
						if (is.finite(res$ssq_b_2) && res$ssq_b_2 > 0) return(list(b = c(NA, res$b), ssq_b_2 = res$ssq_b_2, converged = res$converged))
					}
					fallback = private$cauchit_polr_fallback()
					if (!is.null(fallback)){
						if (estimate_only) return(list(b = fallback$b, ssq_b_2 = NA_real_, converged = TRUE))
						return(c(fallback, list(converged = TRUE)))
					}
					list(b = c(NA, res$b), ssq_b_2 = if (estimate_only) NA_real_ else res$ssq_b_2, converged = res$converged)
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
					if (!is.null(mod$converged) && !mod$converged) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0
				}
			)
			attempt$fit
		}
	)
)
