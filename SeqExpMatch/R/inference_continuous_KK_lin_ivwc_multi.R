#' Lin (2013) IVWC Inference for KK Designs with Continuous Responses
#'
#' @description
#' Fits a Lin (2013) style covariate-adjusted estimator for KK
#' matching-on-the-fly designs with continuous responses by combining matched-pair
#' and reservoir estimators via inverse-variance weighting.
#'
#' The reservoir component uses the usual Lin regression with centered covariates
#' and treatment-by-centered-covariate interactions. The matched-pair component
#' uses pair differences together with pair covariate differences and centered pair
#' means, which is the algebraically equivalent paired Lin regression.
#'
#' @export
InferenceContinMultiKKLinIVWC = R6::R6Class("InferenceContinMultiKKLinIVWC",
	inherit = InferenceKKPassThroughCompound,
	public = list(

		#' @description
		#' Initialize the inference object.
		#'
		#' @param des_obj A completed KK \code{DesignSeqOneByOne} object with a continuous
		#'   response.
		#' @param num_cores The number of CPU cores to use for bootstrap and
		#'   randomization inference.
		#' @param verbose Whether to print progress messages.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = DesignSeqOneByOneKK14$new(n = 20, response_type = "continuous")
		#' for (t in 1:20) {
		#' 	x_t = data.frame(x1 = rnorm(1), x2 = rnorm(1))
		#' 	w_t = seq_des$add_subject_to_experiment_and_assign(x_t)
		#' 	seq_des$add_subject_response(t, x_t$x1 + 0.5 * w_t + rnorm(1))
		#' }
		#' seq_des_inf = InferenceContinMultiKKLinIVWC$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "continuous")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Computes the IVWC Lin estimate of the treatment effect.
		#'
		#' @return The estimated treatment effect.
		compute_treatment_estimate = function(){
			if (is.null(private$cached_values$beta_T_reservoir) && is.null(private$cached_values$beta_T_matched)){
				private$compute_estimate_from_matched_and_reservoir(
					private$lin_for_matched_pairs,
					private$lin_for_reservoir
				)
			}
			if (is.null(private$cached_values$beta_hat_T)){
				ssq_m = private$cached_values$ssq_beta_T_matched
				ssq_r = private$cached_values$ssq_beta_T_reservoir
				private$cached_values$beta_hat_T =
					if (private$only_matches()){
						private$cached_values$beta_T_matched
					} else if (private$only_reservoir()){
						private$cached_values$beta_T_reservoir
					} else if (!is.finite(ssq_m) || ssq_m <= 0) {
						private$cached_values$beta_T_reservoir
					} else if (!is.finite(ssq_r) || ssq_r <= 0) {
						private$cached_values$beta_T_matched
					} else {
						w_star = ssq_r / (ssq_r + ssq_m)
						w_star * private$cached_values$beta_T_matched + (1 - w_star) * private$cached_values$beta_T_reservoir
					}
			}
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a 1 - \code{alpha} confidence interval.
		#'
		#' @param alpha The confidence level in the computed confidence interval is
		#'   1 - \code{alpha}. The default is 0.05.
		#'
		#' @return A confidence interval for the treatment effect.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared_for_inference()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes a two-sided p-value.
		#'
		#' @param delta The null treatment effect. Defaults to 0.
		#'
		#' @return The approximate p-value.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared_for_inference()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		compute_basic_match_data = function(){
			if (is.null(private$X)){
				private$X = private$get_X()
			}
			private$cached_values$KKstats = .compute_kk_lin_basic_match_data(
				X = private$X,
				n = private$n,
				y = private$y,
				w = private$w,
				m_vec = private$m
			)
		},

		shared_for_inference = function(){
			if (is.null(private$cached_values$beta_hat_T)){
				self$compute_treatment_estimate()
			}
			ssq_beta_hat_T =
				if (private$only_matches()){
					ssq_m = private$cached_values$ssq_beta_T_matched
					if (!is.finite(ssq_m) || ssq_m <= 0) {
						if (is.null(private$cached_values$ssq_beta_T_reservoir)) {
							private$lin_for_reservoir()
						}
						private$cached_values$ssq_beta_T_reservoir
					} else {
						ssq_m
					}
				} else if (private$only_reservoir()){
					ssq_r = private$cached_values$ssq_beta_T_reservoir
					if (!is.finite(ssq_r) || ssq_r <= 0) {
						if (is.null(private$cached_values$ssq_beta_T_matched) && private$cached_values$KKstats$m > 0) {
							private$lin_for_matched_pairs()
						}
						private$cached_values$ssq_beta_T_matched
					} else {
						ssq_r
					}
				} else {
					ssq_m = private$cached_values$ssq_beta_T_matched
					ssq_r = private$cached_values$ssq_beta_T_reservoir
					if (!is.finite(ssq_m) || ssq_m <= 0) {
						ssq_r
					} else if (!is.finite(ssq_r) || ssq_r <= 0) {
						ssq_m
					} else {
						ssq_m * ssq_r / (ssq_m + ssq_r)
					}
				}

			private$cached_values$s_beta_hat_T = sqrt(ssq_beta_hat_T)
			private$cached_values$is_z = TRUE
		},

		only_matches = function(){
			private$cached_values$KKstats$nRT <= 2 || private$cached_values$KKstats$nRC <= 2
		},

		only_reservoir = function(){
			private$cached_values$KKstats$m == 0
		},

		reduce_design_matrix_preserving_required = function(X_full, required){
			qr_X = qr(X_full)
			target_rank = qr_X$rank
			candidate_order = c(required, setdiff(qr_X$pivot, required))
			keep = integer(0)

			for (j in candidate_order){
				trial_keep = c(keep, j)
				trial_rank = qr(X_full[, trial_keep, drop = FALSE])$rank
				if (trial_rank > length(keep)){
					keep = trial_keep
				}
				if (length(keep) >= target_rank){
					break
				}
			}

			sort(unique(keep))
		},

		fit_hc2_ols = function(X_full, y, required, j_target){
			keep = private$reduce_design_matrix_preserving_required(X_full, required)
			if (!(j_target %in% keep)) return(NULL)

			X_fit = X_full[, keep, drop = FALSE]
			j_fit = match(j_target, keep)
			if (nrow(X_fit) <= ncol(X_fit)) return(NULL)

			mod = stats::lm.fit(X_fit, y)
			coef_hat = as.numeric(mod$coefficients)
			if (length(coef_hat) != ncol(X_fit) || any(!is.finite(coef_hat))) return(NULL)

			post_fit = tryCatch(
				ols_hc2_post_fit_cpp(X_fit, y, coef_hat, j_fit),
				error = function(e) NULL
			)
			if (is.null(post_fit)) return(NULL)

			beta = post_fit$beta_hat
			ssq = post_fit$ssq_hat
			if (!is.finite(beta) || !is.finite(ssq) || ssq <= 0) return(NULL)

			list(beta = beta, ssq = ssq)
		},

		lin_for_matched_pairs = function(){
			KKstats = private$cached_values$KKstats
			yd = KKstats$y_matched_diffs
			m = length(yd)
			if (m == 0){
				private$cached_values$beta_T_matched = NA_real_
				private$cached_values$ssq_beta_T_matched = NA_real_
				return(invisible(NULL))
			}

			xbar = colMeans(private$get_X())
			Xd = as.matrix(KKstats$X_matched_diffs_full)
			Xm_centered = sweep(as.matrix(KKstats$X_matched_means_full), 2, xbar, "-")
			X_full = cbind(1, Xd, Xm_centered)
			j_target = 1L

			fit = private$fit_hc2_ols(X_full, yd, required = 1L, j_target = j_target)
			if (is.null(fit)) {
				private$cached_values$beta_T_matched = if (m >= 1) mean(yd) else NA_real_
				private$cached_values$ssq_beta_T_matched = if (m >= 2) var(yd) / m else NA_real_
			} else {
				private$cached_values$beta_T_matched = fit$beta
				private$cached_values$ssq_beta_T_matched = fit$ssq
			}
		},

		lin_for_reservoir = function(){
			KKstats = private$cached_values$KKstats
			y_r = KKstats$y_reservoir
			w_r = KKstats$w_reservoir
			X_r = as.matrix(KKstats$X_reservoir)
			if (length(y_r) == 0 || !any(w_r == 1) || !any(w_r == 0)){
				private$cached_values$beta_T_reservoir = NA_real_
				private$cached_values$ssq_beta_T_reservoir = NA_real_
				return(invisible(NULL))
			}

			xbar = colMeans(private$get_X())
			Xc = sweep(X_r, 2, xbar, "-")
			X_int = Xc * w_r
			X_full = cbind(1, w_r, Xc, X_int)
			j_target = 2L

			fit = private$fit_hc2_ols(X_full, y_r, required = c(1L, 2L), j_target = j_target)
			if (is.null(fit)) {
				y_rT = y_r[w_r == 1]
				y_rC = y_r[w_r == 0]
				nRT = length(y_rT)
				nRC = length(y_rC)
				private$cached_values$beta_T_reservoir = mean(y_rT) - mean(y_rC)
				private$cached_values$ssq_beta_T_reservoir =
					if (nRT > 1 && nRC > 1) var(y_rT) / nRT + var(y_rC) / nRC else NA_real_
			} else {
				private$cached_values$beta_T_reservoir = fit$beta
				private$cached_values$ssq_beta_T_reservoir = fit$ssq
			}
		}
	)
)
