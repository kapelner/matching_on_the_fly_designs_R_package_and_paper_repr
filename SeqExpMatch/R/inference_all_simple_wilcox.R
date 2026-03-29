#' Wilcoxon Rank-Sum Inference for Any Sequential Design
#'
#' @description
#' Non-parametric inference based on the two-sample Wilcoxon rank-sum test
#' (Mann-Whitney U) for any sequential experimental design (Bernoulli, Efron, iBCRD,
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
InferenceAllSimpleWilcox = R6::R6Class("InferenceAllSimpleWilcox",
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj  A completed \code{DesignSeqOneByOne} object.
		#' @param num_cores    Number of CPU cores for parallel bootstrap/randomization. Default 1.
		#' @param verbose      Whether to print progress messages. Default \code{FALSE}.
		#' @examples
		#' set.seed(1)
		#' x_dat <- data.frame(
		#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
		#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
		#' )
		#' seq_des <- DesignSeqOneByOneBernoulli$new(n = nrow(x_dat), response_type = "continuous",
#'   verbose = FALSE)
		#' for (i in seq_len(nrow(x_dat))) {
		#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
		#' }
		#' seq_des$add_all_subject_responses(c(1.2, 0.9, 1.5, 1.8, 2.1, 1.7, 2.6, 2.2))
		#' infer <- InferenceAllSimpleWilcox$new(seq_des, verbose = FALSE)
		#' infer
		#'
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			res_type = des_obj$get_response_type()
			if (res_type == "incidence"){
				stop(
					"Wilcoxon rank-sum inference is not implemented for incidence (binary) ",
					"responses: the Hodges-Lehmann estimator is degenerate (almost always 0) ",
					"on 0/1 data. Use InferenceAllSimpleMeanDiff or a clogit estimator instead."
				)
			}
			assertResponseType(res_type, c("continuous", "count", "proportion", "survival", "ordinal"))
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			if (private$any_censoring){
				stop(
					"Wilcoxon rank-sum inference does not support censored survival data. ",
					"Use InferenceSurvivalGehanWilcox for censored survival outcomes."
				)
			}
		},

		#' @description
		#' Returns the Hodges-Lehmann pseudo-median of all pairwise treatment-minus-control
		#' differences.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a \eqn{1 - \alpha} asymptotic confidence interval based on the normal
		#' approximation for the HL estimator.
		#' @param alpha Significance level; default 0.05 gives a 95\% CI.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Returns a two-sided p-value for \eqn{H_0: \Delta = \delta} using the asymptotic
		#' normal approximation for the HL estimator.
		#' @param delta Null value; default 0.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		hl_point_estimate = function(y_vals, w_vals){
			wilcox_hl_point_estimate_cpp(as.numeric(y_vals), as.integer(w_vals))
		},

		compute_fast_bootstrap_distr = function(B, ...) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			# KK designs use design-aware resampling not available via these args; fall back to R loop.
			if (private$is_KK) return(NULL)

			# Simple (non-KK) bootstrap: args = (max_resample_attempts, n, y, dead, w)
			args = list(...)
			max_resample_attempts = args[[1]]
			n = args[[2]]
			y = args[[3]]
			dead = args[[4]]
			w = args[[5]]

			indices_mat = matrix(-1L, nrow = n, ncol = B)

			for (b in seq_len(B)) {
				attempt = 1L
				repeat {
					i_b = sample_int_replace_cpp(n, n)
					w_b = w[i_b]
					if (any(w_b == 1, na.rm = TRUE) && any(w_b == 0, na.rm = TRUE)) {
						indices_mat[, b] = i_b - 1L
						break
					}
					attempt = attempt + 1L
					if (attempt > max_resample_attempts) break
				}
			}

			compute_wilcox_hl_bootstrap_parallel_cpp(
				as.numeric(private$y),
				as.integer(private$w),
				indices_mat,
				private$num_cores
			)
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)

			# Optimization: w_mat is already pre-computed in generate_permutations
			w_mat = permutations$w_mat
			nsim = ncol(w_mat)

			y_sim = as.numeric(y)
			
			# Map transform_responses to transform_code
			t_code = 0L # none
			if (transform_responses == "log") {
				t_code = 1L
			} else if (transform_responses == "logit") {
				t_code = 2L
			} else if (transform_responses == "log1p") {
				t_code = 3L
			}
			
			res = compute_wilcox_hl_distr_parallel_cpp(y_sim, w_mat, as.numeric(delta), t_code, private$num_cores)
			return(res)
		},

		compute_treatment_estimate_during_randomization_inference = function(){
			if (is.null(private$custom_randomization_statistic_function)) {
				private$hl_point_estimate(private$y, private$w)
			} else {
				private$custom_randomization_statistic_function()
			}
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			yT = private$y[private$w == 1]
			yC = private$y[private$w == 0]

			if (length(yT) == 0L || length(yC) == 0L){
				private$cached_values$beta_hat_T   = NA_real_
			if (estimate_only) return(invisible(NULL))
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

			beta = private$hl_point_estimate(private$y, private$w)
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
