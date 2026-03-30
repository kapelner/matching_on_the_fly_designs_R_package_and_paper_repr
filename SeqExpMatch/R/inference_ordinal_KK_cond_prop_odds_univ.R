#' Univariate Conditional Proportional-Odds Inference for KK Designs
#'
#' Fits a conditional proportional-odds model for ordinal responses under a KK
#' matching-on-the-fly design. Each matched pair is treated as its own stratum,
#' and the cumulative-logit model is fit by expanding the ordinal response into
#' binary threshold indicators and then applying conditional logistic regression
#' to each pair-threshold stratum. Reservoir subjects are conceptually singleton
#' strata and therefore contribute no conditional likelihood information; they are
#' omitted from the expanded fit for numerical stability.
#'
#' @details
#' This estimator targets a common treatment log-odds shift across the cumulative
#' logits, analogous to a matched-pairs conditional logistic regression when the
#' ordinal response has only two levels.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneKK14$new(n = nrow(x_dat), response_type = "ordinal",
#' verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalUnivKKCondPropOddsRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
InferenceOrdinalUnivKKCondPropOddsRegr = R6::R6Class(
	"InferenceOrdinalUnivKKCondPropOddsRegr",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize a univariate conditional proportional-odds inference object for
		#' a completed KK design with an ordinal response.
		#' @param	des_obj		A DesignSeqOneByOne object (must be a KK design) whose entire n subjects
		#' are assigned and whose ordinal response y is recorded.
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
		#' @param	verbose			Whether to print progress messages. Default is \code{FALSE}.
		#' @param make_fork_cluster Whether to use a fork cluster for parallelization.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design.")
			}
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Returns the estimated treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the asymptotic confidence interval.
		#' @param alpha Significance level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the asymptotic p-value.
		#' @param delta Null value.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				stop("Testing non-zero delta is not yet implemented.")
			}
		}
	),

	private = list(
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL)) 
			}
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			# Standard conditional logit on pairs only
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			m_vec = private$m
			m = max(m_vec, na.rm = TRUE)
			if (m == 0){
				private$cached_values$beta_hat_T   = NA_real_
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			y_ord = as.integer(factor(private$y, ordered = TRUE))
			K = max(y_ord)
			if (K < 2L){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			n_thresholds = K - 1L
			
			# Filter to matched only
			matched_idx = which(m_vec > 0)
			y_m = y_ord[matched_idx]
			w_m = private$w[matched_idx]
			strata_m = m_vec[matched_idx]
			n_m = length(matched_idx)

			y_stack = integer(n_m * n_thresholds)
			w_stack = integer(n_m * n_thresholds)
			strata_stack = integer(n_m * n_thresholds)

			for (k in seq_len(n_thresholds)){
				idx = ((k - 1L) * n_m + 1L):(k * n_m)
				y_stack[idx] = as.integer(y_m > k)
				w_stack[idx] = w_m
				strata_stack[idx] = strata_m + (k - 1L) * m
			}

			# Fit conditional logistic regression on the stacked data
			mod = clogit_helper(y_stack, data.frame(), w_stack, strata_stack)
			
			if (is.null(mod) || !is.finite(mod$b[1]) || !is.finite(mod$ssq_b_j) || mod$ssq_b_j <= 0){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T   = as.numeric(mod$b[1])
			private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
			private$cached_values$is_z         = TRUE
		}
	)
)
