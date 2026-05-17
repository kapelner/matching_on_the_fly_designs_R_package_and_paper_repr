#' Parametric-Bootstrap-Capable Likelihood Inference
#'
#' Intermediate abstract base for the subset of likelihood-backed inference
#' families that are plausible targets for parametric null-bootstrap
#' likelihood-ratio calibration.
#'
#' This class sits between \code{InferenceAsympLik} and
#' \code{InferenceAsympLikStdModCache} in the hierarchy.  Families with highly
#' bespoke partial-likelihood, quadrature, frailty, copula, or custom
#' combined-likelihood geometry remain direct children of
#' \code{InferenceAsympLik} and do not pass through here.
#'
#' @keywords internal
InferenceParamBootstrap = R6::R6Class("InferenceParamBootstrap",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = list(
		#' @description Bootstrap-calibrated likelihood-ratio two-sided p-value.
		#'
		#' Fits the null model at \code{delta}, simulates \code{B} datasets from
		#' that fitted null, refits unrestricted and null models on each, and
		#' returns the empirical tail probability of the observed LR statistic.
		#'
		#' @param delta        Null treatment effect. Default 0.
		#' @param B            Number of bootstrap replicates. Default 199.
		#' @param show_progress Logical; show a progress bar. Default \code{FALSE}.
		#' @return A scalar p-value, or \code{NA_real_} if the computation fails.
		compute_lik_ratio_bootstrap_two_sided_pval = function(delta = 0, B = 199, show_progress = FALSE){
			if (!isTRUE(private$supports_lik_ratio_param_bootstrap())){
				stop(
					class(self)[1], " does not support parametric-bootstrap LR calibration. ",
					"Override private$supports_lik_ratio_param_bootstrap() and simulate_under_lik_null().",
					call. = FALSE
				)
			}
			if (should_run_asserts()){
				assertNumeric(delta, len = 1)
				assertCount(B, positive = TRUE)
			}
			spec = private$get_likelihood_test_spec()
			if (is.null(spec)) return(NA_real_)

			eval_obs = private$get_memoized_likelihood_test_eval(
				delta         = delta,
				testing_type  = "lik_ratio",
				spec          = spec,
				include_full_negloglik = TRUE,
				include_null_negloglik = TRUE
			)
			if (isTRUE(eval_obs$invalid)) return(NA_real_)
			if (!is.finite(eval_obs$full_negloglik) || !is.finite(eval_obs$null_negloglik)) return(NA_real_)
			lr_obs = 2 * (eval_obs$null_negloglik - eval_obs$full_negloglik)
			if (!is.finite(lr_obs)) return(NA_real_)
			null_fit = eval_obs$null_fit

			actual_cores = private$effective_parallel_cores("param_bootstrap", self$num_cores)

			run_one_lr = function(b){
				boot_spec = tryCatch(
					private$simulate_under_lik_null(spec, delta, null_fit),
					error = function(e) NULL
				)
				if (is.null(boot_spec)) return(NA_real_)
				full_nll_boot = tryCatch(
					boot_spec$neg_loglik(boot_spec$full_fit),
					error = function(e) NA_real_
				)
				if (!is.finite(full_nll_boot)) return(NA_real_)
				null_fit_boot = tryCatch(
					boot_spec$fit_null(delta),
					error = function(e) NULL
				)
				if (is.null(null_fit_boot)) return(NA_real_)
				null_nll_boot = tryCatch(
					boot_spec$neg_loglik(null_fit_boot),
					error = function(e) NA_real_
				)
				if (!is.finite(null_nll_boot)) return(NA_real_)
				lr_boot = 2 * (null_nll_boot - full_nll_boot)
				if (is.finite(lr_boot)) lr_boot else NA_real_
			}

			lr_boots = if (actual_cores <= 1L) {
				vapply(seq_len(B), run_one_lr, numeric(1))
			} else {
				chunk_n = max(1L, min(as.integer(actual_cores), as.integer(B)))
				chunk_id = ceiling(seq_len(B) / ceiling(B / chunk_n))
				chunks = split(seq_len(B), chunk_id)
				as.numeric(unlist(private$par_lapply(
					chunks,
					function(idxs) vapply(idxs, run_one_lr, numeric(1)),
					n_cores = actual_cores,
					budget = 1L,
					show_progress = show_progress
				), use.names = FALSE))
			}

			finite_lr = lr_boots[is.finite(lr_boots)]
			n_finite = length(finite_lr)
			if (n_finite == 0L) return(NA_real_)
			n_exceed = sum(finite_lr >= lr_obs)
			(1 + n_exceed) / (1 + n_finite)
		},

		#' @description Bootstrap-calibrated likelihood-ratio confidence interval.
		#'
		#' Inverts \code{compute_lik_ratio_bootstrap_two_sided_pval} via a
		#' bracket-and-bisect search seeded with the Wald interval.  Each p-value
		#' evaluation costs \code{B} bootstrap refits, so this method is
		#' substantially more expensive than the p-value alone.
		#'
		#' @param alpha        Significance level. Default 0.05.
		#' @param B            Bootstrap replicates per p-value evaluation. Default 199.
		#' @param show_progress Logical; show a progress bar. Default \code{FALSE}.
		#' @return Named two-element numeric vector with the confidence-interval bounds.
		compute_lik_ratio_bootstrap_confidence_interval = function(alpha = 0.05, B = 199, show_progress = FALSE){
			if (!isTRUE(private$supports_lik_ratio_param_bootstrap())){
				stop(
					class(self)[1], " does not support parametric-bootstrap LR calibration. ",
					"Override private$supports_lik_ratio_param_bootstrap() and simulate_under_lik_null().",
					call. = FALSE
				)
			}
			if (should_run_asserts()){
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				assertCount(B, positive = TRUE)
			}
			est = self$compute_estimate()
			if (!is.finite(est)) return(c(NA_real_, NA_real_))

			pval_fn = function(d) self$compute_lik_ratio_bootstrap_two_sided_pval(d, B = B, show_progress = show_progress)

			se = tryCatch(private$get_standard_error(), error = function(e) NA_real_)
			step = if (is.finite(se) && se > 0) se else max(abs(est), 1)
			step = max(step, 1e-4)
			wald_ci = tryCatch(
				private$compute_wald_confidence_interval_impl(alpha),
				error = function(e) c(NA_real_, NA_real_)
			)
			lower_seed = if (length(wald_ci) >= 1L && is.finite(wald_ci[[1L]])) wald_ci[[1L]] else NA_real_
			upper_seed = if (length(wald_ci) >= 2L && is.finite(wald_ci[[2L]])) wald_ci[[2L]] else NA_real_

			ci_vals = pval_invert_ci_cpp(
				pval_fn    = pval_fn,
				est        = est,
				alpha      = alpha,
				step       = step,
				lower_seed = lower_seed,
				upper_seed = upper_seed
			)
			ci = c(ci_vals[1L], ci_vals[2L])
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		}
	),
	private = list(
		supports_lik_ratio_param_bootstrap = function() FALSE,
		#' Simulate a bootstrap dataset under the fitted null likelihood and return
		#' a minimal spec list for refitting.
		#'
		#' Must be overridden by families that set supports_lik_ratio_param_bootstrap()
		#' to TRUE.  The returned list must contain at least:
		#'   - full_fit:  unrestricted fit on the simulated data
		#'   - fit_null:  function(delta, start) returning a constrained fit
		#'   - neg_loglik: function(fit) returning the neg-log-likelihood
		simulate_under_lik_null = function(spec, delta, null_fit){
			NULL
		}
	)
)
