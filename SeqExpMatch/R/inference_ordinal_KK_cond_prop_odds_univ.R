#' Univariate Conditional Proportional-Odds Inference for KK Designs
#'
#' @description
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
SeqDesignInferenceOrdinalUnivKKCondPropOddsRegr = R6::R6Class("SeqDesignInferenceOrdinalUnivKKCondPropOddsRegr",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize a univariate conditional proportional-odds inference object for
		#' a completed KK design with an ordinal response.
		#' @param	seq_des_obj		A SeqDesign object (must be a KK design) whose entire n subjects
		#' 							are assigned and whose ordinal response y is recorded.
		#' @param	num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference
		#' 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs),
		#' 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
		#' @param	verbose			Whether to print progress messages. Default is \code{FALSE}.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "ordinal")
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Returns the estimated treatment effect.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the MLE-based confidence interval.
		#' @param alpha The confidence level in the computed confidence interval is
		#'   \code{1 - alpha}. The default is 0.05.
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the MLE-based p-value.
		#' @param delta The null difference to test against. For any treatment effect
		#'   at all this is set to zero (the default).
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				stop("TO-DO")
			}
		}
	),

	private = list(
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop("Conditional proportional-odds KK estimator: could not compute a finite standard error.")
			}
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			match_indic = private$match_indic
			if (is.null(match_indic)) match_indic = rep(0L, private$n)
			match_indic[is.na(match_indic)] = 0L

			i_matched = which(match_indic > 0L)
			if (length(i_matched) == 0L){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			pair_id = match_indic[i_matched]
			pair_sizes = table(pair_id)
			complete_pair_ids = as.integer(names(pair_sizes)[pair_sizes == 2L])
			if (length(complete_pair_ids) == 0L){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			keep = pair_id %in% complete_pair_ids
			y_m = as.integer(factor(private$y[i_matched][keep], ordered = TRUE))
			w_m = as.integer(private$w[i_matched][keep])
			pair_id = pair_id[keep]

			order_idx = order(pair_id)
			y_m = y_m[order_idx]
			w_m = w_m[order_idx]
			pair_id = pair_id[order_idx]

			m = length(unique(pair_id))
			K = length(unique(y_m))
			if (K < 2L){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			pair_index = match(pair_id, unique(pair_id))
			n_matched = length(y_m)
			n_thresholds = K - 1L

			y_stack = integer(n_matched * n_thresholds)
			w_stack = integer(n_matched * n_thresholds)
			strata_stack = integer(n_matched * n_thresholds)

			for (k in seq_len(n_thresholds)){
				idx = ((k - 1L) * n_matched + 1L):(k * n_matched)
				y_stack[idx] = as.integer(y_m > k)
				w_stack[idx] = w_m
				strata_stack[idx] = pair_index + (k - 1L) * m
			}

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
