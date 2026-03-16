#' Univariate Conditional Proportional-Odds Combined Inference for KK Designs
#'
#' @description
#' Fits a conditional proportional-odds model for ordinal responses under a KK
#' matching-on-the-fly design, including both matched pairs and reservoir subjects.
#' Each matched pair forms a stratum across thresholds. Each reservoir subject
#' is treated as its own stratum. This is analogous to a stratified proportional odds
#' model where the strata are the pairs and the individual reservoir subjects.
#'
#' @details
#' This estimator uses the conditional likelihood approach by expanding the ordinal
#' response into binary threshold indicators. Matched pairs contribute to the
#' conditional likelihood. Note that singleton strata (reservoir subjects) do not
#' contribute to the conditional logit likelihood directly unless combined with other
#' information or if the model is viewed as a stratified model. In this implementation,
#' we follow the standard conditional logit approach where only strata with variation
#' in both treatment and response contribute.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignKK14$
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
#' infer <- SeqDesignInferenceOrdinalUnivKKCondPropOddsCombinedRegr$
#'   new(seq_des, verbose =
#' FALSE)
#' infer
#'
SeqDesignInferenceOrdinalUnivKKCondPropOddsCombinedRegr = R6::R6Class("SeqDesignInferenceOrdinalUnivKKCondPropOddsCombinedRegr",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize a univariate conditional proportional-odds combined inference object.
		#' @param	seq_des_obj		A SeqDesign object (must be a KK design).
		#' @param	num_cores			Number of CPU cores.
		#' @param	verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "ordinal")
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design.")
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
		#' @param alpha Significance level.
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the MLE-based p-value.
		#' @param delta Null value.
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
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
				stop("Conditional proportional-odds combined KK estimator: could not compute a finite standard error.")
			}
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			match_indic = private$match_indic
			if (is.null(match_indic)) match_indic = rep(0L, private$n)
			match_indic[is.na(match_indic)] = 0L

			# Define strata
			# Matched pairs get their match_indic as stratum ID
			# Reservoir subjects (match_indic == 0) each get a unique ID
			strata_ids = match_indic
			reservoir_idx = which(strata_ids == 0L)
			if (length(reservoir_idx) > 0L){
				strata_ids[reservoir_idx] = max(strata_ids) + seq_along(reservoir_idx)
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
			n = private$n
			num_strata = max(strata_ids)

			y_stack = integer(n * n_thresholds)
			w_stack = integer(n * n_thresholds)
			strata_stack = integer(n * n_thresholds)

			for (k in seq_len(n_thresholds)){
				idx = ((k - 1L) * n + 1L):(k * n)
				y_stack[idx] = as.integer(y_ord > k)
				w_stack[idx] = private$w
				# Each threshold gets its own set of strata to maintain independence across thresholds 
				# if we were doing simple pooling, but for conditional logistic, we treat 
				# (pair, threshold) as the stratum.
				strata_stack[idx] = strata_ids + (k - 1L) * num_strata
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
