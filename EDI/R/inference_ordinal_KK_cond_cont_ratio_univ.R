#' Univariate Conditional Continuation-Ratio Inference for KK Designs
#'
#' Fits a conditional continuation-ratio logit model for ordinal responses under a KK
#' matching-on-the-fly design. Each matched pair is treated as a stratum.
#' Reservoir subjects each form their own unique stratum. The continuation-ratio
#' model logit(P(Y = k | Y >= k)) is fit by expanding the ordinal response into
#' binary continuation trials and applying conditional logistic regression.
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
#'   seq_des$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalUnivKKCondContRatioRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
InferenceOrdinalUnivKKCondContRatioRegr = R6::R6Class(
	"InferenceOrdinalUnivKKCondContRatioRegr",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize a univariate conditional continuation-ratio inference object.
		#' @param	des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param	num_cores			Number of CPU cores.
		#' @param	verbose			Whether to print progress messages.
		#' @param harden Whether to apply robustness measures.
		initialize = function(des_obj,  verbose = FALSE, harden = TRUE){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design.")
			}
			super$initialize(des_obj, verbose, harden)
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

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			# Define strata
			strata_ids = m_vec
			reservoir_idx = which(strata_ids == 0L)
			if (length(reservoir_idx) > 0L){
				strata_ids[reservoir_idx] = max(strata_ids) + seq_along(reservoir_idx)
			}

			y_ord = as.integer(factor(private$y, ordered = TRUE))
			K = max(y_ord)
			if (K < 2L){
				private$cache_nonestimable_estimate("ordinal_cond_cont_ratio_too_few_categories")
				return(invisible(NULL))
			}

			n_alpha = K - 1L
			n = private$n
			num_strata = max(strata_ids)

			# Data expansion for continuation ratio (Rcpp optimized)
			expanded = expand_continuation_ratio_data_cpp(as.integer(y_ord),
				as.integer(private$w), as.integer(strata_ids), as.integer(K))
			
			# Fit conditional logistic regression
			mod = clogit_helper(expanded$y, data.frame(), expanded$w, expanded$strata)
			
			if (is.null(mod) || !is.finite(mod$b[1]) || !is.finite(mod$ssq_b_j) || mod$ssq_b_j <= 0){
				if (private$harden && length(unique(expanded$strata)) > 0){
					private$cached_values$beta_hat_T   = 0
					private$cache_nonestimable_se("ordinal_cond_cont_ratio_standard_error_unavailable")
				} else {
					private$cache_nonestimable_estimate("ordinal_cond_cont_ratio_fit_unavailable")
				}
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T   = as.numeric(mod$b[1])
			private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
			private$cached_values$is_z         = TRUE
		}
	)
)
