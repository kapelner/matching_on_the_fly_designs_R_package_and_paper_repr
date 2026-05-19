#' Jonckheere-Terpstra (JT) Test for Ordinal Responses
#'
#' Exact Jonckheere-Terpstra (JT) rank test for a two-arm ordered alternative with an
#' ordinal response. For treatment versus control, the test statistic is the
#' sum of Mann-Whitney U counts across groups. This class provides the exact
#' distribution-based p-value.
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
#' infer <- InferenceOrdinalJonckheereTerpstraTest$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
InferenceOrdinalJonckheereTerpstraTest = R6::R6Class(
	"InferenceOrdinalJonckheereTerpstraTest",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(
		#' @description Initialize the JT test object.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Returns the estimated treatment effect (JT superiority measure).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Returns the weighted JT superiority estimate for Bayesian-bootstrap re-estimation.
		#' @param subject_or_block_weights Bootstrap weights at the subject/block level.
		#' @param estimate_only If TRUE, skip exact p-value calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			private$cached_values$beta_hat_T = private$weighted_superiority(private$y, private$w, row_weights) - 0.5
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Returns the exact two-sided p-value.
		compute_exact_two_sided_pval_for_treatment_effect = function(){
			private$shared()
			private$cached_values$p_exact
		},
		#' @description Not applicable: JT test is exact, not asymptotic. Returns NA.
		#' @param alpha The significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			c(NA_real_, NA_real_)
		},
		#' @description Not applicable: JT test is exact. Use compute_exact_two_sided_pval_for_treatment_effect().
		#' @param delta The null treatment effect (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			NA_real_
		}
	),
	private = list(
		weighted_superiority = function(y_vals, w_vals, row_weights){
			y_vals = as.numeric(y_vals)
			w_vals = as.integer(w_vals)
			row_weights = as.numeric(row_weights)
			i_t = which(w_vals == 1L & is.finite(y_vals) & is.finite(row_weights) & row_weights > 0)
			i_c = which(w_vals == 0L & is.finite(y_vals) & is.finite(row_weights) & row_weights > 0)
			if (length(i_t) == 0L || length(i_c) == 0L) return(NA_real_)
			y_t = y_vals[i_t]
			y_c = y_vals[i_c]
			w_t = row_weights[i_t]
			w_c = row_weights[i_c]
			comp = outer(y_t, y_c, function(a, b) ifelse(a > b, 1, ifelse(a < b, 0, 0.5)))
			w_pair = outer(w_t, w_c, "*")
			den = sum(w_pair)
			if (!is.finite(den) || den <= 0) return(NA_real_)
			sum(comp * w_pair) / den
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			
			res = exact_jonckheere_terpstra_pval_cpp(as.integer(private$y), as.integer(private$w))
			
			private$cached_values$superiority = res$superiority
			private$cached_values$beta_hat_T = res$superiority - 0.5
			if (estimate_only) return(invisible(NULL))
			private$cached_values$p_exact = res$p_exact
			private$cached_values$p_lower = res$p_lower
			private$cached_values$p_upper = res$p_upper
		}
	)
)
