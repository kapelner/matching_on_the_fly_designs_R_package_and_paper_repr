#' Paired Sign Test Inference for KK Designs with Ordinal Response
#'
#' Fits a paired sign test for ordinal responses under a KK matching-on-the-fly
#' design that stores ordinary matched pairs. For matched pairs, it considers the
#' sign of the within-pair differences. Reservoir subjects are not included in this
#' simple paired test.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneKK21$new(n = nrow(x_dat), response_type = "ordinal",
#' verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalPairedSignTest$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
InferenceOrdinalPairedSignTest = R6::R6Class("InferenceOrdinalPairedSignTest",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = utils::modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize the inference object.
		#' @param  des_obj  	A completed KK matching-on-the-fly design object.
		#' @param  verbose  		Whether to print progress messages.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_passthrough(des_obj)
		},
		#' @description Returns the estimated treatment effect (proportion of pairs where T > C).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Recomputes the paired-sign treatment estimate under
		#'   Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Subject-, block-, cluster-, or matched-set
		#'   bootstrap weights.
		#' @param estimate_only If \code{TRUE}, compute only the weighted point
		#'   estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = as.numeric(private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights))
			m_vec = private$m
			if (is.null(m_vec)) {
				private$compute_basic_match_data()
				m_vec = private$m
			}
			if (is.null(m_vec)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			pair_ids = sort(unique(m_vec[m_vec > 0L]))
			if (!length(pair_ids)) {
				private$cache_nonestimable_estimate("ordinal_paired_sign_test_no_discordant_pairs")
				return(NA_real_)
			}
			pos_w = 0
			neg_w = 0
			for (pid in pair_ids) {
				idx = which(m_vec == pid)
				if (length(idx) < 2L) next
				i_t = idx[private$w[idx] == 1L][1L]
				i_c = idx[private$w[idx] == 0L][1L]
				if (!is.finite(i_t) || !is.finite(i_c)) next
				diff = as.numeric(private$y[i_t]) - as.numeric(private$y[i_c])
				pair_w = mean(row_weights[c(i_t, i_c)])
				if (!is.finite(pair_w) || pair_w <= 0) next
				if (diff > 0) pos_w = pos_w + pair_w
				if (diff < 0) neg_w = neg_w + pair_w
			}
			n_eff = pos_w + neg_w
			if (!is.finite(n_eff) || n_eff <= 0) {
				private$cached_values$beta_hat_T = 0
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(private$cached_values$beta_hat_T)
			}
			p_hat = pos_w / n_eff
			private$cached_values$beta_hat_T = as.numeric(p_hat - 0.5)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Computes the confidence interval for the probability P(T > C).
		#' @param  alpha  				The significance level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes the p-value for the sign test.
		#' @param  delta  				The null difference (must be 0 for sign test).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
				if (delta != 0) stop("Sign test only supports testing against delta = 0.")
			}
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
		}
	)),
	private = utils::modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		supports_likelihood_tests = function() FALSE,
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			
			diffs = private$cached_values$KKstats$y_matched_diffs
			# Sign test on matched pairs: ignore ties (diff == 0)
			pos = sum(diffs > 0)
			neg = sum(diffs < 0)
			n_eff = pos + neg
			
			if (n_eff == 0){
				if (private$harden && length(diffs) > 0){
					# If all pairs are tied, the most natural estimate is 0 (p_hat = 0.5)
					# but we have no variance information.
					private$cached_values$beta_hat_T = 0
					private$cache_nonestimable_se("ordinal_paired_sign_test_no_discordant_pairs")
				} else {
					private$cache_nonestimable_estimate("ordinal_paired_sign_test_no_discordant_pairs")
				}
				return(invisible(NULL))
			}
			
			# Estimate: proportion of non-tied pairs favoring treatment
			p_hat = pos / n_eff
			# Beta is usually defined as p_hat - 0.5 for centered tests
			private$cached_values$beta_hat_T = p_hat - 0.5
			# Standard error for proportion
			se = sqrt(p_hat * (1 - p_hat) / n_eff)
			
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
		}
	))
)
