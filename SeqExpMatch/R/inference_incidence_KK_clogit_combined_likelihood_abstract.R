#' Abstract class for Conditional Logistic Combined-Likelihood Compound Inference
#'
#' Fits a single joint likelihood over all KK design data for incidence responses.
#' The matched-pair component uses the conditional logistic likelihood on discordant
#' pairs (signed within-pair differences); the reservoir component uses the marginal
#' logistic likelihood. Both share the treatment effect beta_T and covariate slopes
#' beta_xs.
#'
#' Combined design matrix (column layout: beta_0, beta_T, beta_xs):
#' Discordant pair rows: (0, t_diff_k, X_diff_k)
#' Reservoir rows:       (1, w_i, X_i)
#' The zero in column 1 for pair rows encodes that beta_0 is absent from the
#' conditional likelihood (eliminated by conditioning on the pair sum).
#' Fitting a single logistic regression on the stacked dataset maximises the
#' combined log-likelihood L = L_cond(pairs) + L_marg(reservoir).
#'
#' @keywords internal
InferenceAbstractKKClogitCombinedLikelihood = R6::R6Class("InferenceAbstractKKClogitCombinedLikelihood",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param num_cores			Number of CPU cores for parallel processing.
		#' @param verbose			Whether to print progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "incidence")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Returns the combined-likelihood estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		}
	),

	private = list(

		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		get_standard_error = function(){
			private$shared_combined_likelihood(estimate_only = FALSE)
			private$cached_values$s_beta_hat_T
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		# Fit the combined logistic likelihood over discordant matched-pair differences
		# and reservoir observations with SHARED covariate effects beta_xs.
		shared_combined_likelihood = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			p             = if (private$include_covariates()) ncol(private$X) else 0L
			has_reservoir = nRT > 0 && nRC > 0

			# ---- Build combined design matrix ------------------------------------
			X_comb   = NULL
			y_comb   = NULL
			j_beta_T = 2L

			if (m > 0){
				m_vec = private$m
				if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
				m_vec[is.na(m_vec)] = 0L
				i_matched = which(m_vec > 0)
				y_m      = private$y[i_matched]
				w_m      = private$w[i_matched]
				strata_m = m_vec[i_matched]
				X_mat    = if (p > 0L) as.matrix(private$X[i_matched, , drop = FALSE]) else matrix(nrow = length(y_m), ncol = 0L)

				if (has_reservoir){
					y_r    = KKstats$y_reservoir
					w_r    = KKstats$w_reservoir
					X_r    = if (p > 0L) as.matrix(KKstats$X_reservoir) else matrix(nrow = length(y_r), ncol = 0L)
					design = build_kk_combined_clogit_design_cpp(
						as.double(y_m), as.double(w_m), X_mat, as.integer(strata_m),
						as.double(y_r), as.double(w_r), X_r
					)
					X_comb   = design$X_comb
					y_comb   = design$y_comb
					j_beta_T = 2L
				} else {
					res = collect_discordant_pairs_cpp(
						as.double(y_m), as.double(w_m), X_mat, as.integer(strata_m)
					)
					if (res$nd > 0){
						X_comb   = if (p > 0L) cbind(res$t_diffs, res$X_diffs) else matrix(res$t_diffs, ncol = 1L)
						y_comb   = res$y_01
						j_beta_T = 1L
					}
				}
			} else if (has_reservoir){
				y_r    = KKstats$y_reservoir
				w_r    = KKstats$w_reservoir
				X_comb = if (p > 0L) cbind(1, w_r, as.matrix(KKstats$X_reservoir)) else cbind(1, w_r)
				y_comb = y_r
			}

			if (is.null(X_comb)){
				private$cached_values$beta_hat_T   = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			# ---- QR-reduce to full rank, preserving beta_T column ---------------
			qr_X = qr(X_comb)
			if (qr_X$rank < ncol(X_comb)){
				keep = qr_X$pivot[seq_len(qr_X$rank)]
				if (!(j_beta_T %in% keep)) keep[qr_X$rank] = j_beta_T
				keep     = sort(keep)
				X_comb   = X_comb[, keep, drop = FALSE]
				j_beta_T = which(keep == j_beta_T)
			}

			# ---- Fit combined logistic ------------------------------------------
			mod = tryCatch(
				if (estimate_only) {
					fast_logistic_regression(X_comb, y_comb)
				} else {
					fast_logistic_regression_with_var(X_comb, y_comb, j = j_beta_T)
				},
				error = function(e) NULL
			)
			if (is.null(mod) || !is.finite(mod$b[j_beta_T])){
				private$cached_values$beta_hat_T   = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T   = as.numeric(mod$b[j_beta_T])
			if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(mod$ssq_b_j)
			private$cached_values$is_z         = TRUE
			invisible(NULL)
		}
	)
)
