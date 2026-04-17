#' Exact Binomial Incidence Inference for Matched-Pair Designs
#'
#' Performs exact matched-pair binomial inference for binary outcomes using only
#' discordant matched pairs. This class is available for
#' \code{FixedDesignBinaryMatch} and KK matching-on-the-fly designs. For KK
#' designs, only the matched-pair data are used and the reservoir is ignored.
#'
#' @export
InferenceIncidenceExactBinomial = R6::R6Class("InferenceIncidenceExactBinomial",
	lock_objects = FALSE,
	inherit = InferenceExact,
	public = list(
		#' @description
		#' Initialize exact matched-pair binomial inference for incidence outcomes.
		#' @param des_obj A completed design object.
		#' @param verbose Whether to print progress messages.
		#' @return A new \code{InferenceIncidenceExactBinomial} object.
		initialize = function(des_obj,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj, verbose)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (!private$design_supports_exact_binomial()) {
				stop("Exact binomial incidence inference requires FixedDesignBinaryMatch or KK matching designs.")
			}
			if (is(des_obj, "FixedDesignBinaryMatch")) {
				private$des_obj_priv_int$ensure_bms_computed()
			}
		},

		#' @description
		#' Compute the matched-pair treatment estimate on the log-odds scale.
		#' @param estimate_only Ignored for this estimator.
		#' @return The treatment estimate.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$get_exact_binomial_log_or_estimate()
		}
	),

	private = list(
		default_exact_type = "Binomial",

		resolve_exact_type = function(type){
			if (is.null(type)) type = private$default_exact_type
			if (should_run_asserts()) {
				assertChoice(type, c("Binomial"))
			}
			type
		},

		normalize_exact_inference_args = function(type, args_for_type = NULL){
			if (should_run_asserts()) {
				assertChoice(type, c("Binomial"))
				assertList(args_for_type, null.ok = TRUE)
			}
			utils::modifyList(setNames(list(list()), type), if (is.null(args_for_type)) list() else args_for_type)
		},

		assert_exact_inference_params = function(type, args_for_type){
			if (should_run_asserts()) {
				assertChoice(type, c("Binomial"))
				assertList(args_for_type)
				if (!(type %in% names(args_for_type))) stop("args_for_type must contain a list for ", type)
			}
			args = args_for_type[[type]]
			if (should_run_asserts()) {
				assertList(args)
				assertResponseType(private$des_obj$get_response_type(), "incidence")
				assertNoCensoring(private$any_censoring)
			}
			if (!private$design_supports_exact_binomial()) {
				stop("Exact binomial incidence inference requires FixedDesignBinaryMatch or KK matching designs.")
			}
			stats = private$get_exact_binomial_stats()
			if (should_run_asserts()) {
				if (stats$m <= 0L) {
					stop("Exact binomial incidence inference requires at least one matched pair.")
				}
				if (stats$d_plus + stats$d_minus <= 0L) {
					stop("Exact binomial incidence inference requires at least one discordant matched pair.")
				}
			}
			invisible(args)
		},

		compute_exact_confidence_interval_by_type = function(type, alpha, args_for_type){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				private$assert_exact_inference_params(type, args_for_type)
			}
			switch(type,
				Binomial = private$ci_exact_binomial(alpha)
			)
		},

		compute_exact_two_sided_pval_for_treatment_effect_by_type = function(type, delta, args_for_type){
			if (should_run_asserts()) {
				assertNumeric(delta, len = 1)
				private$assert_exact_inference_params(type, args_for_type)
			}
			switch(type,
				Binomial = private$pval_exact_binomial(delta)
			)
		},

		design_supports_exact_binomial = function(){
			is(private$des_obj, "FixedDesignBinaryMatch") || is(private$des_obj, "DesignSeqOneByOneKK14")
		},

		pval_exact_binomial = function(delta_0){
			stats = private$get_exact_binomial_stats()
			zhang_exact_binom_pval_cpp(stats$d_plus, stats$d_minus, delta_0)
		},

		ci_exact_binomial = function(alpha){
			stats = private$get_exact_binomial_stats()
			d_total = stats$d_plus + stats$d_minus
			ci_prob = stats::binom.test(stats$d_plus, d_total, conf.level = 1 - alpha)$conf.int
			ci = stats::qlogis(ci_prob)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},

		get_exact_binomial_log_or_estimate = function(){
			stats = private$get_exact_binomial_stats()
			if (stats$m <= 0L) return(NA_real_)
			log((stats$d_plus + 0.5) / (stats$d_minus + 0.5))
		},

		get_exact_binomial_stats = function(){
			if (!is.null(private$cached_values$incidence_exact_binomial_stats)) {
				return(private$cached_values$incidence_exact_binomial_stats)
			}

			if (is(private$des_obj, "FixedDesignBinaryMatch")) {
				private$des_obj_priv_int$ensure_bms_computed()
			}

			m_vec = private$des_obj_priv_int$m
			if (should_run_asserts()) {
				if (is.null(m_vec)) {
					stop("Matching structure is unavailable for exact binomial incidence inference.")
				}
			}
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L
			KKstats = compute_zhang_match_data_cpp(private$w, m_vec, private$y, private$get_X())
			stats = list(
				m = as.integer(KKstats$m),
				d_plus = as.integer(KKstats$d_plus),
				d_minus = as.integer(KKstats$d_minus)
			)
			private$cached_values$incidence_exact_binomial_stats = stats
			stats
		}
	)
)
