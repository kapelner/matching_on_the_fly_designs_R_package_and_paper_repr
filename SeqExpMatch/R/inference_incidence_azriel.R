#' Azriel Blocked Incidence Inference
#'
#' Unadjusted blocked-design incidence inference using the simple mean-difference
#' point estimate with a block-stratified standard error.
#'
#' @export
InferenceIncidAzriel = R6::R6Class("InferenceIncidAzriel",
	lock_objects = FALSE,
	inherit = InferenceAllSimpleMeanDiff,
	public = list(
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "incidence")
			if (!private$design_is_supported_blocking(des_obj)) {
				stop(class(self)[1], " requires a blocking design.")
			}
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
			if (is.null(private$get_block_ids())) {
				stop(class(self)[1], " requires block identifiers from the completed design.")
			}
		}
	),

	private = list(
		m = NULL,

		design_is_supported_blocking = function(des_obj){
			is(des_obj, "DesignSeqOneByOneKK14") ||
				is(des_obj, "FixedDesignBlocking") ||
				is(des_obj, "DesignSeqOneByOneSPBR") ||
				is(des_obj, "DesignSeqOneByOneRandomBlockSize")
		},

		get_block_ids = function(){			
			if (is.null(private$m)) {
				private$m = tryCatch(private$des_obj$get_m(), error = function(e) private$des_obj_priv_int$m)
				strata_cols = private$des_obj_priv_int$strata_cols
				Xraw = private$des_obj_priv_int$Xraw
				if (!is.null(strata_cols) && length(strata_cols) > 0L && nrow(Xraw) == length(private$y)) {
					strata_keys = vapply(seq_len(nrow(Xraw)), function(i) {
						vals = vapply(strata_cols, function(col) {
							val = Xraw[i, ][[col]]
							if (is.na(val)) "NA" else as.character(val)
						}, character(1))
						paste(vals, collapse = "|")
					}, character(1))
					private$m = match(strata_keys, unique(strata_keys))
				}
			}
			if (is.null(private$m)) private$m = NULL
			private$m = as.integer(private$m)
			if (length(private$m) != length(private$y)) private$m = NULL
			private$m
		},

		get_standard_error = function(){
			if (!is.null(private$cached_values$azriel_s_beta_hat_T)) {
				return(private$cached_values$azriel_s_beta_hat_T)
			}
			if (is.null(private$m)){
				return(NA_real_)
			}
			private$cached_values$azriel_s_beta_hat_T = compute_azriel_block_se_cpp(
				private$des_obj_priv_int$y,
				private$m,
				private$des_obj_priv_int$n
			)
			private$cached_values$azriel_s_beta_hat_T
		},

		get_degrees_of_freedom = function(){
			NA_real_
		}
	)
)
