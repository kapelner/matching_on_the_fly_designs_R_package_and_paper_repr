#' Abstract class for all-subject marginal incidence inference in KK designs
#'
#' @keywords internal
InferenceAbstractKKMarginalIncid = R6::R6Class("InferenceAbstractKKMarginalIncid",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param verbose A flag indicating whether messages should be displayed.
		initialize = function(des_obj,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "incidence")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
		}
	),

	private = list(
		get_covariate_names = function(){
			X = private$get_X()
			p = ncol(X)
			x_names = colnames(X)
			if (is.null(x_names)){
				x_names = paste0("x", seq_len(p))
			}
			x_names
		},

		get_cluster_ids = function(){
			des_priv = private$des_obj_priv_int
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec_int = as.integer(m_vec)
			m_vec_int[is.na(m_vec_int)] = 0L

			# Normalize design's m_vec the same way
			des_m = des_priv$m
			if (is.null(des_m)) des_m = rep(NA_integer_, private$n)
			des_m_int = as.integer(des_m)
			des_m_int[is.na(des_m_int)] = 0L

			# Check design-level cache (only when m_vec matches design's m_vec)
			if (!is.null(des_priv$kk_cluster_id) && identical(m_vec_int, des_m_int)){
				return(des_priv$kk_cluster_id)
			}
			# Check inference-level cache (for bootstrap resamples)
			if (!is.null(private$cached_values$kk_cluster_id) &&
				identical(m_vec_int, private$cached_values$kk_cluster_id_m_vec)){
				return(private$cached_values$kk_cluster_id)
			}

			cluster_id = compute_kk_cluster_ids_cpp(m_vec_int)

			# Store at design level if this is the original m_vec
			if (identical(m_vec_int, des_m_int)){
				des_priv$kk_cluster_id = cluster_id
				des_priv$kk_cluster_id_m_vec = m_vec_int
			} else {
				private$cached_values$kk_cluster_id = cluster_id
				private$cached_values$kk_cluster_id_m_vec = m_vec_int
			}
			cluster_id
		}
	)
)
