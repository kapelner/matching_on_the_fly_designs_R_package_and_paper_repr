#' Abstract class for all-subject marginal incidence inference in KK designs
#'
#' @keywords internal
InferenceAbstractKKMarginalIncid = R6::R6Class("InferenceAbstractKKMarginalIncid",
	inherit = InferenceKKPassThrough,
	public = list(

		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "incidence")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose)
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
			cluster_id = private$cached_values$cluster_id
			if (!is.null(cluster_id)) return(cluster_id)

			m_vec = private$m
			if (is.null(m_vec)){
				m_vec = rep(0L, private$n)
			}
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L
			cluster_id = compute_kk_cluster_ids_cpp(m_vec)
			private$cached_values$cluster_id = cluster_id
			cluster_id
		}
	)
)
