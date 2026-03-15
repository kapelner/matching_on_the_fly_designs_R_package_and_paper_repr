#' Multivariate Adjacent-Category Logit Inference for Ordinal Responses
#'
#' @description
#' Adjacent-category logit inference for ordinal responses with treatment and
#' observed covariates entering linearly into the adjacent-category logits.
#'
#' @export
SeqDesignInferenceOrdinalMultiAdjCatLogitRegr = R6::R6Class("SeqDesignInferenceOrdinalMultiAdjCatLogitRegr",
	inherit = SeqDesignInferenceOrdinalUniAdjCatLogitRegr,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param	seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param	num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference
		#' 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs),
		#' 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
		#' @param	verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),

	private = list(
		adjacent_category_design_matrix = function(){
			X_full = cbind(private$w, private$get_X())
			qr_X = qr(X_full)
			if (qr_X$rank < ncol(X_full)){
				keep = qr_X$pivot[seq_len(qr_X$rank)]
				if (!(1L %in% keep)) keep[qr_X$rank] = 1L
				keep = sort(keep)
				X_full = X_full[, keep, drop = FALSE]
			}
			colnames(X_full)[1] = "treatment"
			X_full
		},

		generate_mod = function(){
			res = fast_adjacent_category_logit_with_var_cpp(
				X = private$adjacent_category_design_matrix(),
				y = as.numeric(private$y)
			)
			list(
				b = c(NA, res$b),
				ssq_b_2 = res$ssq_b_1
			)
		}
	)
)
