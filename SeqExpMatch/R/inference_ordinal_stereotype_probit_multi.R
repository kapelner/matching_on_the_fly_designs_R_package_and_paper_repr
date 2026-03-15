#' Multivariate Stereotype Probit Inference for Ordinal Responses
#
#' @description
#' Stereotype probit inference for ordinal responses using treatment plus observed covariates.
#'
#' @export
SeqDesignInferenceOrdinalMultiStereotypeProbitRegr = R6::R6Class("SeqDesignInferenceOrdinalMultiStereotypeProbitRegr",
	inherit = SeqDesignInferenceOrdinalUniStereotypeProbitRegr,
	public = list(

		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),

	private = list(
		stereotype_design_matrix = function(){
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
		}
	)
)
