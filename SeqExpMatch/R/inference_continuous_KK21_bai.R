#' Inference based on Maximum Likelihood for KK designs
#'
#' @description
#' Inference for mean difference
#'
#'
#' @export
InferenceBaiAdjustedTKK21 = R6::R6Class("InferenceBaiAdjustedTKK21",
	inherit = InferenceBaiAdjustedT,
	public = list(

	),

	private = list(
	distance = function(avg1, avg2){
		sum(private$des_obj_priv_int$covariate_weights * (avg1 - avg2)^2)
	}
	)
)
