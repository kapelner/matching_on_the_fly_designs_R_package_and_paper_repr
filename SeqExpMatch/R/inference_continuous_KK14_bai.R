#' Inference based on Maximum Likelihood for KK designs
#'
#' @description
#' Inference for mean difference
#'
#'
#' @export
InferenceBaiAdjustedTKK14 = R6::R6Class("InferenceBaiAdjustedTKK14",
	inherit = InferenceBaiAdjustedT,
	public = list(

	),

	private = list(
	distance = function(avg1, avg2){
		sum((avg1 - avg2)^2)
	}
	)
)
