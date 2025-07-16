#' 					"KK21" for Kapelner and Krieger's (2021) CARA Matching on the Fly with Differential Covariate Weights Design
#' 					"KK21stepwise" for Kapelner and Krieger's (2021) CARA Matching on the Fly with Differential Covariate Weights Stepwise Design


#' 					All "KK21" designs further require "num_boot" which is the number of bootstrap samples taken to approximate the subject-distance distribution. 
#' 					If unspecified, default is 500. There is an optional flag "proportion_use_speedup = TRUE" which uses a continuous regression on log(y/(1-y))
#' 					instead of a beta regression each time to generate the weights in KK21 designs. The default is this flag is on.

	if (grepl("KK21", design)){
		private$isKK21 = TRUE
		if (is.null(private$other_params$num_boot)){
			private$other_params$num_boot = 500
		} else {
			assertCount(private$other_params$num_boot, positive = TRUE)
		}
		if (is.null(private$other_params$proportion_use_speedup)){
			private$other_params$proportion_use_speedup = TRUE							
		} else {
			assertFlag(private$other_params$proportion_use_speedup)	
		}
	}
}