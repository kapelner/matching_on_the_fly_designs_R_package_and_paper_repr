#' Inference based on Maximum Likelihood for KK designs  
#'
#' @description
#' Inference for mean difference
#' 
#'
#' @export
SeqDesignInferenceBaiAdjustedTKK21 = R6::R6Class("SeqDesignInferenceBaiAdjustedTKK21",
  inherit = SeqDesignInferenceBaiAdjustedT,
  public = list(
    
    #' @description
    #' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
    #' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
    #' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference 
    #' 							(which is very slow). The default is 1 for serial computation. This parameter is ignored
    #' 							for \code{test_type = "MLE-or-KM-based"}.
    #' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
    #' @param convex      A flag indicating whether the estimator should use a convex combination of the Bai et all
    #'                    matched pairs estimate with the reservoir estimate, or just the Bai et all estimate on by its self.
    #' 
    initialize = function(seq_des_obj, num_cores = 1, verbose = TRUE, convex = FALSE){
      super$initialize(seq_des_obj, num_cores, verbose, convex = convex)
      assertNoCensoring(private$any_censoring)
      assert_class(seq_des_obj, "SeqDesignKK21")
    }
  ),
  
  private = list(
    distance = function(avg1, avg2){
      sum(private$seq_des_obj_priv_int$covariate_weights * (avg1 - avg2)^2)
    }
  )		
)