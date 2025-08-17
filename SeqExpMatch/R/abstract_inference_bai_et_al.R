#' Inference based on Maximum Likelihood for KK designs  
#'
#' @description
#' Inference for mean difference
#' 
#'
SeqDesignInferenceBaiAdjustedT = R6::R6Class("SeqDesignInferenceBaiAdjustedT",
  inherit = SeqDesignInferenceMLEorKMKK,
  public = list(
    
    #' @description
    #' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
    #' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
    #' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference 
    #' 							(which is very slow). The default is 1 for serial computation. This parameter is ignored
    #' 							for \code{test_type = "MLE-or-KM-based"}.
    #' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
    #' @param convex      A flag indicating whether the estimator should use a convex combination of the Bai et al
    #'                    matched pairs estimate with the reservoir estimate, or just the Bai et al estimate by its self.
	#' @param thin		For internal use only. Do not specify. You can thank R6's single constructor-only for this coding noise.
    #' 
    initialize = function(seq_des_obj, num_cores = 1, verbose = TRUE, convex = FALSE, thin = FALSE){
		if (!thin){
	      super$initialize(seq_des_obj, num_cores, verbose)
	      private$convex = convex
		  private$compute_reservoir_and_match_statistics()
	      assertNoCensoring(private$any_censoring)
		}
    },
    
    #' Compute treatment effect
    #'	
    #' @description
    #' Computes the appropriate estimate for compound mean difference across pairs and reservoir
    #' 
    #' @return 	The setting-appropriate (see description) numeric estimate of the treatment effect
    #' 
    #' @examples
    #' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD", response_type = "continuous")
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
    #' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
    #' 
    #' seq_des_inf = SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des)
    #' seq_des_inf$compute_treatment_estimate()
    #' 	
    compute_treatment_estimate = function(){
      if (!private$convex || private$KKstats$nRT <= 1 || private$KKstats$nRC <= 1){ #if er are not using the res in the test, only use the match pairs
        private$cached_values$beta_hat_T = private$KKstats$d_bar	
      } else if (private$KKstats$m == 0){ #sometimes there's no matches
        private$cached_values$beta_hat_T = private$KKstats$r_bar			
      } else {
        if (is.null(private$cached_values$s_beta_hat_T)){
          private$shared()
        }
        w_star_bai = private$KKstats$ssqR / (private$KKstats$ssqR + private$bai_var_d_bar)
        private$cached_values$beta_hat_T = w_star_bai * private$KKstats$d_bar + (1 - w_star_bai) * private$KKstats$r_bar #proper weighting
      }
      private$cached_values$beta_hat_T
    },
    
    #' Compute confidence interval
    #'
    #' @description
    #' Computes a 1-alpha level frequentist confidence interval
    #' 
    #' Here we use the theory that MLE's computed for GLM's are asymptotically normal (except in the case 
    #' of estimat_type "median difference" where a nonparametric bootstrap confidence interval (see the \code{controlTest::quantileControlTest} method)
    #' is employed. Hence these confidence intervals are asymptotically valid and thus approximate for any sample size.
    #' 
    #' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
    #' 
    #' @return 	A (1 - alpha)-sized frequentist confidence interval for the treatment effect
    #' 
    #' @examples
    #' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD")
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
    #' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
    #' 
    #' seq_des_inf = SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des, test_type = "MLE-or-KM-based")
    #' seq_des_inf$compute_confidence_interval()
    #'	
    compute_mle_confidence_interval = function(alpha = 0.05){
      assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
      if (is.null(private$cached_values$beta_hat_T)){
        private$compute_treatment_estimate()
      }
      if (is.null(private$cached_values$s_beta_hat_T)){
        private$shared()
      }		
      private$cached_values$is_z = TRUE
      private$compute_z_or_t_ci_from_s_and_df(alpha)
    },
    
    #' Compute p-value
    #'
    #' @description
    #' Computes a 2-sided p-value
    #'
    #' @param delta	The null difference to test against. For any treatment effect at all this is set to zero (the default).
    #' 
    #' @return 	The approximate frequentist p-value
    #' 
    #' @examples
    #' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD")
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
    #' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
    #' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
    #' 
    #' seq_des_inf = SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des)
    #' seq_des_inf$compute_two_sided_pval_for_treatment_effect()
    #' 		
    compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
      assertNumeric(delta)
      if (is.null(private$cached_values$beta_hat_T)){
        private$compute_treatment_estimate()
      }
      if (is.null(private$cached_values$s_beta_hat_T)){
        private$shared()
      }
      2 * pnorm(
       -abs(private$cached_values$beta_hat_T / private$cached_values$s_beta_hat_T)
      ) #approximate by using N(0, 1) distribution
      
    }
  ),
  
  private = list(
    convex_flag = NULL,
    bai_var_d_bar = NULL,
    
    shared = function(){
      m = private$KKstats$m
      if (m == 0){
        private$cached_values$s_beta_hat_T = ifelse(private$convex, sqrt(private$KKstats$ssqR), 0)
      } else {
	      private$bai_var_d_bar = private$compute_bai_variance_for_pairs() / m
	      private$cached_values$s_beta_hat_T = if (private$convex_flag && private$KKstats$nRT > 1 && private$KKstats$nRC > 1){
	                                              sqrt(
													private$bai_var_d_bar * private$KKstats$ssqR /
	                                                  (private$bai_var_d_bar + private$KKstats$ssqR)
	                                              ) # convex estimator
	                                            } else {
	                                              sqrt(private$bai_var_d_bar) # just bai estimator
	                                            }		
	  }
    },
    
    compute_bai_variance_for_pairs = function(){  
      pairs_df = data.frame(
        pair_id = 1 : private$KKstats$m,
        yT = private$KKstats$yTs_matched,
        yC = private$KKstats$yCs_matched,
        d_i = private$KKstats$y_matched_diff
      )
      
      halves = private$compute_halves()
      
      delta_sq = mean(pairs_df$d_i)^2
      tau_sq = mean(pairs_df$d_i^2)
      
      # lambda_squ^2 term
      lambda_squ = 0
      if (nrow(halves) > 0){
        for (i in 1 : nrow(halves)){
          pair1_id = as.numeric(as.character(halves[i, 1]))
          pair2_id = as.numeric(as.character(halves[i, 3]))
          
          d1 = pairs_df$d_i[pairs_df$pair_id == pair1_id]
          d2 = pairs_df$d_i[pairs_df$pair_id == pair2_id]
          
          lambda_squ = lambda_squ + (d1 * d2)
        }
        lambda_squ = lambda_squ / nrow(halves)
      }
      v_sq = tau_sq - (lambda_squ + delta_sq) / 2
      
      ########################## alternative bai estimate, does not work ####################
      # if(private$aultimate_bai) {
      #   # This implements the estimator from Remark 3.11, Equation (28) in Bai et al. (2022)
      #   muT = mean(pairs_df$yT)
      #   muC = mean(pairs_df$yC)
      #   
      #   # sigma_hat^2(1) and sigma_hat^2(0)
      #   yT_sig_sq = mean((pairs_df$yT - muT)^2)
      #   yC_sig_sq = mean((pairs_df$yC - muC)^2)
      #   
      #   # lambda_hat^2 term
      #   lambda_squ = 0
      #   if(nrow(halves) > 0){
      #     for (i in 1:nrow(halves)){
      #       pair1_id = as.numeric(as.character(halves[i, 1]))
      #       pair2_id = as.numeric(as.character(halves[i, 3]))
      #       
      #       sum1 = pairs_df$yT[pairs_df$pair_id == pair1_id] + pairs_df$yC[pairs_df$pair_id == pair1_id]
      #       sum2 = pairs_df$yT[pairs_df$pair_id == pair2_id] + pairs_df$yC[pairs_df$pair_id == pair2_id]
      #       
      #       lambda_squ = lambda_squ + (sum1 * sum2)
      #     }
      #     # Note: The sum is over floor(m/2) pairs of pairs. The multiplier is 1/floor(m/2).
      #     # Your (2/m) is a good approximation. Let's make it exact.
      #     lambda_squ = lambda_squ / nrow(halves)
      #   }
      #   
      #   # v^2 = sigma_hat^2(1) + sigma_hat^2(0) - 0.5 * (lambda_hat^2 - (mu_hat(1) + mu_hat(0))^2)
      #   v_sq = yT_sig_sq + yC_sig_sq - 0.5 * (lambda_squ - (muT + muC)^2)
      #   
      # } else {
      ######################################################################################
      
      # The variance cannot be negative.
      max(v_sq, 1e-8)
    },
    
    compute_halves = function(){
      m = private$KKstats$m
      if (m < 2) return(data.frame()) # Cannot make pairs of pairs if there's < 2 pairs
      
      X = private$get_X()
      
      # Compute the average covariate vector for each pair ## can be rcpp'ed for sure
      pair_avg = do.call(rbind, lapply(split(1:nrow(X), private$seq_des_obj_priv_int$match_indic), 
                                       function(i) colMeans(X[i, , drop = FALSE])))
      
      # Create a distance matrix between the pair-average covariates
      dist_mat = matrix(data = 0, nrow = m, ncol = m)
      dist_mat = dist_mat + diag(Inf, nrow = m, ncol = m) #set diag equal to inf
      for(i in 1:(m-1)){
        for(j in (i+1):m){
          d = private$distance(pair_avg[i,], pair_avg[j,]) # distance formulas are defined in the daughter classes
          dist_mat[i,j] = d
          dist_mat[j,i] = d
        }
      }
      
      # Use nbpMatching to find the optimal pairing of the pairs to minimize total distance
      dist_obj = suppressWarnings(nbpMatching::distancematrix(dist_mat))
      match_obj = suppressWarnings(nbpMatching::nonbimatch(dist_obj))
      
      halves = match_obj$halves
      # If there's an odd number of pairs, remove the "ghost" match
      if (m %% 2 == 1){
        ghost_row = which(halves[,3] == "ghost")
        if(length(ghost_row) > 0) halves = halves[-ghost_row, ]
      }
      halves
    }
  )		
)