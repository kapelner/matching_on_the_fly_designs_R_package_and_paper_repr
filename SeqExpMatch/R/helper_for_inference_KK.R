#' Helper class for KK designs
#' 
#' @description
#' An abstract R6 Class that aids in tests and confidence intervals all KK sequential designs.
#' 
SeqDesignInferenceHelperKK = R6::R6Class("SeqDesignInferenceHelperKK",
	public = list(
		#' Initialize a helper object
		#' @description
		#' Initialize a helper object
		#'		
		#' @param seq_des_obj	A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param X				A data.table object of the subjects' covariates in the study
		initialize = function(seq_des_obj, X){
			if (!seq_des_obj$.__enclos_env__$private$isKK){
				stop("For this type of inference, you need a KK design.")
			}
			private$seq_des_obj = seq_des_obj
			private$X = X
			private$match_indic = seq_des_obj$.__enclos_env__$private$match_indic
		},
				
		compute_continuous_post_matching_data_KK = function(){	
			if (is.null(private$KKstats)){
				#get matched data
				
				m = max(private$match_indic)
				y_matched_diffs = array(NA, m)
				X_matched_diffs = matrix(NA, nrow = m, ncol = ncol(private$X))
				if (m > 0){
					for (match_id in 1 : m){ #we want to just calculate the diffs inside matches and ignore the reservoir
						yTs = private$seq_des_obj$y[private$seq_des_obj$w == 1 & private$match_indic == match_id]
						yCs = private$seq_des_obj$y[private$seq_des_obj$w == 0 & private$match_indic == match_id]
						y_matched_diffs[match_id] = mean(yTs) - mean(yCs)
						
						xmTs = private$X[private$seq_des_obj$w == 1 & private$match_indic == match_id, ]
						xmCs = private$X[private$seq_des_obj$w == 0 & private$match_indic == match_id, ]
						X_matched_diffs[match_id, ] = mean(xmTs) - mean(xmCs)
					}
				}			
				
				#get reservoir data
				X_reservoir = 	private$X[private$match_indic == 0, ]
				y_reservoir = 	private$seq_des_obj$y[private$match_indic == 0]
				w_reservoir = 	private$seq_des_obj$w[private$match_indic == 0]
				y_reservoir_T = y_reservoir[w_reservoir == 1] #get the reservoir responses from the treatment
				y_reservoir_C = y_reservoir[w_reservoir == 0] #get the reservoir responses from the control
				r_bar = mean(y_reservoir_T) - mean(y_reservoir_C) #compute the classic estimator from the reservoir: ybar_T - ybar_C
				
				#get reservoir sample sizes
				nRT = length(y_reservoir_T) #how many treatment observations are there in the reservoir?
				nRC = length(y_reservoir_C) #how many control observations are there in the reservoir?
				nR = nRT + nRC #how many observations are there in the reservoir?
				
				ssqR = (var(y_reservoir_T) * (nRT - 1) + var(y_reservoir_C) * (nRC - 1)) / (nR - 2) * (1 / nRT + 1 / nRC)
				ssqD_bar = var(y_matched_diffs) / m
				
				w_star = ssqR / (ssqR + ssqD_bar)
				
				private$KKstats = list(
					X_matched_diffs = X_matched_diffs,
					y_matched_diffs = y_matched_diffs,
					X_reservoir = X_reservoir,
					y_reservoir = y_reservoir,
					w_reservoir = w_reservoir,
					nRT = nRT,
					nRC = nRC,
					m = m,
					d_bar = mean(y_matched_diffs),
					ssqD_bar = ssqD_bar,
					r_bar = r_bar,
					ssqR = ssqR,
					w_star = w_star
				)
			}
			private$KKstats
		}
	),
	private = list(	
		seq_des_obj = NULL,
		X = NULL,
		KKstats = NULL,	
		match_indic = NULL,
		match_indic_model_matrix = NULL,
		
		generate_data_frame_with_matching_dummies = function(){
			if (!is.null(private$match_indic_model_matrix)){
				if (!is.null(private$match_indic) & uniqueN(private$match_indic) > 1){
					mm = model.matrix(~ 0 + factor(private$match_indic)) 
					mm = mm[, 2 : (ncol(mm) - 1)]
					private$match_indic_model_matrix = mm
				}
			}
			cbind(data.frame(w = private$seq_des_obj$w), private$X, private$match_indic_model_matrix)
		}
		
	)
)
		