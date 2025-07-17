#' A Sequential Design
#' 
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data intialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#' 
SeqDesignKK21stepwiseabstract = R6::R6Class("SeqDesignKK21stepwiseabstract",
	inherit = SeqDesignKK21abstract,
	public = list(
		#' 
		#' @description
		#' Initialize a sequential experimental design
		#' 
  		#' @param response_type 	The data type of response values which must be one of the following: 
		#' 							"continuous", 
		#' 							"incidence", 
		#' 							"proportion", 
		#' 							"count", 
		#' 							"survival".
		#' 							This package will enforce that all added responses via the \code{add_subject_response} method will be
		#' 							of the appropriate type.
		#' @param prob_T	The probability of the treatment assignment. This defaults to \code{0.5}.
		#' @param include_is_missing_as_a_new_feature	If missing data is present in a variable, should we include another dummy variable for its
		#' 												missingness in addition to imputing its value? If the feature is type factor, instead of creating
		#' 												a new column, we allow missingness to be its own level. The default is \code{TRUE}.
		#' @param n			The sample size (if fixed). Default is \code{NULL} for not fixed.
		#' @param verbose	A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}.
		#' @param num_boot the number of bootstrap samples taken to approximate the subject-distance distribution. Default is \code{NULL} for not 500.
		#' @return 			A new `SeqDesign` object of the specific type
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(response_type = "continuous")
		#'  
		initialize = function(
			response_type, 
			prob_T = 0.5,
			include_is_missing_as_a_new_feature = TRUE, 
			verbose = TRUE,
			n = NULL,
			num_boot = NULL
		){
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, verbose, n, num_boot)
		}
	),
	private = list(
		
		assign_wt = function(){
			wt = 	if (sum(private$match_indic == 0, na.rm = TRUE) == 0 | private$t <= (ncol(private$Xraw) + 2) | private$t <= private$t_0){
						#we're early or the reservoir is empty, so randomize
						#cat("    assign_wt_KK21stepwise CRD t", private$t, "\n")
						private$match_indic[private$t] = 0
						private$assign_wt_CRD()
					} else if (is.null(private$X) | (sum(!is.na(private$y)) < 2 * (ncol(private$X) + 2))){ 
						#This is the number of responses collected before
						#the algorithm begins estimating the covariate-specific weights. If left unspecified this defaults to \code{2 * (p + 2)} i.e. two data points
						#for every glm regression parameter estimated (p covariates, the intercept and the coefficient of the additive treatment effect). The minimum
						#value is p + 2 to allow OLS to estimate. If this number of not met, we default to KK14 matching (which is better than nothing as it matches
						#the covariate distributions at the very least).	
						private$assign_wt_KK14()
					} else {		
						all_subject_data = private$compute_all_subject_data()		
						#1) need to calculate the weights...
						
						#we need to now run the appropriate (based on response_type) stepwise procedure to get the weights
						raw_weights = private[[paste0("compute_weights_KK21stepwise_", private$response_type)]](
							all_subject_data$X_all_with_y_scaled, #to calculate weights, we need to use only the data that has y's!
							all_subject_data$y_all,
							all_subject_data$w_all_with_y_scaled, 
							all_subject_data$dead_all
						)						
						if (any(is.na(raw_weights)) | any(is.infinite(raw_weights)) | any(is.nan(raw_weights)) | any(raw_weights < 0)){
							stop("raw weight values illegal KK21stepwise")
						}
						#ensure the weights are normalized
						private$covariate_weights = raw_weights / sum(raw_weights)
						names(private$covariate_weights) = colnames(all_subject_data$X_all_with_y_scaled)
						#cat("    assign_wt_KK21stepwise using weights t", private$t, "weights", weights, "\n")
						
						#2) now iterate over all items in reservoir and calculate the weighted sqd distiance vs new guy 
						reservoir_indices = which(private$match_indic == 0)
						weighted_features = colnames(all_subject_data$X_all_with_y_scaled)
						x_new = all_subject_data$xt_all_scaled[weighted_features]
						X_all_scaled_col_subset = all_subject_data$X_all_scaled[, weighted_features, drop = FALSE]
						
#						weighted_sqd_distances = array(NA, length(reservoir_indices))
#						for (r in 1 : length(reservoir_indices)){
#							x_r_x_new_delta = x_new - X_all_scaled_col_subset[reservoir_indices[r], ] #we have to use all data here
#							weighted_sqd_distances[r] = x_r_x_new_delta^2 %*% private$covariate_weights			
#						}
						weighted_sqd_distances = compute_weighted_sqd_distances_cpp(
													x_new,
												    X_all_scaled_col_subset,
												    reservoir_indices,
												    private$covariate_weights
												 )						
						#3) find minimum weighted sqd distiance index
						min_weighted_sqd_dist_index = which(weighted_sqd_distances == min(weighted_sqd_distances))
						
						#generate a cutoff for the weighted minimum distance squared based on bootstrap
#						bootstrapped_weighted_sqd_distances = array(NA, private$other_params$num_boot)
#						for (b in 1 : private$other_params$num_boot){
#							two_xs  = X_all_scaled_col_subset[sample.int(private$t, 2), ] #private$X[sample_int_ccrank(private$t, 2, rep(1, (private$t))), ] 
#							delta_x = two_xs[1, ] - two_xs[2, ]
#							bootstrapped_weighted_sqd_distances[b] = delta_x^2 %*% private$covariate_weights
#						}
						
						bootstrapped_weighted_sqd_distances = compute_bootstrapped_weighted_sqd_distances_cpp(
															    X_all_scaled_col_subset,
															    private$covariate_weights,
															    private$t,
															    private$other_params$num_boot
															  )
							
						min_weighted_dsqd_cutoff_sq = quantile(bootstrapped_weighted_sqd_distances, private$other_params$lambda)
						
						#5) Now, does the minimum make the cut?
						if (length(weighted_sqd_distances[min_weighted_sqd_dist_index]) > 1 || length(min_weighted_dsqd_cutoff_sq) > 1){
							min_weighted_sqd_dist_index = min_weighted_sqd_dist_index[1] #if there's a tie, just take the first one
						}
						#  (a) if it's smaller than the threshold, we're in business: match it
						if (weighted_sqd_distances[min_weighted_sqd_dist_index] < min_weighted_dsqd_cutoff_sq){
							new_match_id = max(private$match_indic, na.rm = TRUE) + 1 #the ID of a new match
							private$match_indic[reservoir_indices[min_weighted_sqd_dist_index]] = new_match_id
							private$match_indic[private$t] = new_match_id
							1 - private$w[reservoir_indices[min_weighted_sqd_dist_index]]
						} else { # (b) otherwise, randomize and add it to the reservoir
							private$match_indic[private$t] = 0	
							private$assign_wt_CRD()	
						}
					}
			if (is.na(private$match_indic[private$t])){ #this should never happen
				stop("no match data recorded")
			}
			wt
		},


		
		compute_weights_KK21stepwise = function(Xfull, response_obj, ws, abs_z_compute_fun){
			weights = array(NA, ncol(Xfull))
			j_droppeds = c()
			X_stepwise = matrix(NA, nrow = nrow(Xfull), ncol = 0)
			
			repeat {
				covs_to_try = setdiff(1 : ncol(Xfull), j_droppeds)				
				if (length(covs_to_try) == 0){ #if there's none left, we jet
					break
				}
				abs_approx_zs = array(NA, ncol(Xfull))
				for (j in covs_to_try){
					abs_approx_zs[j] = abs_z_compute_fun(response_obj, cbind(Xfull[, j], X_stepwise, ws))
				}
				j_max = which.max(abs_approx_zs)
				weights[j_max] = abs_approx_zs[j_max]
				j_droppeds = c(j_droppeds, j_max)
				X_stepwise = cbind(X_stepwise, Xfull[, j_max])
			}
			if (any(is.na(weights))){
				stop("boom")					
			}
			weights
		},
		
		compute_weights_KK21stepwise_continuous = function(xs, ys, ws, ...){
			private$compute_weights_KK21stepwise(xs, scale(ys), ws, function(response_obj, covariate_data_matrix){
				
#				ols_mod = lm(response_obj ~ covariate_data_matrix)
#				abs(coef(suppressWarnings(summary(ols_mod)))[2, 3])
				
				#25% SPEEDUP
			    Xmat = cbind(1, covariate_data_matrix)
			    qr_decomp = qr(Xmat)
			    Qmat = qr.Q(qr_decomp)
			    Rmat = qr.R(qr_decomp)
			    Qmatt = t(Qmat)
			    Rmatt = t(Rmat)
			    RtRinv = solve(Rmatt %*% Rmat)
			    b = (RtRinv %*% Rmatt %*% Qmatt %*% response_obj)
			    s_sq_e = sum((response_obj - Xmat %*% b)^2) / (length(response_obj) - ncol(Qmat))
			    abs(b[2]) / sqrt(s_sq_e * RtRinv[2, 2])
			})
		},
		
		compute_weights_KK21stepwise_incidence = function(xs, ys, ws, ...){
			private$compute_weights_KK21stepwise(xs, ys, ws, function(response_obj, covariate_data_matrix){
				logistic_regr_mod = suppressWarnings(glm(response_obj ~ covariate_data_matrix, family = "binomial"))
				abs(coef(summary_glm_lean(logistic_regr_mod))[2, 3])
			})
		},
		
		compute_weights_KK21stepwise_count = function(xs, ys, ws, ...){	
			private$compute_weights_KK21stepwise(xs, ys, ws, function(response_obj, covariate_data_matrix){
				negbin_regr_mod = robust_negbinreg(response_obj ~ ., cbind(data.frame(response_obj = response_obj), covariate_data_matrix))
				abs(coef(summary_glm_lean(negbin_regr_mod))[2, 3])
			})	
		},
		
		compute_weights_KK21stepwise_proportion = function(xs, ys, ws, ...){
			if (!private$other_params$proportion_use_speedup){
				tryCatch({
					weight = 	private$compute_weights_KK21stepwise(xs, ys, ws, function(response_obj, covariate_data_matrix){					
									beta_regr_mod = suppressWarnings(betareg::betareg(response_obj ~ ., data = cbind(data.frame(response_obj = response_obj), covariate_data_matrix)))
									summary_beta_regr_mod = coef(summary(beta_regr_mod)) 
									tab = 	if (!is.null(summary_beta_regr_mod$mean)){ #beta model
												summary_beta_regr_mod$mean 
											} else if (!is.null(summary_beta_regr_mod$mu)){ #extended-support xbetax model
												summary_beta_regr_mod$mu
											}
									ifelse(nrow(tab) >= 2, abs(tab[2, 3]), NA)
								})
					if (!is.na(weight)){
						return(weight)
					}
				}, error = function(e){})
			}
			#if that didn't work, let's just use the continuous weights on a transformed proportion
			ys[ys == 0] = .Machine$double.eps
			ys[ys == 1] = 1 - .Machine$double.eps
			private$compute_weights_KK21stepwise_continuous(xs, log(ys / (1 - ys)), ws, ...)
		},
		
		compute_weights_KK21stepwise_survival = function(xs, ys, ws, deaths){		
			private$compute_weights_KK21stepwise(xs, survival::Surv(ys, deaths), ws, function(response_obj, covariate_data_matrix){
				#sometims the weibull is unstable... so try other distributions... this doesn't matter since we are just trying to get weights
				#and we are not relying on the model assumptions
				for (dist in c("weibull", "lognormal", "loglogistic")){
					surv_regr_mod = robust_survreg_with_surv_object(response_obj, covariate_data_matrix, dist = dist)
					if (is.null(surv_regr_mod)){
						break
					}
					summary_surv_regr_mod = suppressWarnings(summary(surv_regr_mod)$table)
					if (any(is.nan(summary_surv_regr_mod))){
						break
					}
					weight = ifelse(nrow(summary_surv_regr_mod) >= 2, abs(summary_surv_regr_mod[2, 3]), NA)
					#1 - summary(weibull_regr_mod)$table[2, 4]
					if (!is.na(weight)){
						return(weight)
					}
				}	
				#if that didn't work, default to OLS and log the survival times... again... this doesn't matter since we are just trying to get weights
				#and we are not relying on the model assumptions
				ols_mod = lm(log(as.numeric(response_obj)[1 : length(response_obj)]) ~ covariate_data_matrix)
				abs(coef(suppressWarnings(summary(ols_mod)))[2, 3])
			})	
		}		
	)
)