#' A Sequential Design
#' 
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data intialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#' 
#' @export
SeqDesignKK21 = R6::R6Class("SeqDesignKK21",
	inherit = SeqDesignKK14,
	public = list(
		#' 
		#' @description
		#' Initialize a matching-on-the-fly sequential experimental design which matches based on Kapelner and Krieger (2021) with
		#' option to use matching parameters of Morrison and Owen (2025)
		#'
		#' @param response_type 	The data type of response values which must be one of the following: 
		#' 							"continuous"(the default), 
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
		#' @param lambda   The quantile cutoff of the subject distance distribution for determining matches. If unspecified and \code{morrison = FALSE}, default is 10\%.
		#' @param t_0_pct  The percentage of total sample size n where matching begins. If unspecified and \code{morrison = FALSE}, default is 35\%.
		#' @param morrison 	Default is \code{FALSE} which implies matching via the KK14 algorithm using \code{lambda} and \code{t_0_pct} matching.
		#'					If \code{TRUE}, we use Morrison and Owen (2025)'s formula for \code{lambda} which differs in the fixed n versus variable n
		#'					settings and matching begins immediately with no wait for a certain reservoir size like in KK14.
		#' @param p			The number of covariate features. Must be specified when \code{morrison = TRUE} otherwise do not specify this argument.
		#' @param num_boot the number of bootstrap samples taken to approximate the subject-distance distribution. Default is \code{NULL} for not 500.
		#' @param proportion_use_speedup 	Should we speed up the estimation of the weights in the response = proportion case via a continuous regression on log(y/(1-y))
		#' 									instead of a beta regression each time? This is at the expense of the weights being less accurate. Default is \code{TRUE}.
		#' @param thin		For internal use only. Do not specify. You can thank R6's single constructor-only for this coding noise.

		#' @return 			A new `SeqDesignKK21` object
		#' 
		#' @examples
		#' seq_des = SeqDesignKK21$new(response_type = "continuous")
		#'  
		initialize = function(
			response_type = "continuous",  
			prob_T = 0.5,
			include_is_missing_as_a_new_feature = TRUE, 
			n = NULL,
			verbose = FALSE,
			lambda = NULL,
			t_0_pct = NULL,
			morrison = FALSE,
			p = NULL,
			num_boot = NULL,
			proportion_use_speedup = TRUE,
			thin = FALSE
		){
			if (!thin){
				super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, lambda, t_0_pct, morrison, p)
				if (is.null(num_boot)){
					num_boot = 500
				} else {
					assertCount(num_boot, positive = TRUE)
				}
				private$num_boot = num_boot
				assertFlag(proportion_use_speedup)	
				private$proportion_use_speedup = proportion_use_speedup						
			}
		}
	),
	private = list(
		covariate_weights = NULL,
		num_boot = NULL,
		proportion_use_speedup = NULL,
		
		duplicate = function(){
			d = super$duplicate()
			d$.__enclos_env__$private$covariate_weights = private$covariate_weights
			d$.__enclos_env__$private$num_boot = private$num_boot
			d$.__enclos_env__$private$proportion_use_speedup = private$proportion_use_speedup
			d
		},
		
		assign_wt = function(){
			wt = 	if (private$too_early_to_match()){
						#we're early or the reservoir is empty, so randomize
						#cat("    assign_wt", class(self)[1], " t", private$t, "\n")
						private$match_indic[private$t] = 0
						private$assign_wt_CRD()
					} else if (is.null(private$X) | (sum(!is.na(private$y)) < 2 * (ncol(private$X) + 2))){ 
						#This is the number of responses collected before
						#the algorithm begins estimating the covariate-specific weights. If left unspecified this defaults to \code{2 * (p + 2)} i.e. two data points
						#for every glm regression parameter estimated (p covariates, the intercept and the coefficient of the additive treatment effect). The minimum
						#value is p + 2 to allow OLS to estimate. If this number of not met, we default to KK14 matching (which is better than nothing as it matches
						#the covariate distributions at the very least).	
						super$assign_wt()
					} else {
						all_subject_data = private$compute_all_subject_data()
						#1) need to calculate the weights - this is different between KK21 and KK21stepwise
						raw_weights = private$compute_weights(all_subject_data)						

						if (any(is.na(raw_weights)) | any(is.infinite(raw_weights)) | any(is.nan(raw_weights)) | any(raw_weights < 0)){
							stop("raw weight values illegal in design ", class(self)[1])
						}
						#ensure the weights are normalized
						private$covariate_weights = raw_weights / sum(raw_weights)
						names(private$covariate_weights) = colnames(all_subject_data$X_all_with_y_scaled)
						#cat("    assign_wt_KK21 using sorted weights t", private$t, "weights", sort(weights), "\n")
						
						#2) now iterate over all items in reservoir and calculate the weighted sqd distiance vs new guy 
						reservoir_indices = which(private$match_indic == 0)
						weighted_features = colnames(all_subject_data$X_all_with_y_scaled)
						x_new = all_subject_data$xt_all_scaled[weighted_features]
						X_all_scaled_col_subset = all_subject_data$X_all_scaled[, weighted_features, drop = FALSE]
						
#						weighted_sqd_distances = array(NA, length(reservoir_indices))
#						for (r in 1 : length(reservoir_indices)){
#							x_r_x_new_delta = x_new - X_all_scaled_col_subset[reservoir_indices[r], ]
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
#						bootstrapped_weighted_sqd_distances = array(NA, private$num_boot)
#						for (b in 1 : private$num_boot){
#							two_xs  = X_all_scaled_col_subset[sample.int(private$t, 2), ] #private$X[sample_int_ccrank(private$t, 2, rep(1, (private$t))), ]
#							delta_x = two_xs[1, ] - two_xs[2, ]
#							bootstrapped_weighted_sqd_distances[b] = delta_x^2 %*% private$covariate_weights
#						}

						bootstrapped_weighted_sqd_distances = compute_bootstrapped_weighted_sqd_distances_cpp(
															    X_all_scaled_col_subset,
															    private$covariate_weights,
															    private$t,
															    private$num_boot
															  )

						min_weighted_dsqd_cutoff_sq = quantile(bootstrapped_weighted_sqd_distances, private$compute_lambda())
						
						#5) Now, does the minimum make the cut?
						if (length(weighted_sqd_distances[min_weighted_sqd_dist_index]) > 1 || length(min_weighted_dsqd_cutoff_sq) > 1){
							min_weighted_sqd_dist_index = min_weighted_sqd_dist_index[1] #if there's a tie, just take the first one
						}
						#  (a) if it's smaller than the threshold, we're in business: match it
						if (weighted_sqd_distances[min_weighted_sqd_dist_index] < min_weighted_dsqd_cutoff_sq){
							new_match_id = max(private$match_indic, na.rm = TRUE) + 1 #the ID of a new match
							private$match_indic[reservoir_indices[min_weighted_sqd_dist_index]] = new_match_id
							private$match_indic[private$t] = new_match_id
							#assign opposite
							1 - private$w[reservoir_indices[min_weighted_sqd_dist_index]]
						# (b) otherwise, randomize and add it to the reservoir
						} else { 
							private$match_indic[private$t] = 0	
							private$assign_wt_CRD()
						}
					}
			if (is.na(private$match_indic[private$t])){ #this should never happen
				stop("no match data recorded")
			}
			wt
		},
		
		compute_weights = function(all_subject_data){
			raw_weights = array(NA, all_subject_data$rank_all_with_y_scaled)
			for (j in 1 : all_subject_data$rank_all_with_y_scaled){
				raw_weights[j] = private[[paste0("compute_weight_KK21_", private$response_type)]](
					all_subject_data$X_all_with_y_scaled, 
					all_subject_data$y_all, 
					all_subject_data$dead_all, 
					j
				) 
			}
			raw_weights
		},
		
		compute_weight_KK21_continuous = function(xs_to_date, ys_to_date, deaths_to_date, j){
			if (nrow(xs_to_date) == 1){
				.Machine$double.eps
			} else {
				mod = fast_ols_with_var_cpp(cbind(1, xs_to_date[, j, drop = FALSE]), ys_to_date)
				abs(mod$b[2] / sqrt(mod$ssq_b_j))
			}
#			ols_mod = lm(ys_to_date ~ xs_to_date[, j])
#			summary_ols_mod = suppressWarnings(coef(summary(ols_mod)))
#			
#			
#			#1 - coef(summary(logistic_regr_mod))[2, 4]
#			#if there was only one row, then this feature was all one unique value... so send back a weight of nada
#			ifelse(nrow(summary_ols_mod) >= 2, abs(summary_ols_mod[2, 3]), )
		},
		
		compute_weight_KK21_incidence = function(xs_to_date, ys_to_date, deaths_to_date, j){
			if (nrow(xs_to_date) == 1){
				.Machine$double.eps
			} else {
				mod = fast_logistic_regression_with_var(cbind(1, xs_to_date[, j, drop = FALSE]), ys_to_date)
				abs(mod$b[2] / sqrt(mod$ssq_b_2))
			}			
#			tryCatch({
#				logistic_regr_mod = suppressWarnings(glm(ys_to_date ~ xs_to_date[, j], family = "binomial"))
#			}, error = function(e){return(.Machine$double.eps)}) #sometimes these glm's blow up and we don't really care that much
#			summary_logistic_regr_mod = coef(summary_glm_lean(logistic_regr_mod))
#			#1 - coef(summary(logistic_regr_mod))[2, 4]
#			#if there was only one row, then this feature was all one unique value... so send back a weight of nada
#			ifelse(nrow(summary_logistic_regr_mod) >= 2, abs(summary_logistic_regr_mod[2, 3]), .Machine$double.eps)
		},
		
		compute_weight_KK21_count = function(xs_to_date, ys_to_date, deaths_to_date, j){
			tryCatch({
				negbin_regr_mod = suppressWarnings(MASS::glm.nb(y ~ x, data = data.frame(x = xs_to_date[, j], y = ys_to_date)))
				summary_negbin_regr_mod = coef(summary_glm_lean(negbin_regr_mod))
				#1 - coef(summary(negbin_regr_mod))[2, 4]
				#if there was only one row, then this feature was all one unique value... so send back a weight of nada
				return(ifelse(nrow(summary_negbin_regr_mod) >= 2, abs(summary_negbin_regr_mod[2, 3]), .Machine$double.eps))
			}, error = function(e){}) #sometimes these glm's blow up and we don't really care that much
			.Machine$double.eps #otherwise weight it nada
		},
		
		compute_weight_KK21_proportion = function(xs_to_date, ys_to_date, deaths_to_date, j){
			if (!private$proportion_use_speedup){
				tryCatch({
					beta_regr_mod = suppressWarnings(betareg::betareg(y ~ x, data = data.frame(x = xs_to_date[, j], y = ys_to_date)))
					summary_beta_regr_mod = coef(summary(beta_regr_mod)) 
					#1 - coef(summary(beta_regr_mod))$mean[2, 4]
					tab = 	if (!is.null(summary_beta_regr_mod$mean)){ #beta model
								summary_beta_regr_mod$mean 
							} else if (!is.null(summary_beta_regr_mod$mu)){ #extended-support xbetax model
								summary_beta_regr_mod$mu
							}
					if (nrow(tab) >= 2){
						return(abs(tab[2, 3]))
					}
				}, error = function(e){}) #sometimes these glm's blow up and we don't really care that much
			}
			#if that didn't work, let's just use the continuous weights on a transformed proportion
			ys_to_date[ys_to_date == 0] = .Machine$double.eps
			ys_to_date[ys_to_date == 1] = 1 - .Machine$double.eps
			private$compute_weight_KK21_continuous(xs_to_date, log(ys_to_date / (1 - ys_to_date)), deaths_to_date, j)
		},
		
		compute_weight_KK21_survival = function(xs_to_date, ys_to_date, deaths_to_date, j){
			surv_obj = survival::Surv(ys_to_date, deaths_to_date)
			#sometims the weibull is unstable... so try other distributions... this doesn't matter since we are just trying to get weights
			#and we are not relying on the model assumptions
			for (dist in c("weibull", "lognormal", "loglogistic")){
				surv_regr_mod = robust_survreg_with_surv_object(surv_obj, xs_to_date[, j], dist = dist)
				if (is.null(surv_regr_mod)){
					break
				}
				summary_surv_regr_mod = suppressWarnings(summary(surv_regr_mod)$table)
				if (any(is.nan(summary_surv_regr_mod))){
					break
				}
				weight = ifelse(nrow(summary_surv_regr_mod) >= 2, abs(summary_surv_regr_mod[2, 3]), NA)
				#1 - summary(weibull_regr_mod)$table[2, 4]
				if (!is.na(weight)){ #sometimes these glm's blow up and we don't really care that much
					return(weight)
				}
			}
			#if that didn't work, default to OLS and log the survival times... again... this doesn't matter since we are just trying to get weights
			#and we are not relying on the model assumptions being true
			private$compute_weight_KK21_continuous(xs_to_date, log(ys_to_date), deaths_to_date, j)
		}	
		
	)
)