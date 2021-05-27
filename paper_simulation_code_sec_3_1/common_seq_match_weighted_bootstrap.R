NUM_BOOT_DSQD_DIST = 500

#initialize the indicator treatment vector
indic_T = array(NA, n) 
#initialize the reservoir
match_indic = array(-1, n) #0 indicates reservoir

min_t0_based_on_matching_pct = round(t_0_matching_pct * n)
#now we're going to go through and do the matching

for (t in 1 : n){
	#if there is nothing in the reservoir, randomize and add it to reservoir
	if (length(match_indic[match_indic == 0]) == 0 || t <= p || t <= min_t0_based_on_matching_pct){
		indic_T[t] = rbinom(1, 1, prob_trt)
		match_indic[t] = 0
	} else {
		#1) need to calculate the weights
		# (a) get old x's and y's (we assume that x's are orthogonal for now)
		xs_to_date = x_s[1 : (t - 1), ]
		ys_to_date = y[1 : (t - 1)]
		# (b) run simple correlations to get Rsq's and use them as the relative weights
		weights = array(NA, p)
		for (j in 1 : p){
			weights[j] = cor(xs_to_date[, j], ys_to_date)^2
		}
		# (c) now we need to scale them
		weights = weights / sum(weights)
		# cat("nsim", nsim, "t", t, "weights", weights, "\n")

		#2) now iterate over all items in reservoir and calculate the weighted sqd distiance vs new guy 
		reservoir_indices = which(match_indic == 0)
		x_star = x_s[t, ]
		weighted_sqd_distances = array(NA, length(reservoir_indices))
		for (r in 1 : length(reservoir_indices)){
			delta_x = x_star - x_s[reservoir_indices[r], ]
			weighted_sqd_distances[r] = delta_x^2 %*% weights			
		}
		#3) find minimum weighted sqd distiance index
		min_weighted_sqd_dist_index = which(weighted_sqd_distances == min(weighted_sqd_distances))
		
		#generate a cutoff for the weighted minimum distance squared based on bootstrap
		bootstrapped_weighted_sqd_distances = array(NA, NUM_BOOT_DSQD_DIST)
		for (b in 1 : NUM_BOOT_DSQD_DIST){
			two_xs  = x_s[sample.int(t, 2), ]
			delta_x = two_xs[1, ] - two_xs[2, ]
			bootstrapped_weighted_sqd_distances[b] = delta_x^2 %*% weights
		}
		
		min_weighted_dsqd_cutoff_sq = quantile(bootstrapped_weighted_sqd_distances, prob_match_cutoff_lambda)
		
		#5) Now, does the minimum make the cut?
		if (length(weighted_sqd_distances[min_weighted_sqd_dist_index]) > 1 || length(min_weighted_dsqd_cutoff_sq) > 1){
			min_weighted_sqd_dist_index = min_weighted_sqd_dist_index[1] #if there's a tie, just take the first one
		}
		#  (a) if it's smaller than the threshold, we're in business: match it
		if (weighted_sqd_distances[min_weighted_sqd_dist_index] < min_weighted_dsqd_cutoff_sq){
			match_num = max(match_indic) + 1
			match_indic[reservoir_indices[min_weighted_sqd_dist_index]] = match_num
			match_indic[t] = match_num
			indic_T[t] = 1 - indic_T[reservoir_indices[min_weighted_sqd_dist_index]]
		# (b) otherwise, randomize and add it to the reservoir
		} else { 
			indic_T[t] = rbinom(1, 1, prob_trt)
			match_indic[t] = 0		
		}
	}
	
	#realize the y for subject t
	sys.source(paste("create_response_", response_model, ".R", sep = ""), envir = environment())
}

# cat("nsim", nsim, "t", t, "weights", weights, "\n")
sys.source("common_matching_post_processing.R", envir = environment())