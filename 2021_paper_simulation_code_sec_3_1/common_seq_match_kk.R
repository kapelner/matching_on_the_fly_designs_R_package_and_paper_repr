#do the standard matching procedure

#initialize the indicator treatment vector
indic_T = array(NA, n) 
#initialize the reservoir
match_indic = array(-1, n) #0 indicates reservoir

min_t0_based_on_matching_pct = round(t_0_matching_pct * n)

#now we're going to go through and do the matching

for (t0 in 1 : n){
  # cat("for\n")
	#if there is nothing in the reservoir, randomize and add it to reservoir
	if (length(match_indic[match_indic == 0]) == 0 || t0 <= min_t0_based_on_matching_pct){
	  
	  # cat("if\n")
		indic_T[t0] = rbinom(1, 1, prob_trt)
		match_indic[t0] = 0
	} else {
	  # cat("else\n")
		#first calculate the threshold we're operating at
		xs_to_date = x_s[1 : t0, ]
		S_xs_inv = solve(var(xs_to_date))
		F_crit =  qf(prob_match_cutoff_lambda, p, t0 - p)
		T_cutoff_sq = p * (n - 1) / (n - p) * F_crit
		#now iterate over all items in reservoir and take the minimum distance x
		reservoir_indices = which(match_indic == 0)
		x_star = x_s[t0, ]
		sqd_distances = array(NA, length(reservoir_indices))
		for (r in 1 : length(reservoir_indices)){
			sqd_distances[r] = 1 / 2 * 
					t(x_star - x_s[reservoir_indices[r], ]) %*%
					S_xs_inv %*%
					(x_star - x_s[reservoir_indices[r], ])			
		}
#		cat(paste("t", t, "sqd_distances", paste(sqd_distances, collapse = ", "), "T_cutoff_sq", T_cutoff_sq, "\n"))
		
		#find minimum distance index
		min_sqd_dist_index = which(sqd_distances == min(sqd_distances))
		if (length(sqd_distances[min_sqd_dist_index]) > 1 || length(T_cutoff_sq) > 1){
			min_sqd_dist_index = min_sqd_dist_index[1] #if there's a tie, just take the first one
		}
		#if it's smaller than the threshold, we're in business: match it
		if (sqd_distances[min_sqd_dist_index] < T_cutoff_sq){
			match_num = max(match_indic) + 1
			match_indic[reservoir_indices[min_sqd_dist_index]] = match_num
			match_indic[t0] = match_num
			indic_T[t0] = 1 - indic_T[reservoir_indices[min_sqd_dist_index]]
		} else { #otherwise, randomize and add it to the reservoir
			indic_T[t0] = rbinom(1, 1, prob_trt)
			match_indic[t0] = 0
		}
	}
}
#create response vector
sys.source(paste("create_response_", response_model, ".R", sep = ""), envir = environment())
sys.source("common_matching_post_processing.R", envir = environment())