if (nRT <= 1 || nRC <= 1){
	beta_hat_T = d_bar
	pct_only_matches = 1
	pct_only_reservoir = 0	
	yMT = Xymatched[Xymatched$indic_T == 1, ]$y
	yCT = Xymatched[Xymatched$indic_T == 0, ]$y
	matched_correlation = cor(yMT, yCT)
} else if (m == 0){
	beta_hat_T = r_bar
	pct_only_matches = 0
	pct_only_reservoir = 1	
} else {
	beta_hat_T = w_star * d_bar + (1 - w_star) * r_bar
	pct_only_matches = 0
	pct_only_reservoir = 0	
	yMT = Xymatched[Xymatched$indic_T == 1, ]$y
	yCT = Xymatched[Xymatched$indic_T == 0, ]$y
	matched_correlation = cor(yMT, yCT)	
}

## now we have to monte-carlo the exact test
Xyleft_copy = Xyleft
b_T_sims = array(NA, Nsim_exact_test)
for (nsim_exact_test in 1 : Nsim_exact_test){
	
	#### unconditional permuting - flip coins for each (equivalent to conditional permuting)
	trt_vec_multiple = (rbinom(max(Xy$match_indic), 1, prob_trt) - 0.5) * 2
	ydiffs_copy = ydiffs * trt_vec_multiple
	d_bar_samp = mean(ydiffs_copy)
	ssqD_bar_samp = var(ydiffs_copy) / length(ydiffs_copy)
	
	### conditional permuting for the reservoir
	Xyleft_copy$indic_T = sample(c(rep(1, nRT), rep(0, nRC)))	
	YleftT_samp = Xyleft_copy[Xyleft_copy$indic_T == 1, ]$y
	YleftC_samp = Xyleft_copy[Xyleft_copy$indic_T == 0, ]$y
	r_bar_samp = mean(YleftT_samp) - mean(YleftC_samp)
	ssqR_samp = (var(YleftT_samp) * (nRT - 1) + var(YleftC_samp) * (nRC - 1)) / (nR - 2) * (1 / nRT + 1 / nRC)	
	
	#now compute b_T_sim from the permuted stats
	w_star_samp = ssqR_samp / (ssqR_samp + ssqD_bar_samp)
	
	if (m == 0){
		b_T_sims[nsim_exact_test] = r_bar_samp
	} else if (nRT <= 2 || nRC <= 2){
		b_T_sims[nsim_exact_test] = d_bar_samp
	} else {
		b_T_sims[nsim_exact_test] = w_star_samp * d_bar_samp + (1 - w_star_samp) * r_bar_samp
	}	
}

#this is the empirical two-sided p-value based on simulation
pval = sum(abs(b_T_sims) > abs(beta_hat_T)) / Nsim_exact_test
